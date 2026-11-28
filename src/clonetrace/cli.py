"""CloneTrace command-line interface."""

import csv
import json
import math
from collections.abc import Sequence
from pathlib import Path
from typing import Annotated, Literal, cast

import typer
from pydantic import ValidationError

from clonetrace import __version__
from clonetrace.cluster import (
    MutationClusterAssignment,
    PyCloneObservation,
    PyCloneVIError,
    fit_binomial_mixture,
    run_pyclonevi,
)
from clonetrace.ingest.cbioportal import ingest_cbioportal
from clonetrace.ingest.cohort import build_cohort_table
from clonetrace.ingest.gdc import ingest_gdc
from clonetrace.ingest.geo import ingest_geo
from clonetrace.qc import estimate_ccf_from_vaf
from clonetrace.schemas import MutationRecord
from clonetrace.simulate import generate_scenario_grid, simulate

app = typer.Typer(
    name="clonetrace",
    help="Infer and analyze tumor clonal architecture.",
    no_args_is_help=True,
)

_MUTATION_COLUMNS = (
    "sample",
    "gene",
    "chrom",
    "pos",
    "ref",
    "alt",
    "ref_count",
    "alt_count",
)
_QC_ANNOTATION_COLUMNS = {"purity", "total_cn", "multiplicity"}
_PYCLONEVI_ANNOTATION_COLUMNS = {
    "major_cn",
    "minor_cn",
    "normal_cn",
    "tumour_content",
}
_CLUSTER_OUTPUT_COLUMNS = (
    "mutation_id",
    "sample",
    "gene",
    "chrom",
    "pos",
    "ref",
    "alt",
    "ref_count",
    "alt_count",
    "method",
    "cluster_id",
    "posterior_probability",
    "cluster_probabilities",
    "cellular_prevalence",
    "cellular_prevalence_std",
)


def _read_mutation_rows(
    path: Path,
    *,
    required_annotations: set[str] | None = None,
) -> list[tuple[MutationRecord, dict[str, str | None], int]]:
    allowed_columns = (
        set(_MUTATION_COLUMNS) | _QC_ANNOTATION_COLUMNS | _PYCLONEVI_ANNOTATION_COLUMNS
    )
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = reader.fieldnames
        if not headers:
            raise typer.BadParameter(f"{path}: TSV must include a header")
        if any(not header for header in headers) or len(headers) != len(set(headers)):
            raise typer.BadParameter(f"{path}: TSV has empty or duplicate column names")
        missing = set(_MUTATION_COLUMNS) - set(headers)
        if missing:
            raise typer.BadParameter(
                f"{path}: TSV is missing required columns: {', '.join(sorted(missing))}"
            )
        unknown = set(headers) - allowed_columns
        if unknown:
            raise typer.BadParameter(
                f"{path}: TSV has unsupported columns: {', '.join(sorted(unknown))}"
            )
        missing_annotations = (required_annotations or set()) - set(headers)
        if missing_annotations:
            raise typer.BadParameter(
                f"{path}: TSV is missing required annotation columns: "
                f"{', '.join(sorted(missing_annotations))}"
            )

        rows: list[tuple[MutationRecord, dict[str, str | None], int]] = []
        seen_mutations: set[tuple[str, str, int, str, str]] = set()
        for line, row in enumerate(reader, start=2):
            if None in row:
                raise typer.BadParameter(f"{path}: line {line} has extra values")
            if all(value is None or not value.strip() for value in row.values()):
                continue
            try:
                mutation = MutationRecord.model_validate(
                    {column: row.get(column) for column in _MUTATION_COLUMNS}
                )
            except ValidationError as error:
                raise typer.BadParameter(
                    f"{path}: line {line} has an invalid mutation record: {error}"
                ) from error
            identity = (mutation.sample, mutation.chrom, mutation.pos, mutation.ref, mutation.alt)
            if identity in seen_mutations:
                raise typer.BadParameter(
                    f"{path}: line {line} duplicates mutation "
                    f"{mutation.sample}:{mutation.chrom}:{mutation.pos}:{mutation.ref}>{mutation.alt}"
                )
            seen_mutations.add(identity)
            rows.append((mutation, row, line))

    if not rows:
        raise typer.BadParameter(f"{path}: TSV contains no mutation records")
    return rows


def _mutation_id(mutation: MutationRecord) -> str:
    return f"{mutation.sample}|{mutation.chrom}|{mutation.pos}|{mutation.ref}>{mutation.alt}"


def _parse_annotation_number(
    row: dict[str, str | None],
    column: str,
    *,
    path: Path,
    line: int,
) -> float:
    value = row.get(column)
    if value is None or not value.strip():
        raise typer.BadParameter(f"{path}: line {line}: {column} must be a number")
    try:
        parsed = float(value)
    except ValueError as error:
        raise typer.BadParameter(
            f"{path}: line {line}: {column} must be a number"
        ) from error
    if not math.isfinite(parsed):
        raise typer.BadParameter(f"{path}: line {line}: {column} must be finite")
    return parsed


def _parse_multiplicity(row: dict[str, str | None], *, path: Path, line: int) -> int:
    value = _parse_annotation_number(row, "multiplicity", path=path, line=line)
    if not value.is_integer():
        raise typer.BadParameter(f"{path}: line {line}: multiplicity must be an integer")
    return int(value)


def _write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _cluster_output_rows(
    assignments: Sequence[MutationClusterAssignment],
    *,
    method: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for assignment in assignments:
        mutation = assignment.mutation
        rows.append(
            {
                "mutation_id": _mutation_id(mutation),
                "sample": mutation.sample,
                "gene": mutation.gene,
                "chrom": mutation.chrom,
                "pos": mutation.pos,
                "ref": mutation.ref,
                "alt": mutation.alt,
                "ref_count": mutation.ref_count,
                "alt_count": mutation.alt_count,
                "method": method,
                "cluster_id": assignment.cluster_id,
                "posterior_probability": assignment.posterior_probability,
                "cluster_probabilities": (
                    json.dumps(assignment.cluster_probabilities, separators=(",", ":"))
                    if assignment.cluster_probabilities is not None
                    else ""
                ),
                "cellular_prevalence": (
                    assignment.cellular_prevalence
                    if assignment.cellular_prevalence is not None
                    else ""
                ),
                "cellular_prevalence_std": (
                    assignment.cellular_prevalence_std
                    if assignment.cellular_prevalence_std is not None
                    else ""
                ),
            }
        )
    return rows


@app.callback(invoke_without_command=True)
def main(
    version: bool = typer.Option(
        False,
        "--version",
        "-v",
        help="Show the CloneTrace version.",
    ),
) -> None:
    """Run CloneTrace."""
    if version:
        typer.echo(f"clonetrace {__version__}")


@app.command("ingest")
def ingest(
    source: str = typer.Option(..., "--source", help="Data source: gdc"),
    project: str = typer.Option("", "--project", help="GDC project ID, for example TCGA-BRCA"),
    study: str = typer.Option("", "--study", help="cBioPortal study ID"),
    accession: str = typer.Option("", "--accession", help="GEO accession, for example GSE12345"),
    processed_dir: str = typer.Option(
        "", "--processed-dir", help="Directory of processed GEO files"
    ),
    license_url: str = typer.Option(
        "", "--license-url", help="Verified processed-data license URL"
    ),
    cache_dir: str = typer.Option(
        "data/cache",
        "--cache-dir",
        help="Response cache for supported sources; GDC source files are never cached",
    ),
    output_dir: str = typer.Option(
        "data/processed",
        "--output-dir",
        help="Canonical output directory",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Refresh cached response data for cBioPortal or GEO",
    ),
) -> None:
    """Download and harmonize data from an approved open source."""

    source_name = source.lower()
    if source_name == "geo":
        if not accession:
            raise typer.BadParameter("--accession is required for --source geo")
        if not processed_dir:
            raise typer.BadParameter("--processed-dir is required for --source geo")
        if not license_url:
            raise typer.BadParameter("--license-url is required for --source geo")
        geo_result = ingest_geo(
            accession=accession,
            processed_dir=processed_dir,
            license_url=license_url,
            cache_dir=cache_dir,
            output_dir=output_dir,
        )
        typer.echo(
            f"Ingested {geo_result.accession}: "
            f"{geo_result.mutation_count} mutations, "
            f"{geo_result.segment_count} copy-number segments, "
            f"{geo_result.clinical_count} clinical records"
        )
        typer.echo(f"Manifest: {geo_result.manifest}")
        return
    if source_name == "cbioportal":
        if not study:
            raise typer.BadParameter("--study is required for --source cbioportal")
        cbio_result = ingest_cbioportal(
            study=study,
            cache_dir=cache_dir,
            output_dir=output_dir,
            force=force,
        )
        typer.echo(
            f"Ingested {cbio_result.study}: "
            f"{cbio_result.mutation_count} mutations, "
            f"{cbio_result.clinical_count} clinical records, "
            f"{cbio_result.cna_count} CNA calls"
        )
        typer.echo(f"Manifest: {cbio_result.manifest}")
        return
    if source_name == "gdc":
        if not project:
            raise typer.BadParameter("--project is required for --source gdc")
        result = ingest_gdc(
            project=project,
            output_dir=output_dir,
        )
        typer.echo(
            f"Ingested {result.project}: "
            f"{result.mutation_count} mutations, "
            f"{result.segment_record_count} copy-number segments, "
            f"{result.clinical_count} clinical records"
        )
        typer.echo(f"Cohort table: {result.output_dir / 'cohort.tsv'}")
        typer.echo(f"Manifest: {result.manifest}")
        return
    raise typer.BadParameter("supported sources are gdc, cbioportal, and geo")


@app.command("cohort")
def cohort(
    input_dir: str = typer.Option("data/processed", "--input-dir"),
    output_dir: str = typer.Option("data/processed", "--output-dir"),
    source: str = typer.Option("combined", "--source"),
) -> None:
    """Build a harmonized cohort table from canonical ingested TSV files."""

    path = build_cohort_table(input_dir=input_dir, output_dir=output_dir, source=source)
    typer.echo(f"Cohort table: {path}")


@app.command("simulate")
def simulate_command(
    output_dir: str = typer.Option("data/simulated", "--output-dir"),
    seed: int = typer.Option(1, "--seed"),
    depth: int = typer.Option(100, "--depth"),
    purity: float = typer.Option(0.8, "--purity"),
    clone_count: int = typer.Option(3, "--clones"),
    topology: str = typer.Option("branching", "--topology"),
    mutation_count: int = typer.Option(40, "--mutations"),
    cn_event_probability: float = typer.Option(0.25, "--cn-event-probability"),
) -> None:
    """Generate a deterministic simulated sample and ground truth."""

    if topology not in {"linear", "branching"}:
        raise typer.BadParameter("--topology must be linear or branching")
    result = simulate(
        output_dir=output_dir,
        seed=seed,
        depth=depth,
        purity=purity,
        clone_count=clone_count,
        topology=cast(Literal["linear", "branching"], topology),
        mutation_count=mutation_count,
        cn_event_probability=cn_event_probability,
    )
    typer.echo(
        f"Simulated {result.mutation_count} mutations across "
        f"{result.clone_count} clones"
    )
    typer.echo(f"Ground truth: {result.ground_truth_path}")


@app.command("simulate-grid")
def simulate_grid_command(
    output_dir: str = typer.Option("data/simulated/grid", "--output-dir"),
    seed: int = typer.Option(1, "--seed"),
    replicates: int = typer.Option(20, "--replicates"),
    mutation_count: int = typer.Option(40, "--mutations"),
) -> None:
    """Generate the reproducible benchmark scenario grid."""

    result = generate_scenario_grid(
        output_dir=output_dir,
        seed=seed,
        replicates=replicates,
        mutation_count=mutation_count,
    )
    typer.echo(
        f"Generated {result.scenario_count * result.replicate_count} simulations "
        f"from {result.scenario_count} scenarios"
    )
    typer.echo(f"Scenario index: {result.index_path}")
    typer.echo(f"Scenario manifest: {result.manifest_path}")


@app.command("qc")
def qc_command(
    input_path: Annotated[
        Path,
        typer.Option("--input", help="Mutation TSV with per-row QC annotations"),
    ],
    output_path: Annotated[Path, typer.Option("--output", help="Output QC TSV")],
) -> None:
    """Estimate mutation CCF from VAF, purity, copy number, and multiplicity."""

    if input_path.resolve() == output_path.resolve():
        raise typer.BadParameter("--output must differ from --input")
    rows = _read_mutation_rows(
        input_path,
        required_annotations=_QC_ANNOTATION_COLUMNS,
    )
    output_rows: list[dict[str, object]] = []
    above_one = 0
    for mutation, row, line in rows:
        depth = mutation.ref_count + mutation.alt_count
        if depth == 0:
            raise typer.BadParameter(
                f"{input_path}: line {line}: total read depth must be positive"
            )
        purity = _parse_annotation_number(row, "purity", path=input_path, line=line)
        total_cn = _parse_annotation_number(row, "total_cn", path=input_path, line=line)
        multiplicity = _parse_multiplicity(row, path=input_path, line=line)
        try:
            ccf = estimate_ccf_from_vaf(
                vaf=mutation.alt_count / depth,
                purity=purity,
                total_cn=total_cn,
                multiplicity=multiplicity,
            )
        except ValueError as error:
            raise typer.BadParameter(f"{input_path}: line {line}: {error}") from error
        above_one += ccf > 1
        output_rows.append(
            {
                "mutation_id": _mutation_id(mutation),
                **mutation.model_dump(mode="python"),
                "purity": purity,
                "total_cn": total_cn,
                "multiplicity": multiplicity,
                "observed_vaf": mutation.alt_count / depth,
                "estimated_ccf": ccf,
                "ccf_above_one": str(ccf > 1).lower(),
            }
        )

    _write_tsv(
        output_path,
        (
            "mutation_id",
            *_MUTATION_COLUMNS,
            "purity",
            "total_cn",
            "multiplicity",
            "observed_vaf",
            "estimated_ccf",
            "ccf_above_one",
        ),
        output_rows,
    )
    typer.echo(f"Wrote {len(output_rows)} QC records to {output_path}")
    if above_one:
        typer.echo(
            f"Review required: {above_one} estimated CCF value(s) exceed 1; none were clipped."
        )


@app.command("cluster")
def cluster_command(
    input_path: Annotated[Path, typer.Option("--input", help="Single-sample mutation TSV")],
    output_path: Annotated[
        Path,
        typer.Option("--output", help="Standardized cluster assignment TSV"),
    ],
    method: Annotated[
        str,
        typer.Option("--method", help="Clustering method: baseline or pyclonevi"),
    ],
    max_clusters: int = typer.Option(6, "--max-clusters", min=1),
    pyclonevi_executable: str = typer.Option("pyclone-vi", "--pyclonevi-executable"),
    num_clusters: int = typer.Option(10, "--num-clusters", min=2),
    density: str = typer.Option("binomial", "--density"),
    seed: int = typer.Option(1, "--seed", min=0),
) -> None:
    """Cluster one sample and write normalized mutation assignments."""

    method = method.lower()
    if method not in {"baseline", "pyclonevi"}:
        raise typer.BadParameter(
            "supported methods are baseline and pyclonevi; MOBSTER has no CloneTrace adapter"
        )
    if input_path.resolve() == output_path.resolve():
        raise typer.BadParameter("--output must differ from --input")
    required_annotations = (
        _PYCLONEVI_ANNOTATION_COLUMNS if method == "pyclonevi" else set()
    )
    rows = _read_mutation_rows(input_path, required_annotations=required_annotations)
    mutations = [mutation for mutation, _row, _line in rows]
    samples = {mutation.sample for mutation in mutations}
    if len(samples) != 1:
        raise typer.BadParameter("cluster accepts one sample per run")

    if method == "baseline":
        try:
            baseline_result = fit_binomial_mixture(mutations, max_clusters=max_clusters)
        except ValueError as error:
            raise typer.BadParameter(str(error)) from error
        assignments = baseline_result.assignments
        cluster_count = baseline_result.cluster_count
    else:
        observations: list[PyCloneObservation] = []
        for mutation, row, line in rows:
            try:
                observations.append(
                    PyCloneObservation(
                        mutation=mutation,
                        major_cn=_parse_annotation_number(
                            row, "major_cn", path=input_path, line=line
                        ),
                        minor_cn=_parse_annotation_number(
                            row, "minor_cn", path=input_path, line=line
                        ),
                        normal_cn=_parse_annotation_number(
                            row, "normal_cn", path=input_path, line=line
                        ),
                        tumour_content=_parse_annotation_number(
                            row, "tumour_content", path=input_path, line=line
                        ),
                    )
                )
            except ValueError as error:
                raise typer.BadParameter(f"{input_path}: line {line}: {error}") from error
        run_dir = output_path.parent / f"{output_path.stem}_pyclonevi"
        try:
            pyclonevi_result = run_pyclonevi(
                observations,
                run_dir,
                executable=pyclonevi_executable,
                num_clusters=num_clusters,
                density=cast(Literal["binomial", "beta-binomial"], density),
                seed=seed,
            )
        except (PyCloneVIError, ValueError) as error:
            raise typer.BadParameter(str(error)) from error
        assignments = pyclonevi_result.assignments
        cluster_count = pyclonevi_result.cluster_count

    output_rows = _cluster_output_rows(assignments, method=method)
    _write_tsv(output_path, _CLUSTER_OUTPUT_COLUMNS, output_rows)
    typer.echo(
        f"Wrote {len(output_rows)} {method} assignments across "
        f"{cluster_count} clusters to {output_path}"
    )
