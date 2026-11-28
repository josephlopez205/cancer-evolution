"""CloneTrace command-line interface."""

from typing import Literal, cast

import typer

from clonetrace import __version__
from clonetrace.ingest.cbioportal import ingest_cbioportal
from clonetrace.ingest.cohort import build_cohort_table
from clonetrace.ingest.gdc import ingest_gdc
from clonetrace.ingest.geo import ingest_geo
from clonetrace.simulate import simulate

app = typer.Typer(
    name="clonetrace",
    help="Infer and analyze tumor clonal architecture.",
    no_args_is_help=True,
)


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
