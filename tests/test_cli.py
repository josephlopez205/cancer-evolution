import csv

import pytest
from typer.testing import CliRunner

from clonetrace import __version__
from clonetrace.cli import app
from clonetrace.cluster import (
    MutationClusterAssignment,
    PyCloneVIClusterSummary,
    PyCloneVIResult,
)
from clonetrace.simulate.simulator import simulate

runner = CliRunner()


def test_help() -> None:
    result = runner.invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "Infer and analyze tumor clonal architecture" in result.stdout


def test_version() -> None:
    result = runner.invoke(app, ["--version"])

    assert result.exit_code == 0
    assert result.stdout.strip() == f"clonetrace {__version__}"


def test_ingest_rejects_unknown_source() -> None:
    result = runner.invoke(
        app,
        ["ingest", "--source", "tcga", "--project", "TCGA-BRCA"],
    )

    assert result.exit_code != 0


def test_ingest_cbioportal_requires_study() -> None:
    result = runner.invoke(
        app,
        ["ingest", "--source", "cbioportal"],
    )

    assert result.exit_code != 0


def test_ingest_geo_requires_license_and_processed_dir() -> None:
    result = runner.invoke(
        app,
        ["ingest", "--source", "geo", "--accession", "GSE12345"],
    )

    assert result.exit_code != 0


def test_cohort_command_writes_table(tmp_path) -> None:
    result = runner.invoke(
        app,
        [
            "cohort",
            "--input-dir",
            str(tmp_path),
            "--output-dir",
            str(tmp_path),
            "--source",
            "gdc",
        ],
    )

    assert result.exit_code == 0
    assert (tmp_path / "cohort.tsv").exists()


def test_simulate_command_writes_ground_truth(tmp_path) -> None:
    result = runner.invoke(
        app,
        [
            "simulate",
            "--output-dir",
            str(tmp_path),
            "--seed",
            "11",
            "--depth",
            "50",
            "--purity",
            "0.5",
            "--clones",
            "3",
            "--topology",
            "linear",
            "--mutations",
            "15",
        ],
    )

    assert result.exit_code == 0
    assert (tmp_path / "mutations.tsv").exists()
    assert (tmp_path / "copy_number_segments.tsv").exists()
    assert (tmp_path / "ground_truth.json").exists()
    assert (tmp_path / "manifest.json").exists()


def test_simulate_grid_command_writes_index(tmp_path) -> None:
    result = runner.invoke(
        app,
        [
            "simulate-grid",
            "--output-dir",
            str(tmp_path),
            "--seed",
            "5",
            "--replicates",
            "1",
            "--mutations",
            "10",
        ],
    )

    assert result.exit_code == 0
    assert (tmp_path / "scenario_index.tsv").exists()
    assert (tmp_path / "scenario_manifest.json").exists()


def test_qc_and_cluster_commands_smoke_on_simulated_sample(tmp_path) -> None:
    simulation_dir = tmp_path / "simulation"
    simulated = simulate(
        simulation_dir,
        seed=17,
        depth=200,
        purity=0.8,
        clone_count=3,
        mutation_count=20,
        cn_event_probability=0,
    )
    qc_input = tmp_path / "mutations_with_qc.tsv"
    with simulated.mutations_path.open("r", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source, delimiter="\t")
        fieldnames = [*reader.fieldnames, "purity", "total_cn", "multiplicity"]
        rows = [
            {**row, "purity": "0.8", "total_cn": "2", "multiplicity": "1"}
            for row in reader
        ]
    with qc_input.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    qc_output = tmp_path / "qc.tsv"
    qc_result = runner.invoke(
        app,
        ["qc", "--input", str(qc_input), "--output", str(qc_output)],
    )
    assert qc_result.exit_code == 0, qc_result.stdout
    with qc_output.open("r", encoding="utf-8", newline="") as handle:
        qc_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(qc_rows) == 20
    assert all(float(row["estimated_ccf"]) > 0 for row in qc_rows)

    cluster_output = tmp_path / "clusters.tsv"
    cluster_result = runner.invoke(
        app,
        [
            "cluster",
            "--input",
            str(simulated.mutations_path),
            "--output",
            str(cluster_output),
            "--method",
            "baseline",
            "--max-clusters",
            "3",
        ],
    )
    assert cluster_result.exit_code == 0, cluster_result.stdout
    with cluster_output.open("r", encoding="utf-8", newline="") as handle:
        cluster_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(cluster_rows) == 20
    assert all(row["method"] == "baseline" for row in cluster_rows)
    assert all(0 <= float(row["posterior_probability"]) <= 1 for row in cluster_rows)


def test_qc_command_preserves_ccf_estimates_above_one(tmp_path) -> None:
    input_path = tmp_path / "mutation.tsv"
    input_path.write_text(
        "sample\tgene\tchrom\tpos\tref\talt\tref_count\talt_count\t"
        "purity\ttotal_cn\tmultiplicity\n"
        "sample-1\tTP53\t1\t10\tA\tC\t40\t60\t0.8\t2\t1\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "qc.tsv"

    result = runner.invoke(
        app,
        ["qc", "--input", str(input_path), "--output", str(output_path)],
    )

    assert result.exit_code == 0, result.stdout
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert float(row["estimated_ccf"]) == pytest.approx(1.5)
    assert row["ccf_above_one"] == "true"
    assert "none were clipped" in result.stdout


def test_cluster_commands_use_the_same_output_schema_for_each_method(
    tmp_path,
    monkeypatch,
) -> None:
    input_path = tmp_path / "mutations.tsv"
    input_path.write_text(
        "sample\tgene\tchrom\tpos\tref\talt\tref_count\talt_count\t"
        "major_cn\tminor_cn\tnormal_cn\ttumour_content\n"
        "sample-1\tTP53\t1\t10\tA\tC\t70\t30\t2\t1\t2\t0.75\n"
        "sample-1\tPIK3CA\t1\t20\tG\tT\t65\t35\t2\t1\t2\t0.75\n",
        encoding="utf-8",
    )
    baseline_output = tmp_path / "baseline.tsv"
    baseline_result = runner.invoke(
        app,
        [
            "cluster",
            "--input",
            str(input_path),
            "--output",
            str(baseline_output),
            "--method",
            "baseline",
        ],
    )
    assert baseline_result.exit_code == 0, baseline_result.stdout

    def fake_run_pyclonevi(observations, output_dir, **kwargs):
        assignments = tuple(
            MutationClusterAssignment(
                mutation=observation.mutation,
                cluster_id=0,
                posterior_probability=0.9,
                cluster_probabilities=None,
                cellular_prevalence=0.7,
                cellular_prevalence_std=0.05,
            )
            for observation in observations
        )
        return PyCloneVIResult(
            sample="sample-1",
            method="pyclonevi",
            cluster_count=1,
            assignments=assignments,
            clusters=(PyCloneVIClusterSummary(0, len(assignments), 0.7),),
            input_path=output_dir / "input.tsv",
            fit_path=output_dir / "fit.h5",
            results_path=output_dir / "results.tsv",
        )

    monkeypatch.setattr("clonetrace.cli.run_pyclonevi", fake_run_pyclonevi)
    pyclone_output = tmp_path / "pyclonevi.tsv"
    pyclone_result = runner.invoke(
        app,
        [
            "cluster",
            "--input",
            str(input_path),
            "--output",
            str(pyclone_output),
            "--method",
            "pyclonevi",
        ],
    )
    assert pyclone_result.exit_code == 0, pyclone_result.stdout
    with baseline_output.open("r", encoding="utf-8", newline="") as handle:
        baseline_headers = next(csv.reader(handle, delimiter="\t"))
    with pyclone_output.open("r", encoding="utf-8", newline="") as handle:
        pyclone_headers = next(csv.reader(handle, delimiter="\t"))
    assert pyclone_headers == baseline_headers
    assert len(pyclone_headers) == len(set(pyclone_headers))
    with pyclone_output.open("r", encoding="utf-8", newline="") as handle:
        pyclone_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(pyclone_rows) == 2
    assert pyclone_rows[0]["cellular_prevalence"] == "0.7"


def test_cluster_pyclonevi_requires_explicit_annotations(tmp_path) -> None:
    input_path = tmp_path / "mutations.tsv"
    input_path.write_text(
        "sample\tgene\tchrom\tpos\tref\talt\tref_count\talt_count\n"
        "sample-1\tTP53\t1\t10\tA\tC\t70\t30\n",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "cluster",
            "--input",
            str(input_path),
            "--output",
            str(tmp_path / "clusters.tsv"),
            "--method",
            "pyclonevi",
        ],
    )

    assert result.exit_code != 0
    assert "major_cn, minor_cn, normal_cn, tumour_content" in result.output
