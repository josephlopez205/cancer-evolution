from typer.testing import CliRunner

from clonetrace import __version__
from clonetrace.cli import app

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
