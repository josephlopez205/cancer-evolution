"""CloneTrace command-line interface."""

import typer

from clonetrace import __version__
from clonetrace.ingest.gdc import ingest_gdc

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
    project: str = typer.Option(..., "--project", help="GDC project ID, for example TCGA-BRCA"),
    cache_dir: str = typer.Option("data/gdc", "--cache-dir", help="GDC download cache"),
    output_dir: str = typer.Option(
        "data/processed",
        "--output-dir",
        help="Canonical output directory",
    ),
    force: bool = typer.Option(False, "--force", help="Redownload cached files"),
) -> None:
    """Download and harmonize data from an approved open source."""

    if source.lower() != "gdc":
        raise typer.BadParameter("only --source gdc is supported")
    result = ingest_gdc(
        project=project,
        cache_dir=cache_dir,
        output_dir=output_dir,
        force=force,
    )
    typer.echo(
        f"Ingested {result.project}: "
        f"{result.mutation_count} mutations, "
        f"{result.segment_record_count} copy-number segments"
    )
    typer.echo(f"Manifest: {result.manifest}")
