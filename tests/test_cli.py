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
