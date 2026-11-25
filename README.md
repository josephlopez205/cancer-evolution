# CloneTrace

CloneTrace is a reproducible, open-source Python pipeline for inferring clonal architecture and evolutionary trajectories from public tumor sequencing data.

The project is in early development. Source selection and the implementation roadmap are maintained in [PLAN.md](PLAN.md).

## Development

CloneTrace requires Python 3.11 or newer.

```powershell
python -m pip install uv
uv sync --locked --extra dev
uv run pytest
uv run ruff check .
uv run pre-commit run --all-files
```

The command-line entry point can be checked with:

```powershell
uv run clonetrace --help
```

## Ingest open GDC data

```powershell
uv run clonetrace ingest --source gdc --project TCGA-BRCA
```

This caches open-tier MAF, copy-number, and clinical files under `data/gdc/`, writes canonical
mutation and segment tables under `data/processed/`, and verifies file size and MD5 checksums.

## Ingest a public cBioPortal study

```powershell
uv run clonetrace ingest --source cbioportal --study brca_metabric
```

This caches the public study metadata, mutation, discrete CNA, and clinical responses and writes
canonical mutation and clinical TSV tables plus a copy-number-call TSV.
