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

GDC open-tier file contents are fetched through the public API and processed in memory; raw files
are not cached or persisted locally. Canonical mutation, segment, clinical, and cohort tables are
written under `data/processed/`, and source file checksums are recorded in the manifest. Clinical
fields are fetched from the public GDC Cases API.

## Ingest a public cBioPortal study

```powershell
uv run clonetrace ingest --source cbioportal --study brca_metabric
```

This caches the public study metadata, mutation, discrete CNA, and clinical responses and writes
canonical mutation and clinical TSV tables plus a copy-number-call TSV.

## Load a processed GEO dataset

```powershell
uv run clonetrace ingest --source geo --accession GSE12345 `
  --processed-dir data/external/gse12345 `
  --license-url https://example.org/verified-license
```

The GEO loader only accepts processed MAF, segment, and clinical TSV files. It rejects raw BAM,
CRAM, and SRA files and requires an absolute URL for the dataset's verified license.

## Build a cohort table

```powershell
uv run clonetrace cohort --input-dir data/processed --output-dir data/processed --source gdc
```
