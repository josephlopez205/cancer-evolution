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

## Generate a deterministic simulated sample

```powershell
uv run clonetrace simulate --seed 42 --depth 100 --purity 0.8 --clones 3 --topology branching --mutations 40
```

The command writes canonical mutation and copy-number TSV files plus `ground_truth.json` and a SHA-256 manifest under `data/simulated/`.

## Generate the benchmark scenario grid

```powershell
uv run clonetrace simulate-grid --seed 1 --replicates 20 --mutations 40
```

This generates all depth, purity, clone-count, and topology combinations specified in P2-T2,
with 20 seeded replicates per scenario. The grid index and aggregate manifest are written under
`data/simulated/grid/`.

## Preview the Phase 2 simulation viewer

The early scientific viewer displays a bundled, simulator-generated sample with its clone tree,
mutation read counts, and copy-number events. The values are ground truth for a simulated sample,
not inferred clinical results.

```powershell
cd viewer
npm install
npm run dev
```

See [viewer/README.md](viewer/README.md) for validation commands and the current data limitations.

## Bayesian binomial clustering baseline

The Python API `clonetrace.cluster.fit_binomial_mixture` clusters read counts from one sample.
It uses Beta priors for component VAFs, a symmetric Dirichlet prior for mixture weights, and
variational posterior assignment probabilities. Candidate cluster counts are compared using the
variational evidence lower bound. Component means are observed VAFs; copy-number and purity
adjustments remain separate QC inputs; see the QC and clustering CLI commands below.

## QC and clustering CLI

The QC command reads canonical mutation columns plus per-mutation `purity`, `total_cn`, and
`multiplicity` annotations. It writes observed VAF and estimated CCF without clipping estimates
above 1:

```powershell
clonetrace qc --input mutations_with_qc.tsv --output qc.tsv
```

Cluster one sample with the built-in baseline:

```powershell
clonetrace cluster --input mutations.tsv --output assignments.tsv --method baseline
```

The optional PyClone-VI method additionally requires `major_cn`, `minor_cn`, `normal_cn`, and
`tumour_content` columns on every mutation row. Its external input, fit, and result files are
retained beside the normalized assignments TSV. Both methods write the same assignment columns;
MOBSTER remains an external integration route without a CloneTrace CLI adapter.

## PyClone-VI integration

The optional `clonetrace.cluster.run_pyclonevi` wrapper invokes the external `pyclone-vi`
executable in two stages and normalizes its TSV output into CloneTrace mutation assignments.
Install PyClone-VI separately using its
[upstream installation instructions](https://github.com/Roth-Lab/pyclone-vi#installation);
it is not a dependency of the CloneTrace Python environment.

Each input mutation must be paired with explicit major/minor tumor copy number, normal copy
number, and tumor content via `PyCloneObservation`. The wrapper does not guess segment overlap,
sex-chromosome copy number, or purity defaults. It retains the prepared input, HDF5 fit, and
normalized-tool TSV in the requested output directory. The external tool's own license and
installation requirements apply.

A seeded MOBSTER smoke fit has completed on simulated data, but MOBSTER is still an optional
external R integration rather than a runnable CloneTrace wrapper. See
[docs/mobster_integration.md](docs/mobster_integration.md) for the run record and adapter limits.
