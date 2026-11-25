# PLAN.md: CloneTrace

> A reproducible, open-source pipeline for inferring clonal architecture and evolutionary trajectories from public tumor sequencing data, plus the research paper that benchmarks it.

This file is the source of truth for the project. Agents and humans should read it before starting work, tick off tasks as they are completed, and log deviations in the **Decision Log** at the bottom.

---

## 1. Project summary

| Item | Value |
|---|---|
| Working title | CloneTrace: consensus clonal inference and driver timing from open-access tumor data |
| Type | Research paper + software package |
| Language | Python 3.11+ |
| License | MIT (code), CC-BY-4.0 (docs and derived tables) |
| Duration | 16 weeks |
| Data policy | **Open data only.** No controlled-access (dbGaP, EGA) data. No COSMIC bulk downloads. |
| Primary deliverables | `clonetrace` pip package and CLI, workflow, Docker image, React web viewer, preprint |

### Research question
Can a lightweight, fully open pipeline recover subclonal structure and evolutionary timing (clonal vs. subclonal drivers, branching vs. linear evolution) from public tumor data, and how does it compare with established tools?

### Hypotheses
1. **H1:** Driver mutations within a cancer type differ in timing (early clonal vs. late subclonal).
2. **H2:** Subclonal diversity (clone count, diversity index) is associated with outcome or treatment response.
3. **H3:** A consensus across several clustering methods is more stable than any single method.

### Non-goals
- Calling variants from raw reads (BAMs are controlled-access; start from public MAF and segment files).
- Single-cell clonal inference (future work).
- Clinical or diagnostic claims.

---

## 2. Working rules for agents

1. **Read this file first**, then work on the lowest-numbered unchecked task whose dependencies are complete.
2. **One task, one branch, one PR.** Branch name: `feat/<task-id>-<slug>`.
3. **Plan before coding.** For each task, produce an implementation plan listing files to touch, tests to add, and acceptance checks, then wait for approval if the task is marked `[REVIEW]`.
4. **Tests are required.** No module is done without unit tests; no pipeline stage is done without a test on simulated data.
5. **Never fabricate data or results.** If a download fails or a source is unavailable, record it in the Decision Log and stop; do not substitute synthetic data unless the task says so.
6. **Verify licenses and access terms** before ingesting any new data source. Record the license, access date, and URL in `data/SOURCES.md`.
7. **Pin everything.** Use a lockfile; record tool versions and random seeds in every output.
8. **Report uncertainty.** Every inferred clone, tree, and timing estimate must carry a confidence or posterior value.
9. **Do not commit large data.** Raw and processed data live in `data/` (git-ignored) and are re-creatable from scripts.
10. **Keep outputs deterministic** where possible (fixed seeds, sorted outputs).

Suggested agent artifacts: an implementation plan per task, a task list mirroring Section 6, and a walkthrough (screenshots or logs) attached to each PR showing that acceptance criteria pass.

---

## 3. Tech stack

| Layer | Choice |
|---|---|
| Packaging | `uv` or `poetry`, `pyproject.toml`, `src/` layout |
| CLI | `typer` |
| Data | `pandas` or `polars`, `pyarrow`/Parquet, `pydantic` schemas |
| Clustering | PyClone-VI, MOBSTER (via wrapper), baseline Bayesian mixture (own implementation) |
| Tree inference | Pairtree (wrapper) and a minimal custom method as fallback |
| Stats | `scipy`, `statsmodels`, `lifelines`, `scikit-learn` |
| Workflow | Snakemake |
| Viewer data contract | Versioned JSON bundle exported by `clonetrace export`, validated by a JSON Schema generated from the `pydantic` models |
| Viewer (React) | React 18, TypeScript, Vite, Tailwind CSS, `d3-hierarchy` and `d3-scale` for clone trees, Plotly (`react-plotly.js`) or Recharts for charts, TanStack Table for cohort tables |
| Optional API | FastAPI (`clonetrace serve`) for local use with large cohorts; the viewer works without it from static JSON |
| Quality | `pytest`, `hypothesis`, `ruff`, `mypy`, `pre-commit`; frontend: `vitest`, React Testing Library, Playwright, `eslint`, `prettier` |
| CI | GitHub Actions |
| Docs | `mkdocs-material` |
| Containers | Docker (multi-stage), published to GHCR |

---

## 4. Repository layout

```
clonetrace/
├── PLAN.md
├── README.md
├── LICENSE
├── pyproject.toml
├── Snakefile
├── Dockerfile
├── src/clonetrace/
│   ├── cli.py
│   ├── schemas.py
│   ├── ingest/        # gdc.py, cbioportal.py, geo.py
│   ├── qc/            # purity.py, ploidy.py, vaf_cn.py
│   ├── cluster/       # pyclonevi.py, mobster.py, baseline.py
│   ├── consensus/     # coassociation.py, stability.py
│   ├── tree/          # pairtree.py, minimal.py
│   ├── timing/        # driver_class.py, mutation_time.py
│   ├── simulate/      # simulator.py, scenarios.py
│   ├── benchmark/     # metrics.py, run.py
│   ├── report/        # html.py, templates/
│   ├── export/        # bundle.py, json_schema.py (viewer data contract)
│   └── serve/         # api.py (optional FastAPI)
├── viewer/            # React + TypeScript app (Vite)
│   ├── package.json
│   ├── src/
│   │   ├── components/   # CloneTree, ClonePrevalence, DriverTimeline, CohortTable
│   │   ├── pages/        # Cohort, Sample, Benchmark, About
│   │   ├── lib/          # data loading, schema types, color scales
│   │   └── types/        # generated from the JSON Schema
│   ├── public/sample-bundle/  # small bundled demo data
│   └── tests/
├── tests/
├── workflows/         # Snakemake rules and configs
├── notebooks/         # tutorial and figure notebooks
├── docs/
├── data/              # git-ignored; SOURCES.md is tracked
├── paper/             # manuscript, figures, bibliography
└── .github/workflows/
```

---

## 5. Open data sources

| ID | Source | Use | Notes |
|---|---|---|---|
| D1 | GDC / TCGA (open tier) | Masked somatic MAFs, copy-number segments, clinical data | Open without login. Do not touch BAMs. |
| D2 | cBioPortal public studies | Curated mutation, CNA, clinical data via API | Check study-level citation requirements. |
| D3 | PCAWG (open tier) | Consensus calls, published timing for comparison | Open tier only. |
| D4 | GEO / SRA / Zenodo | Published multi-region or longitudinal studies | Verify processed calls are public and licensed. |
| D5 | DepMap / CCLE | Cell-line context for driver validation | Open. |
| D6 | Simulated data | Ground-truth benchmarking | Generated in-repo. |

Access terms change. Re-verify each source's license before ingestion and record it in `data/SOURCES.md`.

---

## 6. Phases and tasks

Legend: `[ ]` not started, `[x]` done, `[REVIEW]` needs human approval of the plan before implementation.

### Phase 0: Scoping (Weeks 1–2)

- [x] **P0-T1** Literature and tool landscape review (PyClone-VI, MOBSTER, Pairtree, CONIPHER, others). Output: `paper/notes/landscape.md` with a comparison table.
- [x] **P0-T2** `[REVIEW]` Select one primary cancer type and one validation cohort. Output: rationale in the Decision Log.
- [x] **P0-T3** Data inventory and license check for D1–D5. Output: `data/SOURCES.md`.
- [x] **P0-T4** Scaffold repository: `pyproject.toml`, `src/` layout, pre-commit, CI skeleton, README stub.

**Done when:** cancer type chosen, all sources licensed and listed, CI runs lint and an empty test suite green.

### Phase 1: Data acquisition and harmonization (Weeks 3–4)

- [x] **P1-T1** Define canonical schemas in `schemas.py`: mutations (sample, gene, chrom, pos, ref, alt, ref_count, alt_count), copy number (segments, total/minor CN), clinical (sample, survival, therapy).
- [x] **P1-T2** `ingest/gdc.py`: access open-tier MAF, allele-specific CN segments, and clinical fields through the GDC API by project ID; do not persist source files locally; checksum and produce canonical/cohort tables.
- [x] **P1-T3** `ingest/cbioportal.py`: pull study data via the public API.
- [x] **P1-T4** `ingest/geo.py`: loader for a processed multi-region or longitudinal dataset.
- [x] **P1-T5** Cohort tables and data dictionary. Output: `docs/data_dictionary.md`.

**Done when:** `clonetrace ingest --source gdc --project <ID>` reproduces mutation, copy-number, clinical, and cohort tables through the API without persisting source files locally; schema validation tests pass; no controlled-access files are referenced.

### Phase 2: Simulation benchmark (Weeks 5–6)

- [ ] **P2-T1** `simulate/simulator.py`: generate clone trees, mutation assignments, read counts with configurable depth, purity, clone count, and copy-number events; save ground truth alongside data.
- [ ] **P2-T2** Scenario grid: depth {30x, 60x, 100x, 200x}, purity {0.3, 0.5, 0.8}, clones {1–6}, topology {linear, branching}, 20 replicates each, seeded.
- [ ] **P2-T3** `benchmark/metrics.py`: clone-count error, adjusted Rand index, tree edit distance (or ancestor-descendant accuracy), calibration of reported uncertainty.

**Done when:** the simulator is deterministic by seed; metrics are tested against hand-computed toy cases.

### Phase 3: Core pipeline (Weeks 7–9)

- [ ] **P3-T1** `qc`: purity and ploidy handling (ingest estimates; optional re-estimation), VAF adjustment for copy number and multiplicity.
- [ ] **P3-T2** `cluster/baseline.py`: Bayesian binomial mixture baseline (own implementation).
- [ ] **P3-T3** `cluster/pyclonevi.py`: wrapper with standardized input and output.
- [ ] **P3-T4** `cluster/mobster.py`: wrapper, or document the integration route if R interop is required.
- [ ] **P3-T5** CLI commands `qc` and `cluster`; unit tests plus a smoke test on simulated data.

**Done when:** each clustering method runs end-to-end on one simulated dataset and one real sample through the CLI, with identical output schemas.

### Phase 4: Consensus, trees, timing (Weeks 10–11)

- [ ] **P4-T1** `consensus`: co-association matrix across methods, consensus clustering, per-mutation stability score.
- [ ] **P4-T2** `[REVIEW]` `tree`: Pairtree wrapper and a minimal fallback (sum-rule / ancestry-ordering) method.
- [ ] **P4-T3** `timing`: classify driver mutations as clonal or subclonal; estimate relative mutation timing where copy-number gain permits; flag samples where timing is not identifiable.
- [ ] **P4-T4** Snakemake workflow connecting all stages; config-driven.

**Done when:** `snakemake` runs the full pipeline on one cohort from a clean checkout; outputs include clusters, consensus, tree, and timing with uncertainty.

### Phase 5: Validation and benchmarking (Weeks 12–13)

- [ ] **P5-T1** Run all methods plus consensus on the full simulation grid; produce metric tables and figures.
- [ ] **P5-T2** Compare against published clonal calls (e.g., PCAWG or published cohort results). Report concordance and disagreements.
- [ ] **P5-T3** Sensitivity analysis to depth, purity, and exclusion thresholds; document sample exclusion criteria.
- [ ] **P5-T4** Test H3: is consensus more stable than single methods (bootstrap over mutations, rank-based tests)?

**Done when:** every figure in the Results section is regenerated by a single command, and negative results are reported alongside positive ones.

### Phase 6: Biological analysis (Week 14)

- [ ] **P6-T1** H1: driver timing across the cohort; compare known drivers (e.g., *TP53*, *PIK3CA*, *KRAS*, depending on cancer type) against expectations from the literature.
- [ ] **P6-T2** H2: associate heterogeneity metrics with survival and treatment response using Kaplan–Meier and Cox models, adjusting for stage, age, and purity; apply multiple-testing correction.
- [ ] **P6-T3** Bootstrap confidence intervals on clone calls and downstream effect sizes.

**Done when:** all analyses are scripted, seeded, and accompanied by an explicit statement of what bulk single-timepoint data can and cannot support.

### Phase 7: Viewer, docs, packaging (Week 15)

- [ ] **P7-T1** `report`: per-sample HTML report (clone fractions, tree, driver timing) and cohort summary.
- [ ] **P7-T2** `export`: `clonetrace export` writes a versioned JSON bundle (samples, clusters, consensus, trees, timing, cohort statistics, benchmark results) and a JSON Schema generated from the `pydantic` models. Include a small demo bundle in `viewer/public/sample-bundle/`.
- [ ] **P7-T3** `[REVIEW]` `viewer`: scaffold the React + TypeScript app (Vite, Tailwind, ESLint, Vitest), generate TypeScript types from the JSON Schema, and implement bundle loading with schema-version checks and clear error states.
- [ ] **P7-T4** `viewer` views:
  - **Cohort page:** sortable and filterable sample table, heterogeneity summaries, Kaplan–Meier plots by clone-count group.
  - **Sample page:** interactive clone tree (hover for mutations and driver genes), clone prevalence plot, VAF cluster scatter with consensus stability coloring, driver timing timeline.
  - **Benchmark page:** metric plots across the simulation grid (depth, purity, clone count).
  - Every view shows uncertainty (posterior or stability) and a notice when timing is not identifiable.
- [ ] **P7-T5** `serve` (optional): `clonetrace serve` launches a local FastAPI backend plus the built viewer for large cohorts; the viewer must still work without it from static JSON.
- [ ] **P7-T6** Documentation site, tutorial notebook run end-to-end in CI on a small bundled example, and a "using the viewer" guide with screenshots.
- [ ] **P7-T7** Dockerfile (Python pipeline plus built viewer), GHCR publishing, PyPI test release (include the built viewer assets in the wheel), and optional static deployment of the demo viewer to GitHub Pages.

**Done when:** a new user can install from PyPI or Docker, run `clonetrace export` then open the viewer, and reproduce the tutorial in under 30 minutes; the viewer loads both the demo bundle and a freshly exported bundle without errors.

### Phase 8: Manuscript and release (Week 16)

- [ ] **P8-T1** Draft manuscript in `paper/` (Introduction, Methods, Results, Discussion, Availability).
- [ ] **P8-T2** Final figures, supplementary tables, and a data and code availability statement.
- [ ] **P8-T3** Tag `v1.0.0`, archive on Zenodo for a DOI, post preprint to bioRxiv.
- [ ] **P8-T4** Target venue shortlist: *Bioinformatics*, *NAR Genomics and Bioinformatics*, *GigaScience*.

**Done when:** a tagged release exists with a DOI, and a stranger can reproduce the headline figure from the README.

---

## 7. Testing and verification

| Level | What | Tool |
|---|---|---|
| Unit | Schemas, VAF adjustment, metrics, consensus math | `pytest`, `hypothesis` |
| Integration | Each CLI command on a tiny bundled dataset | `pytest` with fixtures |
| End-to-end | Full Snakemake run on simulated data | CI (scheduled and on release) |
| Contract | Exported bundle validates against the generated JSON Schema; viewer types regenerate without diff | `pytest`, `jsonschema`, CI diff check |
| Frontend unit | Components (clone tree, tables, plots) render from fixture bundles | `vitest`, React Testing Library |
| Frontend E2E | Load demo bundle, open a sample, interact with the tree, check error states | Playwright |
| Regression | Pinned benchmark scores on a small grid | CI with tolerance bands |
| Reproducibility | Docker build and tutorial notebook execution | CI |

Every PR must include: passing CI, new or updated tests, and a short note describing which acceptance criteria were checked.

---

## 8. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Low coverage or poor purity estimates in public data | Filter by depth and purity thresholds; report exclusions |
| Scope creep | One cancer type and one modality first; generalize only after Phase 5 |
| Limited novelty (similar pipelines exist) | Make the contribution the consensus method, rigorous benchmarking, and usability |
| Upstream tool or API changes | Pin versions, cache downloads, add API contract tests |
| License or access changes | Re-check `data/SOURCES.md` at Phase 1 and Phase 8 |
| R or Python interop friction (MOBSTER) | Isolate in a container; allow the method to be optional |
| Over-interpreting bulk data | State identifiability limits; flag low-confidence samples |
| Python and viewer data drift | Generate TypeScript types from the JSON Schema; fail CI on schema or type diffs; version the bundle format |
| Large cohorts slow the browser viewer | Lazy-load per-sample files, paginate tables, offer the optional `serve` backend |
| Frontend scope creep | Ship three views (Cohort, Sample, Benchmark) for v1.0; defer editing, uploads, and accounts |

---

## 9. Definition of done (project level)

- [ ] All Phase 0–8 tasks are checked.
- [ ] Every analysis in the paper reruns from a clean clone with documented commands.
- [ ] Open-data-only constraint verified against `data/SOURCES.md`.
- [ ] Package installs from PyPI and Docker; tutorial passes in CI.
- [ ] Preprint posted; release archived with a DOI.

---

## 10. Open decisions

- [x] Primary cancer type: breast invasive carcinoma (TCGA-BRCA), chosen from open-data availability and cohort size.
- [x] Validation cohort: METABRIC breast cancer (cBioPortal `brca_metabric`).
- [ ] Whether MOBSTER is included in v1.0 or deferred.
- [ ] Minimum coverage and purity thresholds for sample inclusion.
- [ ] Final package name and PyPI availability.
- [ ] Charting library for the viewer (Plotly vs. Recharts) and whether clone trees use plain D3 or a graph library.
- [ ] Whether to host a public demo viewer (GitHub Pages) with the bundled demo data.

## 11. Decision log

| Date | Decision | Reason | Decided by |
|---|---|---|---|
| 2026-11-26 | Viewer is a React + TypeScript app, not Streamlit. Data flows through a versioned JSON bundle; FastAPI backend is optional. | Project requirement; a static-first design keeps the viewer deployable and reproducible without a Python server. | Project owner |
| 2026-11-26 | Selected breast invasive carcinoma (TCGA-BRCA) as the primary cancer type and METABRIC (`brca_metabric`) as the independent validation cohort. | Open GDC availability: TCGA-BRCA has 1098 cases, 992 Masked Somatic Mutation file records, 11,056 open copy-number file records, and 2,288 open clinical file records; TCGA-LUAD has 585 cases and TCGA-COAD 461. cBioPortal reports 1084 TCGA PanCancer Atlas breast samples and 2509 METABRIC breast samples, versus 566 LUAD and 594 COADREAD. Thus BRCA has the largest mutation, copy-number, clinical, and independent validation-cohort availability among the candidates. | Project owner, based on GDC and cBioPortal API queries |
| 2026-11-26 | Recorded an initial open-data inventory for D1-D5 in `data/SOURCES.md`. | D1-D2 have identified open/public primary sources; D3 remains blocked until a specific open-tier PCAWG record is verified, D4 until an exact processed study/record is selected and licensed, and D5 until a specific release and terms are checked. This prevents public portal availability from being mistaken for a verified download license. | Project owner |
| 2026-11-26 | Scaffolded the CloneTrace Python package. | P0-T4 selected the `src/` package layout, Typer CLI, Pydantic schema module, pytest suite, Ruff/mypy pre-commit checks, GitHub Actions CI, MIT license, and README. This establishes a minimal installable/testable baseline before data and pipeline implementation. | Project owner |
| 2026-11-26 | Added canonical Pydantic schemas for P1-T1. | Mutation, copy-number segment, survival, and clinical-record models validate required identifiers and numeric bounds; chromosome and allele normalization is consistent; unknown fields are rejected. | Project owner |
| 2026-11-26 | Corrected and verified API-only GDC ingestion for P1-T2. | The ingester filters to open masked MAF and allele-specific segment file types, streams source content through the GDC API into memory without caching raw data locally, fetches public clinical fields through the Cases API, and writes canonical mutation, segment, clinical, and cohort tables with source checksums in a manifest. | Project owner |
| 2026-11-26 | Implemented cBioPortal public-study ingestion for P1-T3. | The ingester reads study metadata and molecular profiles, fetches detailed mutations and discrete CNA calls, paginates clinical data, caches API responses, exports canonical mutation/clinical TSVs and CNA calls, and writes a SHA-256 manifest. | Project owner |
| 2026-11-26 | Implemented a processed GEO/longitudinal loader for P1-T4. | The loader requires a valid GSE accession and verified license URL, excludes BAM/CRAM/SRA raw reads, caches processed files by SHA-256, exports canonical mutation/segment/clinical tables, and records source provenance in a manifest. | Project owner |
| 2026-11-26 | Added cohort table and data dictionary support for P1-T5. | A source-independent cohort builder derives sample union counts and clinical fields from canonical TSV tables; `docs/data_dictionary.md` documents table columns, units, identifiers, and open-data boundaries. | Project owner |
