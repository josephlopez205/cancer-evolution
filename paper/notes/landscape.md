# Literature and Tool Landscape

**Status:** P0-T1 initial review  
**Accessed:** 2026-11-26  
**Scope:** Clustering, subclonal reconstruction, and phylogenetic inference for this project's bulk-tumor pipeline.

## Comparison table

| Tool | Best fit for | Primary method | Inputs | Outputs | Strengths | Limits / integration concerns | Evidence |
|---|---|---|---|---|---|---|---|
| PyClone-VI | Fast primary mutation clustering | Variational inference; beta-binomial default in the repository example | Mutation copy number, read counts, and optional binomial/beta-binomial model inputs | Cluster labels and cellular prevalence estimates | Scalable, Python-native, configurable restarts, direct fit/write-results CLI | No phylogenetic tree inference; performance and uncertainty depend on copy-number and read-count quality | [PyClone-VI repository](https://github.com/Roth-Lab/pyclone-vi); [paper, BMC Bioinformatics 2020](https://doi.org/10.1186/s12859-020-03919-2) |
| MOBSTER | Independent baseline clustering with evolutionary structure | Maximum-likelihood Dirichlet mixture models plus population-genetics theory; bootstrap confidence assessment | Bulk WGS mutation read/count data and clone-frequency information | Mixture components, bootstrap-supported parameter estimates, visualizations | Adds explicit subclonal-growth model and bootstrapping; useful contrast to generic Bayesian mixture clustering | R package; bootstrap and R interoperability add runtime and packaging cost; public API and fit object must be adapted to canonical schemas | [MOBSTER repository](https://github.com/caravagnalab/mobster); [Nature Genetics 2020](https://doi.org/10.1038/s41588-020-0675-5) |
| Pairtree | Multi-sample phylogenetic inference | Pairwise evolutionary relationship tensor followed by sampling and likelihood-scored phylogenies | Mutation or cluster presence, cellular prevalence estimates, and sample relationships | Probable trees and posterior-like support over reconstructions | Designed for multiple tissue samples and explicit uncertainty; jointly balances prevalence observations and ancestry constraints | Native C build and Unix-oriented tooling; more complex data contract than a one-sample minimal tree; availability/maintenance must be pinned | [Pairtree repository](https://github.com/morrislab/pairtree); [Cancer Biology & Discovery 2022](https://doi.org/10.1158/2643-3230.BCD-21-0092); [STAR Protocol 2022](https://pubmed.ncbi.nlm.nih.gov/36129821/) |
| CONIPHER | End-to-end comparison and multi-region reference | R implementation for mutation clustering and phylogenetic tree reconstruction; can use PyClone v0.13.1 for clustering | Mutation table and sample information | Clustering and tree-building outputs; can run both stages with a wrapper | Established end-to-end R workflow, useful published protocol and wrapper | R/Conda dependencies; combines a different clustering method and tree method, complicating method-specific comparisons; versions must be isolated and pinned | [CONIPHER repository](https://github.com/McGranahanLab/CONIPHER); [Nature Protocols 2023](https://doi.org/10.1038/s41596-023-00913-9) |

## Other methods and concepts considered

- **PyClone / classic PyClone variants:** older Bayesian single-sample approaches with mature comparative use, but PyClone-VI is the more suitable current Python wrapper target because it provides a scalable implementation and direct CLI.
- **Generic Bayesian mixture models:** retain as an in-repo baseline for interpretability and to test whether consensus gains survive a weaker, dependency-free method.
- **Sum-rule and ancestry-ordering methods:** useful as a deterministic fallback, especially for a single bulk sample, but cannot establish phylogenetic uncertainty as completely as a posterior-sampling or constraint-based method.
- **Cell-line and controlled-access references:** excluded from the primary open-data build. They may inform methods, but ingestion must remain restricted to open-tier data under recorded licenses.

## Recommended role in CloneTrace

1. **MOBSTER** and **PyClone-VI** form two independent clustering approaches, while `cluster/baseline.py` supplies a dependency-free third view.
2. **Pairtree** is the primary candidate for tree inference because it can propagate prevalence uncertainty and exploit multiple samples.
3. The minimal custom tree method remains required as a single-sample fallback and to keep the pipeline usable without native/optional dependencies.
4. **CONIPHER** is a candidate external benchmark and end-to-end reference, not a required core dependency for the first implementation. Any comparison should pin its R, package, wrapper, and PyClone versions.
5. Consensus should operate on mutation or cluster co-association across methods rather than comparing point-estimate cluster IDs, because cluster labels are not stable across methods.

## Open integration risks

- Copy-number correction and read-count definitions differ across tools; adapters must preserve a canonical mutation and segment representation.
- Pairtree can require native compilation and is UNIX-oriented. The packaged pipeline should expose it as an optional or isolated tool where necessary.
- MOBSTER is R-native. Its integration route must be validated separately and can be deferred if R packaging or bootstrap extraction is not reproducible enough.
- Published methods may report uncertainty, but only some quantities are directly comparable. Clone stability, tree posterior support, bootstrap confidence, and posterior probability must remain explicitly typed rather than being collapsed into an unsupported confidence score.
- Timing is not identifiable from every bulk sample. Any unadjusted VAF or point estimate must be retained only with a clear identifiability flag.

## Selection criteria for P0-T2

The final cancer type should be selected using these observations:

- Public open-tier MAF, segment, and clinical coverage.
- Sufficient purity and sequencing-depth annotation.
- Availability of known cancer drivers and a manageable driver panel.
- A validation cohort with independent published clonal calls or established benchmark labels.
- Multiple-sample or longitudinal data if available, because it is more informative than a single bulk sample for phylogenies.
- License, citation, and controlled-access constraints verified before acquisition.

## Source and license notes

This review documents software and publication sources. It does not establish usage terms for an implementation or dataset. Before ingesting D1-D5, verify the source license and access tier in `data/SOURCES.md`; no controlled-access data is included in this review.
