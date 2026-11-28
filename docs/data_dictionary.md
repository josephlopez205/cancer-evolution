# CloneTrace Data Dictionary

This dictionary defines the canonical cohort tables produced by Phase 1 ingestion.
All tabular files are UTF-8 TSV files. Missing values are left blank.

## `mutations.tsv`

Canonical mutation records are one row per somatic mutation observation.

| Column | Type | Required | Definition |
|---|---|---:|---|
| `sample` | string | yes | Sample or aliquot identifier |
| `gene` | string | yes | HGNC or source-provided gene symbol |
| `chrom` | string | yes | Chromosome, normalized without a `chr` prefix |
| `pos` | integer | yes | 1-based genomic position |
| `ref` | string | yes | Reference allele in uppercase |
| `alt` | string | yes | Alternate allele in uppercase |
| `ref_count` | integer | yes | Reference read count; non-negative |
| `alt_count` | integer | yes | Alternate read count; non-negative |

## `copy_number_segments.tsv`

Canonical segment records are one row per copy-number segment.

| Column | Type | Required | Definition |
|---|---|---:|---|
| `sample` | string | yes | Sample or aliquot identifier |
| `chrom` | string | yes | Chromosome, normalized without a `chr` prefix |
| `start` | integer | yes | Segment start position; zero-based or source-native but not negative |
| `end` | integer | yes | Segment end position; greater than `start` |
| `total_cn` | float | yes | Estimated total copy number; non-negative |
| `minor_cn` | float | yes | Estimated minor allele copy number; no greater than `total_cn` |

## `purity_ploidy.tsv`

Optional sample-level purity/ploidy estimates supplied by an external method or source. At least
one of `purity` and `ploidy` must be provided on each row. Re-estimation is not performed by this
input loader.

| Column | Type | Required | Definition |
|---|---|---:|---|
| `sample` | string | yes | Sample or aliquot identifier |
| `purity` | float or blank | yes | Tumor purity estimate in the range [0, 1] |
| `ploidy` | float or blank | yes | Positive average tumor ploidy estimate |
| `source` | string | yes | Provenance, such as the source file or tool output |
| `method` | string or blank | no | Estimation method, when known |

Purity-adjusted cancer cell fraction uses the diploid-normal model:
`CCF = VAF * (purity * total_cn + 2 * (1 - purity)) / (purity * multiplicity)`.
The inverse expected-VAF calculation uses the same assumptions. CCF estimates above 1 are retained
and should be interpreted as sampling noise or a model mismatch, not silently clipped.

## `copy_number_calls.tsv`

Curated discrete gene-level CNA calls, currently supplied by cBioPortal.

| Column | Type | Required | Definition |
|---|---|---:|---|
| `sample` | string | yes | Sample identifier |
| `gene` | string | yes | HGNC or source-provided gene symbol |
| `entrez_gene_id` | integer or blank | no | Entrez gene identifier |
| `cna_value` | integer or blank | no | Source discrete CNA value |
| `molecular_profile_id` | string or blank | no | cBioPortal molecular profile |

## `clinical.tsv`

Canonical clinical records are one row per sample.

| Column | Type | Required | Definition |
|---|---|---:|---|
| `sample` | string | yes | Sample or patient identifier |
| `survival_time` | float or blank | no | Overall survival time in the source's native unit |
| `survival_status` | integer, `event`, `censored`, or blank | no | `1`/event versus `0`/censored |
| `therapy` | comma-separated string | no | Normalized therapy labels |

The source-specific unit for `survival_time` must be preserved and documented in the source manifest.
GDC-derived `survival_time` is in days (`days_to_death` or `days_to_last_follow_up`); the loader does
not convert it.
P1-T5 does not infer survival status where the source does not provide a valid outcome.

## `cohort.tsv`

The harmonized sample-level cohort table is derived from the canonical tables above.

| Column | Type | Required | Definition |
|---|---|---:|---|
| `sample` | string | yes | Unique sample identifier in this cohort |
| `source` | string | yes | Source label, such as `gdc`, `cbioportal`, or `geo` |
| `mutation_count` | integer | yes | Number of canonical mutation rows for the sample |
| `segment_count` | integer | yes | Number of segment or gene-CNA rows for the sample |
| `survival_time` | float or blank | no | Source survival time |
| `survival_status` | integer, `event`, `censored`, or blank | no | Survival outcome |
| `therapy` | comma-separated string | no | Normalized therapy labels |

## Manifests

Source manifests remain the provenance record for downloaded or loaded files. They include source
metadata, accession or project identifiers, source/derived paths, sizes, and available checksums.
Missing or controlled-access files are never treated as successful cohort records.

## Identifier and unit rules

- Keep source identifiers stable; do not silently remap patient and sample IDs.
- Normalize chromosome spelling, but do not normalize genomic coordinates across assemblies.
- GDC and cBioPortal reference genomes must be recorded with source metadata.
- No BAM, CRAM, SRA, controlled-access, or COSMIC bulk input belongs in these tables.
- Sort mutation and segment output deterministically by source encounter order or sample/index order.
