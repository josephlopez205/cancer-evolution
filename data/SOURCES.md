# Data Sources and Access Terms

**Inventory date:** 2026-11-26  
**Policy:** Open data only. Raw reads, BAM/CRAM files, controlled-access files, and COSMIC bulk downloads are out of scope.

This file records the planned source inventory for data sources D1-D5. It is not evidence that a particular record is suitable for ingestion. A source-specific record or study must be rechecked and added here before download.

| ID | Source | Intended use | Candidate identifiers | Access and license status | Verification URL(s) | Reverification |
|---|---|---|---|---|---|---|
| D1 | GDC / TCGA open tier | Masked somatic MAF, allele-specific copy-number segments, clinical data | Primary project: `TCGA-BRCA` | Public GDC API access only. File metadata must have `access=open`; source files are streamed from the API and processed in memory, never cached or persisted locally. BAMs and controlled-access files are excluded. GDC and TCGA publication/use terms should be rechecked before release. | [GDC access policies](https://gdc.cancer.gov/access-data/data-access-policies), [GDC access-process overview](https://gdc.cancer.gov/access-data/data-access-processes-and-tools), [GDC API](https://api.gdc.cancer.gov/) | Before P1-T2 and P8-T3 |
| D2 | cBioPortal public studies | Curated breast mutations, CNA, and clinical data | Primary validation study: `brca_metabric`; primary comparison: `brca_tcga_pan_can_atlas_2018` | Public API access. The cBioPortal software is Apache License 2.0. Study-level source licenses, publication terms, and citation requirements differ; METABRIC and TCGA source terms must be recorded before ingestion. | [cBioPortal API](https://www.cbioportal.org/api), [cBioPortal license](https://github.com/cBioPortal/cbioportal/blob/master/LICENSE), [METABRIC study](https://www.cbioportal.org/study/summary?id=brca_metabric), [TCGA PanCancer Atlas breast study](https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018) | Before P1-T3 and P8-T3 |
| D3 | PCAWG open tier | Consensus calls and published timing for comparison | PCAWG breast cancer subset, after exact open record selection | In scope only for published/open-tier processed calls. The ICGC portal states that current PCAWG data availability is for authorized users; no PCAWG data is to be ingested until a specific open-access record is identified and independently verified. | [ICGC 25K/PCAWG access notice](https://dcc.icgc.org/), [PCAWG publication](https://doi.org/10.1038/s41586-020-1969-6) | Before P5-T2; no download pending verification |
| D4 | GEO / SRA / Zenodo | Processed multi-region or longitudinal breast cancer datasets | No dataset selected; P1-T4 provides a verified local processed-file loader for GEO `GSEnnnn` accession IDs. SRA and Zenodo remain unresolved until exact record/DOI review. | Not approved for ingestion until the exact dataset is selected and a license URL is supplied to the loader. GEO/SRA public availability is not by itself a license or reuse permission. NCBI molecular data generally places no NCBI-level reuse restriction, but submitter or source terms still control; each GEO/SRA study and any Zenodo record must have its own license/citation recorded. Zenodo requires each uploaded item's license conditions to be respected. | [GEO](https://www.ncbi.nlm.nih.gov/geo/), [SRA](https://www.ncbi.nlm.nih.gov/sra), [NCBI data usage policies](https://www.ncbi.nlm.nih.gov/home/about/policies/), [Zenodo terms](https://zenodo.org/terms) | Before P1-T4 and for every new record |
| D5 | DepMap / CCLE | Cell-line context for driver validation | No release selected | Open planned context source, but a specific release and CCLE data license must be recorded before download. The portal provides bulk downloads but access terms and source release documentation must be verified per release. | [DepMap portal](https://depmap.org/portal/), [DepMap data page](https://depmap.org/portal/data_page/?tab=allData), [CCLE](https://depmap.org/portal/cell_line/) | Before P6-T1 or any P5/P6 ingestion |

## Initial availability snapshot

The following counts were obtained from the GDC and cBioPortal APIs on 2026-11-26 and are retained only as selection evidence:

| Cancer | GDC cases | GDC Masked Somatic Mutation records | GDC open copy-number records | GDC open clinical records | cBioPortal PanCancer Atlas samples | cBioPortal METABRIC samples |
|---|---:|---:|---:|---:|---:|---:|
| Breast invasive carcinoma (TCGA-BRCA) | 1098 | 992 | 11,056 | 2,288 | 1084 | 2509 |
| Lung adenocarcinoma (TCGA-LUAD) | 585 | 618 | 5328 | 1146 | 566 | - |
| Colorectal adenocarcinoma (TCGA-COAD/COADREAD) | 461 | 457 | 4655 | 995 | 594 | - |

These counts describe file or API records, not unique eligible patients, and do not replace data-cleaning and sample-exclusion checks.

## D1-D5 release checklist

- [ ] Use only public/open access tiers; verify `access=open` or equivalent source-level status for every file.
- [ ] Select one exact study or release before download.
- [ ] Record the exact URL, DOI or accession, license, citation requirement, access date, release/version, and checksum.
- [ ] Review whether clinical and survival fields permit the intended analysis.
- [ ] Explicitly exclude BAM/CRAM and any controlled-access or embargoed data.
- [ ] Add the source to this table and update the re-verification date.
- [ ] For derived or redistributed data, preserve upstream citations and satisfy study-level attribution terms.

## Decision boundary

P0-T3 establishes an initial source inventory and access boundary. D4 is intentionally unresolved until a specific processed breast cancer dataset is selected, and D3 is excluded until an identifiable open-tier PCAWG record is confirmed. This avoids treating repository or portal availability as proof of license or open access.
