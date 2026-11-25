import json

import pytest

from clonetrace.ingest.cbioportal import (
    ingest_cbioportal,
    parse_clinical,
    parse_discrete_cna,
    parse_mutations,
)
from clonetrace.schemas import ClinicalRecord, MutationRecord


def test_parse_mutations() -> None:
    rows = [
        {
            "sampleId": "MB-0000",
            "gene": {"hugoGeneSymbol": "TP53"},
            "chromosome": "17",
            "startPosition": 7579472,
            "referenceAllele": "C",
            "variantAllele": "T",
            "tumorAltCount": 8,
            "normalRefCount": 40,
        }
    ]

    mutations = parse_mutations(rows)

    assert len(mutations) == 1
    assert isinstance(mutations[0], MutationRecord)
    assert mutations[0].gene == "TP53"
    assert mutations[0].ref_count == 40
    assert mutations[0].alt_count == 8


def test_parse_clinical() -> None:
    rows = [
        {"sampleId": "MB-0000", "clinicalAttributeId": "OS_MONTHS", "value": "36.5"},
        {"sampleId": "MB-0000", "clinicalAttributeId": "OS_STATUS", "value": "1:DECEASED"},
        {"sampleId": "MB-0000", "clinicalAttributeId": "CHEMOTHERAPY", "value": "YES"},
    ]

    records = parse_clinical(rows)

    assert len(records) == 1
    assert isinstance(records[0], ClinicalRecord)
    assert records[0].survival_time == 36.5
    assert records[0].survival_status == "event"
    assert records[0].therapy == ["YES"]


def test_parse_discrete_cna() -> None:
    rows = [
        {
            "sampleId": "MB-0000",
            "gene": {"hugoGeneSymbol": "ERBB2"},
            "entrezGeneId": 2064,
            "value": 2,
            "molecularProfileId": "brca_metabric_cna",
        }
    ]

    records = parse_discrete_cna(rows)

    assert records == [
        {
            "sample": "MB-0000",
            "gene": "ERBB2",
            "entrez_gene_id": 2064,
            "cna_value": 2,
            "molecular_profile_id": "brca_metabric_cna",
        }
    ]


def test_ingest_cbioportal_writes_tables_and_manifest(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cache = tmp_path / "cache"
    output = tmp_path / "processed"

    monkeypatch.setattr(
        "clonetrace.ingest.cbioportal.get_study",
        lambda study: {"studyId": study, "publicStudy": True},
    )
    monkeypatch.setattr(
        "clonetrace.ingest.cbioportal.get_molecular_profiles",
        lambda study: [
            {
                "molecularProfileId": "brca_metabric_mutations",
                "molecularAlterationType": "MUTATION_EXTENDED",
            },
            {
                "molecularProfileId": "brca_metabric_cna",
                "molecularAlterationType": "COPY_NUMBER_ALTERATION",
            },
        ],
    )
    monkeypatch.setattr(
        "clonetrace.ingest.cbioportal.get_mutations",
        lambda study, profile: [
            {
                "sampleId": "MB-0000",
                "gene": {"hugoGeneSymbol": "TP53"},
                "chromosome": "17",
                "startPosition": 7579472,
                "referenceAllele": "C",
                "variantAllele": "T",
                "tumorAltCount": 8,
                "normalRefCount": 40,
            }
        ],
    )
    monkeypatch.setattr(
        "clonetrace.ingest.cbioportal.get_discrete_cna",
        lambda study, profile: [
            {
                "sampleId": "MB-0000",
                "gene": {"hugoGeneSymbol": "ERBB2"},
                "entrezGeneId": 2064,
                "value": 2,
                "molecularProfileId": "brca_metabric_cna",
            }
        ],
    )
    monkeypatch.setattr(
        "clonetrace.ingest.cbioportal.get_clinical_data",
        lambda study: [
            {"sampleId": "MB-0000", "clinicalAttributeId": "OS_MONTHS", "value": "36.5"},
            {"sampleId": "MB-0000", "clinicalAttributeId": "OS_STATUS", "value": "1:DECEASED"},
        ],
    )

    result = ingest_cbioportal("brca_metabric", cache_dir=cache, output_dir=output)

    assert result.mutation_count == 1
    assert result.clinical_count == 1
    assert result.cna_count == 1
    assert (output / "mutations.tsv").exists()
    assert (output / "clinical.tsv").exists()
    assert (output / "copy_number_calls.tsv").exists()
    assert json.loads(result.manifest.read_text(encoding="utf-8"))["source"] == "cBioPortal"
