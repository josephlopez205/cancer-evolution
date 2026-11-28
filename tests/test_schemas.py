import pytest
from pydantic import ValidationError

from clonetrace.schemas import (
    Clinical,
    ClinicalRecord,
    CopyNumberSegment,
    CopyNumberTable,
    Mutation,
    MutationRecord,
    PurityPloidyEstimate,
)


def test_mutation_record_normalizes_canonical_fields() -> None:
    mutation = MutationRecord(
        sample_id=" TCGA-1 ",
        gene=" tp53 ",
        chrom="chr17",
        pos=7579472,
        ref="c",
        alt="g",
        ref_count=40,
        alt_count=10,
    )

    assert mutation.sample == "TCGA-1"
    assert mutation.gene == "tp53"
    assert mutation.chrom == "17"
    assert mutation.ref == "C"
    assert mutation.alt == "G"


def test_mutation_alias_matches_record_class() -> None:
    mutation = Mutation(
        sample="sample-1",
        gene="KRAS",
        chrom="12",
        pos=25398284,
        ref="C",
        alt="A",
        ref_count=20,
        alt_count=5,
    )

    assert isinstance(mutation, MutationRecord)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("pos", 0),
        ("ref_count", -1),
        ("alt_count", -1),
    ],
)
def test_mutation_rejects_invalid_numbers(field: str, value: int) -> None:
    data = {
        "sample": "sample-1",
        "gene": "TP53",
        "chrom": "17",
        "pos": 100,
        "ref": "A",
        "alt": "C",
        "ref_count": 20,
        "alt_count": 5,
    }
    data[field] = value

    with pytest.raises(ValidationError):
        MutationRecord.model_validate(data)


def test_copy_number_segment_normalizes_and_validates() -> None:
    segment = CopyNumberSegment(
        sample="sample-1",
        chrom="chrX",
        start=100,
        end=500,
        total_cn=4.0,
        minor_cn=1.0,
    )

    assert segment.chrom == "X"
    assert CopyNumberTable(records=[segment]).records == [segment]


@pytest.mark.parametrize(
    "values",
    [
        {"start": 500, "end": 500},
        {"start": 100, "end": 500, "total_cn": 1.0, "minor_cn": 2.0},
    ],
)
def test_copy_number_segment_rejects_invalid_bounds(values: dict[str, int | float]) -> None:
    data: dict[str, str | int | float] = {
        "sample": "sample-1",
        "chrom": "1",
        "start": 100,
        "end": 500,
        "total_cn": 2.0,
        "minor_cn": 1.0,
    }
    data.update(values)

    with pytest.raises(ValidationError):
        CopyNumberSegment.model_validate(data)


def test_clinical_record_accepts_flat_and_nested_survival() -> None:
    flat = ClinicalRecord(
        sample="sample-1",
        survival_time=36.5,
        survival_status=1,
        therapy=["tamoxifen", "tamoxifen"],
    )
    assert flat.therapy == ["tamoxifen"]

    nested = Clinical(
        sample_id="sample-1",
        survival={"time": 36.5, "status": "event"},
        therapy=["tamoxifen", "radiotherapy"],
    )

    assert isinstance(nested, ClinicalRecord)
    assert nested.survival_time == 36.5
    assert nested.survival_status == "event"
    assert nested.therapy == ["tamoxifen", "radiotherapy"]


def test_records_reject_unknown_fields() -> None:
    with pytest.raises(ValidationError):
        MutationRecord(
            sample="sample-1",
            gene="TP53",
            chrom="17",
            pos=100,
            ref="A",
            alt="C",
            ref_count=20,
            alt_count=5,
            unknown="value",
        )


def test_purity_ploidy_estimate_rejects_out_of_range_values() -> None:
    with pytest.raises(ValidationError):
        PurityPloidyEstimate(
            sample="sample-1",
            purity=1.1,
            source="estimate.tsv",
        )

    with pytest.raises(ValidationError):
        PurityPloidyEstimate(
            sample="sample-1",
            ploidy=0,
            source="estimate.tsv",
        )
