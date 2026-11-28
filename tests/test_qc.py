import pytest
from pydantic import ValidationError

from clonetrace.qc import (
    estimate_ccf_from_vaf,
    expected_vaf,
    read_purity_ploidy_estimates,
)
from clonetrace.schemas import PurityPloidyEstimate


def test_purity_ploidy_estimate_validates_provenance_and_bounds() -> None:
    estimate = PurityPloidyEstimate(
        sample_id="sample-1",
        purity=0.72,
        ploidy=2.4,
        source="facets-results.tsv",
        method="FACETS",
    )

    assert estimate.sample == "sample-1"
    assert estimate.purity == 0.72
    assert estimate.ploidy == 2.4

    with pytest.raises(ValidationError):
        PurityPloidyEstimate(sample="sample-1", source="unknown")


def test_read_purity_ploidy_estimates_accepts_partial_estimates(tmp_path) -> None:
    path = tmp_path / "purity_ploidy.tsv"
    path.write_text(
        "sample\tpurity\tploidy\tsource\tmethod\n"
        "sample-1\t0.72\t2.4\tfacets.tsv\tFACETS\n"
        "sample-2\t\t3.1\tcBioPortal\t\n",
        encoding="utf-8",
    )

    estimates = read_purity_ploidy_estimates(path)

    assert len(estimates) == 2
    assert estimates[0].purity == 0.72
    assert estimates[0].ploidy == 2.4
    assert estimates[1].purity is None
    assert estimates[1].ploidy == 3.1
    assert estimates[1].method is None


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        (
            "sample\tpurity\tploidy\tsource\n"
            "sample-1\t\t\tunknown\n",
            "at least one of purity or ploidy",
        ),
        (
            "sample\tpurity\tploidy\tsource\n"
            "sample-1\t0.5\t2.1\ttool\n"
            "sample-1\t0.6\t2.2\ttool\n",
            "duplicate purity/ploidy estimate",
        ),
        (
            "sample\tpurity\tsource\n"
            "sample-1\t0.5\ttool\n",
            "missing columns: ploidy",
        ),
        (
            "sample\tpurity\tploidy\tsource\n"
            "sample-1\tnot-a-number\t2.1\ttool\n",
            "purity must be numeric",
        ),
    ],
)
def test_read_purity_ploidy_estimates_rejects_invalid_rows(
    tmp_path,
    contents: str,
    message: str,
) -> None:
    path = tmp_path / "invalid.tsv"
    path.write_text(contents, encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        read_purity_ploidy_estimates(path)


def test_expected_vaf_and_ccf_inverse_for_diploid_sample() -> None:
    vaf = expected_vaf(ccf=0.5, purity=0.8, total_cn=2)

    assert vaf == pytest.approx(0.2)
    assert estimate_ccf_from_vaf(vaf=vaf, purity=0.8, total_cn=2) == pytest.approx(0.5)


def test_vaf_model_accounts_for_multiplicity_after_copy_gain() -> None:
    vaf = expected_vaf(ccf=0.5, purity=0.8, total_cn=4, multiplicity=2)

    assert vaf == pytest.approx(2 / 9)
    assert (
        estimate_ccf_from_vaf(vaf=vaf, purity=0.8, total_cn=4, multiplicity=2)
        == pytest.approx(0.5)
    )


def test_ccf_estimates_above_one_are_not_silently_clipped() -> None:
    assert estimate_ccf_from_vaf(vaf=0.6, purity=0.8, total_cn=2) == pytest.approx(1.5)


@pytest.mark.parametrize(
    ("function", "kwargs"),
    [
        (expected_vaf, {"ccf": 1.1, "purity": 0.8, "total_cn": 2}),
        (expected_vaf, {"ccf": 0.5, "purity": 0.8, "total_cn": 2, "multiplicity": 3}),
        (estimate_ccf_from_vaf, {"vaf": 0.4, "purity": 0, "total_cn": 2}),
        (estimate_ccf_from_vaf, {"vaf": 1.1, "purity": 0.8, "total_cn": 2}),
    ],
)
def test_vaf_functions_reject_invalid_inputs(function, kwargs: dict[str, float]) -> None:
    with pytest.raises(ValueError):
        function(**kwargs)
