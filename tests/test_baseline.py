import pytest

from clonetrace.cluster import fit_binomial_mixture
from clonetrace.schemas import MutationRecord


def _mutations(alt_counts: list[int], sample: str = "sample-1") -> list[MutationRecord]:
    return [
        MutationRecord(
            sample=sample,
            gene=f"GENE{index}",
            chrom="1",
            pos=index + 1,
            ref="A",
            alt="C",
            ref_count=100 - alt_count,
            alt_count=alt_count,
        )
        for index, alt_count in enumerate(alt_counts)
    ]


ZERO_DEPTH_MUTATION = MutationRecord(
    sample="sample-1",
    gene="GENE0",
    chrom="1",
    pos=1,
    ref="A",
    alt="C",
    ref_count=0,
    alt_count=0,
)


def test_binomial_mixture_recovers_two_separated_vaf_groups() -> None:
    mutations = _mutations([9, 10, 11, 12] * 4 + [33, 35, 36, 37] * 4)

    result = fit_binomial_mixture(mutations, max_clusters=4)

    assert result.sample == "sample-1"
    assert result.method == "baseline"
    assert result.cluster_count == 2
    assert result.converged
    assert result.clusters[0].posterior_mean_vaf < result.clusters[1].posterior_mean_vaf
    assert result.clusters[0].mutation_count == 16
    assert result.clusters[1].mutation_count == 16
    assert all(assignment.posterior_probability > 0.9 for assignment in result.assignments)
    assert all(
        sum(assignment.cluster_probabilities) == pytest.approx(1.0)
        for assignment in result.assignments
    )


def test_binomial_mixture_selects_one_cluster_for_unimodal_data() -> None:
    result = fit_binomial_mixture(_mutations([19, 20, 20, 21, 22] * 4), max_clusters=4)

    assert result.cluster_count == 1
    assert result.clusters[0].mutation_count == 20
    assert all(
        assignment.cluster_probabilities == pytest.approx((1.0,))
        for assignment in result.assignments
    )


def test_binomial_mixture_is_deterministic() -> None:
    mutations = _mutations([8, 10, 11, 35, 36, 38])

    assert fit_binomial_mixture(mutations) == fit_binomial_mixture(mutations)


@pytest.mark.parametrize(
    ("mutations", "kwargs", "message"),
    [
        ([], {}, "mutations must not be empty"),
        (_mutations([10, 20], sample="one") + _mutations([30], sample="two"), {}, "same sample"),
        ([ZERO_DEPTH_MUTATION], {}, "zero total read depth"),
        (_mutations([10]), {"max_clusters": 0}, "max_clusters"),
        (_mutations([10]), {"max_iterations": 0}, "max_iterations"),
    ],
)
def test_binomial_mixture_rejects_invalid_inputs(mutations, kwargs, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        fit_binomial_mixture(mutations, **kwargs)
