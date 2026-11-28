"""Purity- and copy-number-aware variant allele fraction calculations."""

import math


def _validate_inputs(
    *,
    purity: float,
    total_cn: float,
    multiplicity: int,
) -> None:
    if not math.isfinite(purity) or not 0 <= purity <= 1:
        raise ValueError("purity must be finite and between 0 and 1")
    if not math.isfinite(total_cn) or total_cn <= 0:
        raise ValueError("total_cn must be finite and greater than 0")
    if isinstance(multiplicity, bool) or not isinstance(multiplicity, int) or multiplicity < 1:
        raise ValueError("multiplicity must be a positive integer")
    if multiplicity > total_cn:
        raise ValueError("multiplicity must not exceed total_cn")


def expected_vaf(
    *,
    ccf: float,
    purity: float,
    total_cn: float,
    multiplicity: int = 1,
) -> float:
    """Calculate expected VAF under a diploid normal-contamination model.

    CCF is the fraction of cancer cells carrying the mutation. Normal cells
    contribute two copies; tumor cells contribute ``total_cn`` copies, of which
    ``multiplicity`` carry the variant.
    """

    _validate_inputs(purity=purity, total_cn=total_cn, multiplicity=multiplicity)
    if not math.isfinite(ccf) or not 0 <= ccf <= 1:
        raise ValueError("ccf must be finite and between 0 and 1")
    denominator = purity * total_cn + 2 * (1 - purity)
    if denominator <= 0:
        raise ValueError("purity and copy number produce an invalid denominator")
    return purity * ccf * multiplicity / denominator


def estimate_ccf_from_vaf(
    *,
    vaf: float,
    purity: float,
    total_cn: float,
    multiplicity: int = 1,
) -> float:
    """Invert the expected-VAF model to estimate cancer cell fraction.

    The returned estimate is not clipped to 1. Values above 1 should remain
    visible as evidence of sampling noise or a misspecified copy-number model.
    """

    _validate_inputs(purity=purity, total_cn=total_cn, multiplicity=multiplicity)
    if not math.isfinite(vaf) or not 0 <= vaf <= 1:
        raise ValueError("vaf must be finite and between 0 and 1")
    if purity == 0:
        raise ValueError("purity must be greater than 0 to estimate CCF from VAF")
    denominator = purity * multiplicity
    numerator = vaf * (purity * total_cn + 2 * (1 - purity))
    return numerator / denominator
