"""Metrics for comparing simulated and inferred clonal results."""

from collections.abc import Sequence
from dataclasses import dataclass
from math import comb


def clone_count_error(true_count: int, predicted_count: int) -> int:
    """Return predicted clone count minus true clone count."""

    if true_count < 0 or predicted_count < 0:
        raise ValueError("clone counts must be non-negative")
    return predicted_count - true_count


def _validate_labels(labels: Sequence[int], name: str) -> None:
    if len(labels) == 0:
        raise ValueError(f"{name} must not be empty")
    if any(label < 0 for label in labels):
        raise ValueError(f"{name} values must be non-negative")


def adjusted_rand_index(true_labels: Sequence[int], predicted_labels: Sequence[int]) -> float:
    """Compute the adjusted Rand index for two partitions."""

    _validate_labels(true_labels, "true_labels")
    _validate_labels(predicted_labels, "predicted_labels")
    if len(true_labels) != len(predicted_labels):
        raise ValueError("label sequences must have the same length")

    n = len(true_labels)
    if n == 1:
        return 1.0

    contingency: dict[tuple[int, int], int] = {}
    true_counts: dict[int, int] = {}
    predicted_counts: dict[int, int] = {}
    for true_label, predicted_label in zip(true_labels, predicted_labels, strict=True):
        key = (true_label, predicted_label)
        contingency[key] = contingency.get(key, 0) + 1
        true_counts[true_label] = true_counts.get(true_label, 0) + 1
        predicted_counts[predicted_label] = predicted_counts.get(predicted_label, 0) + 1

    pairs_in_contingency = sum(comb(count, 2) for count in contingency.values())
    pairs_in_true = sum(comb(count, 2) for count in true_counts.values())
    pairs_in_predicted = sum(
        comb(count, 2) for count in predicted_counts.values()
    )
    pairs_total = comb(n, 2)
    expected_index = (
        pairs_in_true * pairs_in_predicted / pairs_total
        if pairs_total
        else 0.0
    )
    maximum_index = (pairs_in_true + pairs_in_predicted) / 2
    denominator = maximum_index - expected_index
    if denominator == 0:
        equivalent = all(
            (true_labels[i] == true_labels[j])
            == (predicted_labels[i] == predicted_labels[j])
            for i in range(n)
            for j in range(i + 1, n)
        )
        return 1.0 if equivalent else 0.0
    return (pairs_in_contingency - expected_index) / denominator


def _ancestor_pairs(parents: Sequence[int]) -> set[tuple[int, int]]:
    """Return (ancestor, descendant) pairs, including transitive ancestors."""

    if len(parents) == 0:
        raise ValueError("parents must not be empty")

    result: set[tuple[int, int]] = set()
    for descendant in range(len(parents)):
        seen: set[int] = set()
        ancestor = parents[descendant]
        while ancestor != descendant and 0 <= ancestor < len(parents) and ancestor not in seen:
            seen.add(ancestor)
            result.add((ancestor, descendant))
            ancestor = parents[ancestor]
    return result


def ancestor_descendant_accuracy(
    true_parents: Sequence[int],
    predicted_parents: Sequence[int],
) -> float:
    """Return the fraction of ordered node pairs with matching ancestry."""

    if len(true_parents) != len(predicted_parents):
        raise ValueError("parent sequences must have the same length")
    if len(true_parents) == 0:
        raise ValueError("parent sequences must not be empty")
    if len(true_parents) == 1:
        return 1.0

    true_pairs = _ancestor_pairs(true_parents)
    predicted_pairs = _ancestor_pairs(predicted_parents)
    node_count = len(true_parents)
    ordered_pair_count = node_count * (node_count - 1)
    matching = 0
    for ancestor in range(node_count):
        for descendant in range(node_count):
            if ancestor == descendant:
                continue
            pair = (ancestor, descendant)
            matching += int((pair in true_pairs) == (pair in predicted_pairs))
    return matching / ordered_pair_count


@dataclass(frozen=True)
class CalibrationResult:
    """Uncertainty calibration summary."""

    brier_score: float
    expected_calibration_error: float
    bin_count: int


def uncertainty_calibration(
    probabilities: Sequence[float],
    outcomes: Sequence[int],
    *,
    bin_count: int = 10,
) -> CalibrationResult:
    """Measure probability calibration using Brier score and equal-width ECE."""

    if len(probabilities) != len(outcomes):
        raise ValueError("probabilities and outcomes must have the same length")
    if not probabilities:
        raise ValueError("probabilities must not be empty")
    if bin_count < 1:
        raise ValueError("bin_count must be at least 1")
    if any(probability < 0 or probability > 1 for probability in probabilities):
        raise ValueError("probabilities must be between 0 and 1")
    if any(outcome not in (0, 1) for outcome in outcomes):
        raise ValueError("outcomes must contain only 0 or 1")

    brier = sum(
        (probability - outcome) ** 2
        for probability, outcome in zip(probabilities, outcomes, strict=True)
    ) / len(probabilities)

    ece = 0.0
    total = len(probabilities)
    for bin_index in range(bin_count):
        lower = bin_index / bin_count
        upper = (bin_index + 1) / bin_count
        members = [
            (probability, outcome)
            for probability, outcome in zip(probabilities, outcomes, strict=True)
            if probability >= lower
            and (probability < upper or (bin_index == bin_count - 1 and probability <= upper))
        ]
        if not members:
            continue
        mean_probability = sum(probability for probability, _ in members) / len(members)
        mean_outcome = sum(outcome for _, outcome in members) / len(members)
        ece += (len(members) / total) * abs(mean_probability - mean_outcome)

    return CalibrationResult(
        brier_score=brier,
        expected_calibration_error=ece,
        bin_count=bin_count,
    )
