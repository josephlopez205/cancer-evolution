import pytest

from clonetrace.benchmark import (
    adjusted_rand_index,
    ancestor_descendant_accuracy,
    clone_count_error,
    uncertainty_calibration,
)


def test_clone_count_error_is_signed() -> None:
    assert clone_count_error(3, 5) == 2
    assert clone_count_error(3, 1) == -2
    assert clone_count_error(3, 3) == 0


def test_adjusted_rand_index_identical_and_opposite_partitions() -> None:
    labels = [0, 0, 1, 1]
    assert adjusted_rand_index(labels, labels) == pytest.approx(1.0)
    assert adjusted_rand_index(labels, [1, 1, 0, 0]) == pytest.approx(1.0)
    assert adjusted_rand_index([0, 0, 0, 0], [0, 1, 2, 3]) == pytest.approx(0.0)


def test_ancestor_descendant_accuracy_for_toy_trees() -> None:
    assert ancestor_descendant_accuracy([0, 0, 0], [0, 0, 0]) == 1.0
    assert ancestor_descendant_accuracy([0, 0, 0], [0, 0, 1]) == pytest.approx(
        5 / 6
    )


def test_uncertainty_calibration_hand_computed() -> None:
    result = uncertainty_calibration(
        [0.1, 0.2, 0.8, 0.9],
        [0, 1, 1, 1],
        bin_count=2,
    )

    assert result.brier_score == pytest.approx(
        ((0.1 - 0) ** 2 + (0.2 - 1) ** 2 + (0.8 - 1) ** 2 + (0.9 - 1) ** 2)
        / 4
    )
    assert result.expected_calibration_error == pytest.approx(0.25)
    assert result.bin_count == 2


def test_metrics_reject_invalid_lengths_and_probabilities() -> None:
    with pytest.raises(ValueError):
        adjusted_rand_index([0, 1], [0])
    with pytest.raises(ValueError):
        ancestor_descendant_accuracy([0, 0], [0])
    with pytest.raises(ValueError):
        uncertainty_calibration([1.2], [1])
