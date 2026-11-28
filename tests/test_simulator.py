import json

import pytest

from clonetrace.simulate import generate_scenario_grid, simulate


def test_simulator_is_deterministic(tmp_path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"

    result_one = simulate(
        first,
        seed=42,
        depth=60,
        purity=0.5,
        clone_count=4,
        topology="branching",
        mutation_count=25,
        cn_event_probability=0.5,
    )
    simulate(
        second,
        seed=42,
        depth=60,
        purity=0.5,
        clone_count=4,
        topology="branching",
        mutation_count=25,
        cn_event_probability=0.5,
    )

    for name in ("mutations.tsv", "copy_number_segments.tsv", "ground_truth.json"):
        assert (first / name).read_bytes() == (second / name).read_bytes()

    ground_truth = json.loads(result_one.ground_truth_path.read_text(encoding="utf-8"))
    assert ground_truth["format_version"] == "1"
    assert result_one.mutation_count == 25
    assert result_one.clone_count == 4
    assert len(ground_truth["mutation_assignments"]) == 25
    assert len(ground_truth["clone_tree"]["parents"]) == 4
    assert abs(sum(ground_truth["clone_tree"]["clone_fractions"]) - 1.0) < 1e-12


@pytest.mark.parametrize(
    ("topology", "expected_parent"),
    [
        ("linear", [0, 0, 1, 2]),
        ("branching", None),
    ],
)
def test_tree_topology(tmp_path, topology: str, expected_parent: list[int] | None) -> None:
    result = simulate(
        tmp_path,
        seed=7,
        clone_count=4,
        topology=topology,
        mutation_count=12,
    )
    ground_truth = json.loads(result.ground_truth_path.read_text(encoding="utf-8"))
    parents = ground_truth["clone_tree"]["parents"]

    assert parents[0] == 0
    if expected_parent is not None:
        assert parents == expected_parent


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("depth", 0),
        ("purity", 0),
        ("purity", 1.1),
        ("clone_count", 0),
        ("mutation_count", 0),
        ("cn_event_probability", -0.1),
    ],
)
def test_simulator_rejects_invalid_config(tmp_path, field: str, value: float) -> None:
    kwargs = {
        "depth": 100,
        "purity": 0.8,
        "clone_count": 3,
        "mutation_count": 20,
        "cn_event_probability": 0.25,
    }
    kwargs[field] = value

    with pytest.raises(ValueError):
        simulate(tmp_path, **kwargs)


def test_scenario_grid_is_reproducible_and_indexed(tmp_path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"

    result_one = generate_scenario_grid(
        first,
        seed=23,
        replicates=2,
        mutation_count=10,
        depths=(30, 60),
        purities=(0.3, 0.8),
        clone_counts=(1, 3),
        topologies=("linear", "branching"),
    )
    generate_scenario_grid(
        second,
        seed=23,
        replicates=2,
        mutation_count=10,
        depths=(30, 60),
        purities=(0.3, 0.8),
        clone_counts=(1, 3),
        topologies=("linear", "branching"),
    )

    assert result_one.scenario_count == 16
    assert result_one.replicate_count == 2
    assert result_one.index_path.read_bytes() == (second / "scenario_index.tsv").read_bytes()
    assert result_one.manifest_path.read_bytes() == (second / "scenario_manifest.json").read_bytes()

    index_rows = result_one.index_path.read_text(encoding="utf-8").splitlines()
    assert len(index_rows) == 33

    manifest = json.loads(result_one.manifest_path.read_text(encoding="utf-8"))
    assert manifest["total_simulations"] == 32
    assert manifest["scenario_count"] == 16


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("replicates", 0),
        ("depths", ()),
        ("purities", (1.1,)),
        ("clone_counts", (0,)),
        ("topologies", ("unrooted",)),
    ],
)
def test_scenario_grid_rejects_invalid_config(tmp_path, field: str, value: object) -> None:
    kwargs: dict[str, object] = {}
    kwargs[field] = value

    with pytest.raises((ValueError, TypeError)):
        generate_scenario_grid(tmp_path, **kwargs)
