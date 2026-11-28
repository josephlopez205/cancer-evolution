"""Reproducible scenario grids for simulation benchmarking."""

import csv
import json
from dataclasses import asdict, dataclass
from itertools import product
from pathlib import Path
from typing import Literal

from clonetrace.simulate.simulator import SimulationResult, simulate

Topology = Literal["linear", "branching"]

DEFAULT_DEPTHS = (30, 60, 100, 200)
DEFAULT_PURITIES = (0.3, 0.5, 0.8)
DEFAULT_CLONE_COUNTS = tuple(range(1, 7))
DEFAULT_TOPOLOGIES: tuple[Topology, ...] = ("linear", "branching")


@dataclass(frozen=True)
class ScenarioGridConfig:
    """Configuration for a complete simulation scenario grid."""

    seed: int = 1
    replicates: int = 20
    mutation_count: int = 40
    cn_event_probability: float = 0.25
    depths: tuple[int, ...] = DEFAULT_DEPTHS
    purities: tuple[float, ...] = DEFAULT_PURITIES
    clone_counts: tuple[int, ...] = DEFAULT_CLONE_COUNTS
    topologies: tuple[Topology, ...] = DEFAULT_TOPOLOGIES


@dataclass(frozen=True)
class ScenarioGridResult:
    """Paths and dimensions produced by a scenario grid."""

    output_dir: Path
    index_path: Path
    manifest_path: Path
    scenario_count: int
    replicate_count: int


def _validate_grid_config(config: ScenarioGridConfig) -> None:
    if config.seed < 0:
        raise ValueError("seed must be non-negative")
    if config.replicates < 1:
        raise ValueError("replicates must be at least 1")
    if config.mutation_count < 1:
        raise ValueError("mutation_count must be at least 1")
    if not 0 <= config.cn_event_probability <= 1:
        raise ValueError("cn_event_probability must be between 0 and 1")
    if not config.depths or any(depth < 1 for depth in config.depths):
        raise ValueError("depths must contain positive values")
    if not config.purities or any(purity <= 0 or purity > 1 for purity in config.purities):
        raise ValueError("purities must contain values greater than 0 and at most 1")
    if not config.clone_counts or any(clone_count < 1 for clone_count in config.clone_counts):
        raise ValueError("clone_counts must contain positive values")
    if not config.topologies:
        raise ValueError("topologies must not be empty")
    if any(topology not in {"linear", "branching"} for topology in config.topologies):
        raise ValueError("topologies must contain only linear or branching")


def generate_scenario_grid(
    output_dir: str | Path,
    *,
    seed: int = 1,
    replicates: int = 20,
    mutation_count: int = 40,
    cn_event_probability: float = 0.25,
    depths: tuple[int, ...] = DEFAULT_DEPTHS,
    purities: tuple[float, ...] = DEFAULT_PURITIES,
    clone_counts: tuple[int, ...] = DEFAULT_CLONE_COUNTS,
    topologies: tuple[Topology, ...] = DEFAULT_TOPOLOGIES,
) -> ScenarioGridResult:
    """Generate every scenario combination and 20 or more seeded replicates."""

    config = ScenarioGridConfig(
        seed=seed,
        replicates=replicates,
        mutation_count=mutation_count,
        cn_event_probability=cn_event_probability,
        depths=depths,
        purities=purities,
        clone_counts=clone_counts,
        topologies=topologies,
    )
    _validate_grid_config(config)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    index_fields = [
        "scenario_id",
        "replicate",
        "seed",
        "depth",
        "purity",
        "clones",
        "topology",
        "mutations",
        "mutations_path",
        "segments_path",
        "ground_truth_path",
        "manifest_path",
    ]
    rows: list[dict[str, object]] = []
    scenario_index = 0

    combinations = product(
        config.depths,
        config.purities,
        config.clone_counts,
        config.topologies,
    )
    for depth, purity, clone_count, topology in combinations:
        scenario_index += 1
        scenario_id = f"scenario-{scenario_index:04d}"
        for replicate in range(1, config.replicates + 1):
            replicate_seed = (
                config.seed
                + ((scenario_index - 1) * config.replicates)
                + replicate
                - 1
            )
            replicate_dir = (
                output_path
                / scenario_id
                / f"depth-{depth}"
                / f"purity-{purity:g}"
                / f"clones-{clone_count}"
                / topology
                / f"replicate-{replicate:03d}"
            )
            sample = f"{scenario_id}-r{replicate:03d}"
            result: SimulationResult = simulate(
                replicate_dir,
                seed=replicate_seed,
                depth=depth,
                purity=purity,
                clone_count=clone_count,
                topology=topology,
                mutation_count=config.mutation_count,
                cn_event_probability=config.cn_event_probability,
                sample=sample,
            )
            rows.append(
                {
                    "scenario_id": scenario_id,
                    "replicate": replicate,
                    "seed": replicate_seed,
                    "depth": depth,
                    "purity": purity,
                    "clones": clone_count,
                    "topology": topology,
                    "mutations": config.mutation_count,
                    "mutations_path": str(result.mutations_path.relative_to(output_path)),
                    "segments_path": str(result.segments_path.relative_to(output_path)),
                    "ground_truth_path": str(result.ground_truth_path.relative_to(output_path)),
                    "manifest_path": str(result.manifest_path.relative_to(output_path)),
                }
            )

    index_path = output_path / "scenario_index.tsv"
    with index_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=index_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    manifest_path = output_path / "scenario_manifest.json"
    manifest = {
        "format_version": "1",
        "config": asdict(config),
        "scenario_count": scenario_index,
        "replicate_count": config.replicates,
        "total_simulations": scenario_index * config.replicates,
        "index": index_path.name,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return ScenarioGridResult(
        output_dir=output_path,
        index_path=index_path,
        manifest_path=manifest_path,
        scenario_count=scenario_index,
        replicate_count=config.replicates,
    )
