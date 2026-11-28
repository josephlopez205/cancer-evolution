"""Deterministic simulation of clonal mutations and sequencing read counts."""

import csv
import hashlib
import json
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

from clonetrace.schemas import CopyNumberSegment, MutationRecord

Topology = Literal["linear", "branching"]


@dataclass(frozen=True)
class SimulationConfig:
    """Configuration for one simulated tumor sample."""

    seed: int = 1
    depth: int = 100
    purity: float = 0.8
    clone_count: int = 3
    topology: Topology = "branching"
    mutation_count: int = 40
    cn_event_probability: float = 0.25
    sample: str = "SIMULATED-1"


@dataclass(frozen=True)
class SimulationResult:
    """Paths and summary values produced by a simulation."""

    output_dir: Path
    mutations_path: Path
    segments_path: Path
    ground_truth_path: Path
    manifest_path: Path
    mutation_count: int
    clone_count: int


def _validate_config(config: SimulationConfig) -> None:
    if config.depth < 1:
        raise ValueError("depth must be at least 1")
    if not 0 < config.purity <= 1:
        raise ValueError("purity must be greater than 0 and at most 1")
    if config.clone_count < 1:
        raise ValueError("clone_count must be at least 1")
    if config.mutation_count < 1:
        raise ValueError("mutation_count must be at least 1")
    if config.topology not in {"linear", "branching"}:
        raise ValueError("topology must be 'linear' or 'branching'")
    if not 0 <= config.cn_event_probability <= 1:
        raise ValueError("cn_event_probability must be between 0 and 1")
    if not config.sample.strip():
        raise ValueError("sample must not be empty")


def _draw_binomial(rng: random.Random, trials: int, probability: float) -> int:
    if probability <= 0:
        return 0
    if probability >= 1:
        return trials
    successes = 0
    for _ in range(trials):
        if rng.random() < probability:
            successes += 1
    return successes


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def simulate(
    output_dir: str | Path,
    *,
    seed: int = 1,
    depth: int = 100,
    purity: float = 0.8,
    clone_count: int = 3,
    topology: Topology = "branching",
    mutation_count: int = 40,
    cn_event_probability: float = 0.25,
    sample: str = "SIMULATED-1",
) -> SimulationResult:
    """Generate one deterministic simulated sample and its ground truth."""

    config = SimulationConfig(
        seed=seed,
        depth=depth,
        purity=purity,
        clone_count=clone_count,
        topology=topology,
        mutation_count=mutation_count,
        cn_event_probability=cn_event_probability,
        sample=sample,
    )
    _validate_config(config)
    rng = random.Random(seed)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    parents: list[int] = [0]
    for clone_id in range(1, clone_count):
        if topology == "linear":
            parents.append(clone_id - 1)
        else:
            parents.append(rng.randrange(clone_id))

    weights = [rng.uniform(0.2, 1.0) for _ in range(clone_count)]
    weight_total = sum(weights)
    clone_fractions = [weight / weight_total for weight in weights]

    cn_events: list[dict[str, Any]] = []
    segments: list[CopyNumberSegment] = []
    for clone_id in range(1, clone_count):
        if rng.random() >= cn_event_probability:
            continue
        chromosome = str(rng.randint(1, 22))
        start = rng.randint(1_000_000, 20_000_000)
        end = start + rng.randint(10_000_000, 100_000_000)
        total_cn = float(rng.choice([3, 4]))
        minor_cn = float(rng.choice([1, 2]) if total_cn >= 3 else 1)
        if minor_cn > total_cn:
            minor_cn = total_cn
        segment = CopyNumberSegment(
            sample=sample,
            chrom=chromosome,
            start=start,
            end=end,
            total_cn=total_cn,
            minor_cn=minor_cn,
        )
        segments.append(segment)
        cn_events.append(
            {
                "clone": clone_id,
                "chrom": chromosome,
                "start": start,
                "end": end,
                "total_cn": total_cn,
                "minor_cn": minor_cn,
                "parent": parents[clone_id],
            }
        )

    mutation_assignments: list[dict[str, Any]] = []
    mutation_rows: list[dict[str, Any]] = []
    gene_names = ("TP53", "PIK3CA", "KRAS", "GATA3", "PTEN", "ERBB2")
    alleles = (("A", "C"), ("C", "T"), ("G", "A"), ("T", "G"))

    for mutation_id in range(mutation_count):
        clone_id = rng.choices(
            range(clone_count),
            weights=clone_fractions,
            k=1,
        )[0]
        ccf = clone_fractions[clone_id]
        copy_number = 2.0
        for event in cn_events:
            if event["clone"] == clone_id:
                copy_number = float(event["total_cn"])
                break

        denominator = 2 * (1 - purity) + purity * copy_number
        alt_probability = purity * ccf / denominator
        alt_probability = min(max(alt_probability, 0.0), 1.0)
        observed_depth = max(1, int(rng.gauss(depth, depth * 0.05)))
        alt_count = _draw_binomial(rng, observed_depth, alt_probability)
        ref_count = observed_depth - alt_count
        ref, alt = rng.choice(alleles)

        record = MutationRecord(
            sample=sample,
            gene=gene_names[mutation_id % len(gene_names)],
            chrom=str(rng.randint(1, 22)),
            pos=rng.randint(1_000_000, 200_000_000),
            ref=ref,
            alt=alt,
            ref_count=ref_count,
            alt_count=alt_count,
        )
        mutation_rows.append(record.model_dump())
        mutation_assignments.append(
            {
                "mutation_id": mutation_id,
                "clone": clone_id,
                "gene": record.gene,
                "chrom": record.chrom,
                "pos": record.pos,
                "true_ccf": ccf,
            }
        )

    mutations_path = output_path / "mutations.tsv"
    segments_path = output_path / "copy_number_segments.tsv"
    ground_truth_path = output_path / "ground_truth.json"
    manifest_path = output_path / "manifest.json"

    mutation_fields = [
        "sample",
        "gene",
        "chrom",
        "pos",
        "ref",
        "alt",
        "ref_count",
        "alt_count",
    ]
    _write_tsv(mutations_path, mutation_fields, mutation_rows)
    segment_fields = ["sample", "chrom", "start", "end", "total_cn", "minor_cn"]
    _write_tsv(
        segments_path,
        segment_fields,
        [segment.model_dump() for segment in segments],
    )

    ground_truth = {
        "format_version": "1",
        "config": asdict(config),
        "clone_tree": {
            "parents": parents,
            "clone_fractions": clone_fractions,
        },
        "mutation_assignments": mutation_assignments,
        "copy_number_events": cn_events,
    }
    ground_truth_path.write_text(
        json.dumps(ground_truth, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    manifest = {
        "format_version": "1",
        "config": asdict(config),
        "files": {
            mutations_path.name: _sha256(mutations_path),
            segments_path.name: _sha256(segments_path),
            ground_truth_path.name: _sha256(ground_truth_path),
        },
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return SimulationResult(
        output_dir=output_path,
        mutations_path=mutations_path,
        segments_path=segments_path,
        ground_truth_path=ground_truth_path,
        manifest_path=manifest_path,
        mutation_count=mutation_count,
        clone_count=clone_count,
    )
