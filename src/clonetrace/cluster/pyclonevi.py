"""Optional subprocess wrapper for the external PyClone-VI executable."""

import csv
import math
import subprocess
from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from clonetrace.cluster.baseline import MutationClusterAssignment
from clonetrace.schemas import MutationRecord

PyCloneDensity = Literal["binomial", "beta-binomial"]

_INPUT_COLUMNS = (
    "mutation_id",
    "sample_id",
    "ref_counts",
    "alt_counts",
    "major_cn",
    "minor_cn",
    "normal_cn",
    "tumour_content",
)
_RESULT_COLUMNS = (
    "mutation_id",
    "sample_id",
    "cluster_id",
    "cellular_prevalence",
    "cellular_prevalence_std",
    "cluster_assignment_prob",
)


class PyCloneVIError(RuntimeError):
    """Raised when the external PyClone-VI command fails or returns invalid data."""


class PyCloneVINotFoundError(PyCloneVIError):
    """Raised when the external PyClone-VI executable is unavailable."""


@dataclass(frozen=True)
class PyCloneObservation:
    """One canonical mutation with explicit PyClone-VI copy-number inputs."""

    mutation: MutationRecord
    major_cn: float
    minor_cn: float
    normal_cn: float
    tumour_content: float

    def __post_init__(self) -> None:
        for field_name in ("major_cn", "minor_cn", "normal_cn", "tumour_content"):
            value = getattr(self, field_name)
            if not math.isfinite(value):
                raise ValueError(f"{field_name} must be finite")
        if self.major_cn < self.minor_cn or self.minor_cn < 0:
            raise ValueError("copy numbers must satisfy major_cn >= minor_cn >= 0")
        if self.major_cn + self.minor_cn <= 0:
            raise ValueError("major_cn and minor_cn must sum to greater than 0")
        if self.normal_cn <= 0:
            raise ValueError("normal_cn must be greater than 0")
        if not 0 <= self.tumour_content <= 1:
            raise ValueError("tumour_content must be between 0 and 1")
        if self.mutation.ref_count + self.mutation.alt_count <= 0:
            raise ValueError("mutations must have positive total read depth")


@dataclass(frozen=True)
class PyCloneVIClusterSummary:
    """Summary of mutations assigned to a PyClone-VI cluster."""

    cluster_id: int
    mutation_count: int
    mean_cellular_prevalence: float


@dataclass(frozen=True)
class PyCloneVIResult:
    """Standardized PyClone-VI assignments and retained run files."""

    sample: str
    method: str
    cluster_count: int
    assignments: tuple[MutationClusterAssignment, ...]
    clusters: tuple[PyCloneVIClusterSummary, ...]
    input_path: Path
    fit_path: Path
    results_path: Path


def _validate_observations(observations: Sequence[PyCloneObservation]) -> str:
    if not observations:
        raise ValueError("observations must not be empty")
    sample = observations[0].mutation.sample
    if any(item.mutation.sample != sample for item in observations):
        raise ValueError("PyClone-VI wrapper accepts one sample per run")
    return sample


def _write_input(path: Path, observations: Sequence[PyCloneObservation]) -> dict[str, int]:
    path.parent.mkdir(parents=True, exist_ok=True)
    mutation_indices: dict[str, int] = {}
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_INPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for index, item in enumerate(observations):
            mutation_id = f"ct-mutation-{index:08d}"
            mutation_indices[mutation_id] = index
            writer.writerow(
                {
                    "mutation_id": mutation_id,
                    "sample_id": item.mutation.sample,
                    "ref_counts": item.mutation.ref_count,
                    "alt_counts": item.mutation.alt_count,
                    "major_cn": item.major_cn,
                    "minor_cn": item.minor_cn,
                    "normal_cn": item.normal_cn,
                    "tumour_content": item.tumour_content,
                }
            )
    return mutation_indices


def _read_results(
    path: Path,
    observations: Sequence[PyCloneObservation],
    mutation_indices: dict[str, int],
) -> tuple[tuple[MutationClusterAssignment, ...], tuple[PyCloneVIClusterSummary, ...]]:
    assignments_by_index: dict[int, MutationClusterAssignment] = {}
    cluster_ccfs: dict[int, list[float]] = defaultdict(list)

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise PyCloneVIError("PyClone-VI results TSV has no header")
        missing_columns = set(_RESULT_COLUMNS) - set(reader.fieldnames)
        if missing_columns:
            raise PyCloneVIError(
                "PyClone-VI results TSV is missing columns: "
                + ", ".join(sorted(missing_columns))
            )

        for line, row in enumerate(reader, start=2):
            if None in row:
                raise PyCloneVIError(f"PyClone-VI results TSV line {line} has extra values")
            mutation_id = row["mutation_id"]
            if mutation_id not in mutation_indices:
                raise PyCloneVIError(
                    f"PyClone-VI returned an unknown mutation ID on line {line}: {mutation_id!r}"
                )
            mutation_index = mutation_indices[mutation_id]
            if mutation_index in assignments_by_index:
                raise PyCloneVIError(
                    f"PyClone-VI returned duplicate results for mutation {mutation_id!r}"
                )
            observation = observations[mutation_index]
            if row["sample_id"] != observation.mutation.sample:
                raise PyCloneVIError(
                    f"PyClone-VI sample ID mismatch for mutation {mutation_id!r}"
                )

            try:
                cluster_id = int(row["cluster_id"])
                cellular_prevalence = float(row["cellular_prevalence"])
                cellular_prevalence_std = float(row["cellular_prevalence_std"])
                assignment_probability = float(row["cluster_assignment_prob"])
            except (TypeError, ValueError) as error:
                raise PyCloneVIError(
                    f"PyClone-VI returned non-numeric result values on line {line}"
                ) from error
            numeric_values = (
                cellular_prevalence,
                cellular_prevalence_std,
                assignment_probability,
            )
            if (
                cluster_id < 0
                or any(not math.isfinite(value) for value in numeric_values)
                or not 0 <= cellular_prevalence <= 1
                or cellular_prevalence_std < 0
                or not 0 <= assignment_probability <= 1
            ):
                raise PyCloneVIError(f"PyClone-VI returned out-of-range values on line {line}")

            assignments_by_index[mutation_index] = MutationClusterAssignment(
                mutation=observation.mutation,
                cluster_id=cluster_id,
                posterior_probability=assignment_probability,
                cluster_probabilities=None,
                cellular_prevalence=cellular_prevalence,
                cellular_prevalence_std=cellular_prevalence_std,
            )
            cluster_ccfs[cluster_id].append(cellular_prevalence)

    if len(assignments_by_index) != len(observations):
        missing_count = len(observations) - len(assignments_by_index)
        raise PyCloneVIError(
            f"PyClone-VI omitted results for {missing_count} input mutation(s)"
        )

    assignments = tuple(assignments_by_index[index] for index in range(len(observations)))
    clusters = tuple(
        PyCloneVIClusterSummary(
            cluster_id=cluster_id,
            mutation_count=len(ccfs),
            mean_cellular_prevalence=sum(ccfs) / len(ccfs),
        )
        for cluster_id, ccfs in sorted(cluster_ccfs.items())
    )
    return assignments, clusters


def _run_command(
    command: Sequence[str],
    *,
    timeout: float | None,
    stage: str,
) -> None:
    try:
        result = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except FileNotFoundError as error:
        raise PyCloneVINotFoundError(
            "PyClone-VI executable was not found. Install PyClone-VI in a separate "
            "environment (for example with Bioconda) and ensure `pyclone-vi` is on PATH."
        ) from error
    except subprocess.TimeoutExpired as error:
        raise PyCloneVIError(f"PyClone-VI {stage} command exceeded its timeout") from error

    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip() or "no error output"
        raise PyCloneVIError(
            f"PyClone-VI {stage} command failed with exit code "
            f"{result.returncode}: {detail}"
        )


def run_pyclonevi(
    observations: Sequence[PyCloneObservation],
    output_dir: str | Path,
    *,
    executable: str | Sequence[str] = "pyclone-vi",
    num_clusters: int = 10,
    density: PyCloneDensity = "binomial",
    num_restarts: int = 1,
    num_threads: int = 1,
    seed: int = 1,
    timeout: float | None = None,
) -> PyCloneVIResult:
    """Run PyClone-VI and normalize its result TSV into CloneTrace assignments.

    Copy number, normal copy number, and tumor content are required per
    mutation so this wrapper does not invent defaults or infer segment overlap.
    Input, HDF5 fit, and TSV result files are retained in ``output_dir``.
    """

    sample = _validate_observations(observations)
    if num_clusters < 2:
        raise ValueError("num_clusters must be at least 2")
    if num_restarts < 1:
        raise ValueError("num_restarts must be at least 1")
    if num_threads < 1:
        raise ValueError("num_threads must be at least 1")
    if seed < 0:
        raise ValueError("seed must be non-negative")
    if density not in {"binomial", "beta-binomial"}:
        raise ValueError("density must be binomial or beta-binomial")
    if timeout is not None and (not math.isfinite(timeout) or timeout <= 0):
        raise ValueError("timeout must be finite and greater than 0")

    executable_parts = [executable] if isinstance(executable, str) else list(executable)
    if not executable_parts or any(not part for part in executable_parts):
        raise ValueError("executable must contain a command")

    work_dir = Path(output_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    input_path = work_dir / "pyclonevi_input.tsv"
    fit_path = work_dir / "pyclonevi_fit.h5"
    results_path = work_dir / "pyclonevi_results.tsv"
    existing = [path.name for path in (input_path, fit_path, results_path) if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing PyClone-VI output file(s): "
            + ", ".join(existing)
        )

    mutation_indices = _write_input(input_path, observations)
    command_prefix = [str(part) for part in executable_parts]
    _run_command(
        [
            *command_prefix,
            "fit",
            "--in-file",
            str(input_path.resolve()),
            "--out-file",
            str(fit_path.resolve()),
            "--num-clusters",
            str(num_clusters),
            "--density",
            density,
            "--num-restarts",
            str(num_restarts),
            "--num-threads",
            str(num_threads),
            "--seed",
            str(seed),
        ],
        timeout=timeout,
        stage="fit",
    )
    if not fit_path.is_file():
        raise PyCloneVIError("PyClone-VI fit command succeeded but did not create its HDF5 output")
    _run_command(
        [
            *command_prefix,
            "write-results-file",
            "--in-file",
            str(fit_path.resolve()),
            "--out-file",
            str(results_path.resolve()),
        ],
        timeout=timeout,
        stage="write-results-file",
    )
    if not results_path.is_file():
        raise PyCloneVIError(
            "PyClone-VI write-results-file command succeeded but did not create its TSV output"
        )
    assignments, clusters = _read_results(results_path, observations, mutation_indices)
    return PyCloneVIResult(
        sample=sample,
        method="pyclonevi",
        cluster_count=len(clusters),
        assignments=assignments,
        clusters=clusters,
        input_path=input_path,
        fit_path=fit_path,
        results_path=results_path,
    )
