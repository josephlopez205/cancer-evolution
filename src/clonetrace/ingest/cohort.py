"""Cohort-level harmonization helpers."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Literal

from clonetrace.schemas import CohortRecord, CohortTable


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _float_or_none(value: str | None) -> float | None:
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _status_or_none(value: str | None) -> bool | Literal["event", "censored"] | int | None:
    if not value:
        return None
    if value in {"0", "1"}:
        return int(value)
    lower = value.lower()
    if lower == "event":
        return "event"
    if lower == "censored":
        return "censored"
    return None


def build_cohort_table(
    input_dir: Path | str,
    output_dir: Path | str,
    source: str,
) -> Path:
    """Build a deterministic cohort table from canonical ingested TSV files."""

    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    mutation_rows = _read_tsv(input_path / "mutations.tsv")
    segment_rows = _read_tsv(input_path / "copy_number_segments.tsv")
    cna_rows = _read_tsv(input_path / "copy_number_calls.tsv")
    clinical_rows = _read_tsv(input_path / "clinical.tsv")

    samples: set[str] = set()
    mutation_counts: dict[str, int] = {}
    segment_counts: dict[str, int] = {}
    clinical_by_sample: dict[str, dict[str, str]] = {}

    for row in mutation_rows:
        sample = (row.get("sample") or "").strip()
        if sample:
            samples.add(sample)
            mutation_counts[sample] = mutation_counts.get(sample, 0) + 1
    for row in segment_rows + cna_rows:
        sample = (row.get("sample") or "").strip()
        if sample:
            samples.add(sample)
            segment_counts[sample] = segment_counts.get(sample, 0) + 1
    for row in clinical_rows:
        sample = (row.get("sample") or "").strip()
        if sample:
            samples.add(sample)
            clinical_by_sample[sample] = dict(row)

    records: list[CohortRecord] = []
    for sample in sorted(samples):
        clinical = clinical_by_sample.get(sample, {})
        therapy_raw = clinical.get("therapy", "")
        therapy = [item.strip() for item in therapy_raw.split(",") if item.strip()]
        records.append(
            CohortRecord(
                sample=sample,
                source=source,
                mutation_count=mutation_counts.get(sample, 0),
                segment_count=segment_counts.get(sample, 0),
                survival_time=_float_or_none(clinical.get("survival_time")),
                survival_status=_status_or_none(clinical.get("survival_status")),
                therapy=therapy,
            )
        )

    CohortTable(records=records)
    path = output_path / "cohort.tsv"
    fieldnames = [
        "sample",
        "source",
        "mutation_count",
        "segment_count",
        "survival_time",
        "survival_status",
        "therapy",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        if not records:
            handle.write("")
            return path
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for record in records:
            output_row: dict[str, Any] = record.model_dump(mode="python")
            output_row["therapy"] = ",".join(record.therapy)
            writer.writerow(output_row)
    return path
