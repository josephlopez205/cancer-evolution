"""GEO/SRA processed multi-region or longitudinal data ingestion."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any, Literal

from pydantic import ValidationError

from clonetrace.ingest.gdc import parse_maf, parse_segments
from clonetrace.schemas import (
    ClinicalRecord,
    ClinicalTable,
    CopyNumberSegment,
    CopyNumberTable,
    MutationRecord,
    MutationTable,
)


@dataclass(frozen=True)
class GeoIngestResult:
    accession: str
    license_url: str
    cache_dir: Path
    output_dir: Path
    manifest: Path
    mutation_count: int
    segment_count: int
    clinical_count: int


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_geo_accession(accession: str) -> str:
    accession = accession.strip().upper()
    if not re.fullmatch(r"GSE\d+", accession):
        raise ValueError("GEO accession must be in GSEnnnn format")
    return accession


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_clinical(path: Path) -> list[ClinicalRecord]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        records: list[ClinicalRecord] = []
        for row in reader:
            sample = str(row.get("sample") or row.get("sample_id") or "").strip()
            if not sample:
                continue
            therapy_raw = str(row.get("therapy") or "").strip()
            therapy = [item.strip() for item in therapy_raw.split(",") if item.strip()]

            survival_time: float | None = None
            for name in ("survival_time", "OS_MONTHS", "time"):
                raw = str(row.get(name) or "").strip()
                if raw:
                    try:
                        survival_time = float(raw)
                    except ValueError:
                        pass
                    break

            survival_status: bool | Literal["event", "censored"] | int | None = None
            for name in ("survival_status", "OS_STATUS", "status"):
                raw = str(row.get(name) or "").strip()
                if raw:
                    lower = raw.lower()
                    if raw in {"0", "1"}:
                        survival_status = int(raw)
                    elif "deceased" in lower or "event" in lower or lower.startswith("1"):
                        survival_status = "event"
                    else:
                        survival_status = "censored"
                    break

            try:
                records.append(
                    ClinicalRecord(
                        sample=sample,
                        survival_time=survival_time,
                        survival_status=survival_status,
                        therapy=therapy,
                    )
                )
            except ValidationError:
                continue
    return records


def ingest_geo(
    accession: str,
    processed_dir: Path | str,
    license_url: str,
    cache_dir: Path | str = "data/geo",
    output_dir: Path | str = "data/processed",
) -> GeoIngestResult:
    """Load a verified processed GEO dataset without using raw sequencing reads."""

    accession = validate_geo_accession(accession)
    license_url = license_url.strip()
    if not license_url.startswith(("http://", "https://")):
        raise ValueError("license_url must be an absolute http(s) URL")

    source = Path(processed_dir)
    if not source.is_dir():
        raise FileNotFoundError(f"processed GEO directory not found: {source}")

    cache = Path(cache_dir) / accession
    output = Path(output_dir)
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)

    mutations: list[MutationRecord] = []
    segments: list[CopyNumberSegment] = []
    clinical: list[ClinicalRecord] = []
    cached_files: list[dict[str, Any]] = []

    for source_file in sorted(source.rglob("*")):
        if not source_file.is_file():
            continue
        lower = source_file.name.lower()
        if lower.startswith(".") or lower.endswith((".bam", ".cram", ".sra")):
            continue

        destination = cache / source_file.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and (
            destination.stat().st_size != source_file.stat().st_size
            or _sha256(destination) != _sha256(source_file)
        ):
            destination.unlink()
        if not destination.exists():
            shutil.copy2(source_file, destination)

        cached_files.append(
            {
                "source_path": str(source_file),
                "cached_path": str(destination),
                "size": destination.stat().st_size,
                "sha256": _sha256(destination),
            }
        )

        if lower.endswith(".maf") or lower.endswith(".maf.gz"):
            mutations.extend(parse_maf(destination))
        elif lower.endswith(".seg") or lower.endswith(".txt") or ".seg." in lower:
            segments.extend(parse_segments(destination))
        elif lower.endswith(".tsv") and ("clinical" in lower):
            records = parse_clinical(destination)
            clinical.extend(records)

    mutations = list(
        {
            (item.sample, item.gene, item.chrom, item.pos, item.ref, item.alt): item
            for item in mutations
        }.values()
    )
    segments = list(
        {
            (
                item.sample,
                item.chrom,
                item.start,
                item.end,
                item.total_cn,
                item.minor_cn,
            ): item
            for item in segments
        }.values()
    )

    mutation_table = MutationTable(records=mutations)
    segment_table = CopyNumberTable(records=segments)
    clinical_table = ClinicalTable(records=clinical)

    _write_tsv(
        output / "mutations.tsv",
        [record.model_dump(mode="python") for record in mutation_table.records],
    )
    _write_tsv(
        output / "copy_number_segments.tsv",
        [record.model_dump(mode="python") for record in segment_table.records],
    )
    _write_tsv(
        output / "clinical.tsv",
        [
            {
                "sample": record.sample,
                "survival_time": record.survival_time,
                "survival_status": record.survival_status,
                "therapy": ",".join(record.therapy),
            }
            for record in clinical_table.records
        ],
    )

    manifest: dict[str, Any] = {
        "source": "GEO",
        "accession": accession,
        "license_url": license_url,
        "access_date": date.today().isoformat(),
        "raw_reads_excluded": True,
        "controlled_access_excluded": True,
        "files": cached_files,
        "tables": {
            "mutations": str(output / "mutations.tsv"),
            "copy_number_segments": str(output / "copy_number_segments.tsv"),
            "clinical": str(output / "clinical.tsv"),
        },
    }
    manifest_path = output / "geo_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    return GeoIngestResult(
        accession=accession,
        license_url=license_url,
        cache_dir=cache,
        output_dir=output,
        manifest=manifest_path,
        mutation_count=len(mutation_table.records),
        segment_count=len(segment_table.records),
        clinical_count=len(clinical_table.records),
    )
