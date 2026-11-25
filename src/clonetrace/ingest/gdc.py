"""GDC open-tier data ingestion."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import shutil
import urllib.parse
import urllib.request
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from io import TextIOWrapper
from pathlib import Path
from typing import Any, cast

from pydantic import ValidationError

from clonetrace.schemas import (
    CopyNumberSegment,
    CopyNumberTable,
    MutationRecord,
    MutationTable,
)

GDC_API = "https://api.gdc.cancer.gov"
GDC_DATA = f"{GDC_API}/data"


@dataclass(frozen=True)
class GdcFile:
    file_id: str
    file_name: str
    md5sum: str
    size: int
    data_category: str
    data_type: str


@dataclass(frozen=True)
class GdcIngestResult:
    project: str
    cache_dir: Path
    output_dir: Path
    manifest: Path
    maf_count: int
    segment_count: int
    clinical_count: int
    mutation_count: int
    segment_record_count: int


def _post_json(url: str, body: dict[str, Any]) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        data=json.dumps(body).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return cast(dict[str, Any], json.loads(response.read()))


def _get_json(url: str) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return cast(dict[str, Any], json.loads(response.read()))


def _open_filter(project: str, data_category: str) -> dict[str, Any]:
    return {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {"field": "cases.project.project_id", "value": [project]},
            },
            {"op": "in", "content": {"field": "access", "value": ["open"]}},
            {
                "op": "in",
                "content": {"field": "data_category", "value": [data_category]},
            },
        ],
    }


def search_open_files(
    project: str,
    data_category: str,
    page_size: int = 100,
) -> list[GdcFile]:
    """List open-tier GDC files for a project and data category."""

    files: list[GdcFile] = []
    offset = 0
    while True:
        params = urllib.parse.urlencode(
            {
                "filters": json.dumps(_open_filter(project, data_category)),
                "fields": "file_id,file_name,md5sum,size,data_category,data_type",
                "from": offset,
                "size": page_size,
                "sort": "file_name:asc",
            }
        )
        payload = _get_json(f"{GDC_API}/files?{params}")
        data = payload["data"]
        hits = data["hits"]
        for hit in hits:
            files.append(
                GdcFile(
                    file_id=str(hit["file_id"]),
                    file_name=str(hit["file_name"]),
                    md5sum=str(hit.get("md5sum") or ""),
                    size=int(hit.get("size") or 0),
                    data_category=str(hit.get("data_category") or ""),
                    data_type=str(hit.get("data_type") or ""),
                )
            )

        offset += len(hits)
        pagination = data["pagination"]
        if not hits or offset >= int(pagination["total"]):
            break
    return files


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_open_file(file: GdcFile, cache_dir: Path, force: bool = False) -> Path:
    """Download one open GDC file to cache and verify source checksums."""

    cache_dir.mkdir(parents=True, exist_ok=True)
    destination = cache_dir / file.file_id / file.file_name
    destination.parent.mkdir(parents=True, exist_ok=True)

    if destination.exists() and not force:
        if file.size and destination.stat().st_size != file.size:
            raise RuntimeError(f"cached file size mismatch: {destination}")
    else:
        temporary = destination.with_suffix(destination.suffix + ".part")
        request = urllib.request.Request(f"{GDC_DATA}/{file.file_id}")
        with urllib.request.urlopen(request, timeout=300) as response:
            with temporary.open("wb") as handle:
                shutil.copyfileobj(response, handle)
        temporary.replace(destination)

    if file.md5sum:
        actual_md5 = _md5(destination)
        if actual_md5.lower() != file.md5sum.lower():
            raise RuntimeError(
                f"MD5 mismatch for {file.file_name}: expected {file.md5sum}, got {actual_md5}"
            )
    if file.size and destination.stat().st_size != file.size:
        raise RuntimeError(f"size mismatch for {file.file_name}")
    return destination


def _open_text(path: Path) -> TextIOWrapper:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def _first(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = row.get(name)
        if value:
            return value
    return ""


def _int(row: dict[str, str], *names: str) -> int:
    value = _first(row, *names)
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def parse_maf(path: Path) -> list[MutationRecord]:
    """Parse a GDC masked MAF into canonical mutation records."""

    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        mutations: list[MutationRecord] = []
        for row in reader:
            if not row:
                continue
            sample = _first(row, "Tumor_Sample_Barcode", "sample")
            gene = _first(row, "Hugo_Symbol", "gene")
            chrom = _first(row, "Chromosome", "chrom")
            pos = _int(row, "Start_Position", "pos")
            ref = _first(row, "Reference_Allele", "ref")
            alt = _first(row, "Tumor_Seq_Allele2", "Tumor_Seq_Allele1", "alt")
            ref_count = _int(row, "t_ref_count", "ref_count")
            alt_count = _int(row, "t_alt_count", "alt_count")
            if not all((sample, gene, chrom, pos, ref, alt)):
                continue
            try:
                mutations.append(
                    MutationRecord(
                        sample=sample,
                        gene=gene,
                        chrom=chrom,
                        pos=pos,
                        ref=ref,
                        alt=alt,
                        ref_count=ref_count,
                        alt_count=alt_count,
                    )
                )
            except ValidationError:
                continue
    return mutations


def parse_segments(path: Path) -> list[CopyNumberSegment]:
    """Parse GDC segment files into canonical copy-number records."""

    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        segments: list[CopyNumberSegment] = []
        for row in reader:
            sample = _first(row, "GDC_Aliquot_ID", "Sample", "sample", "Sample_name")
            chrom = _first(row, "Chromosome", "chrom")
            start = _int(row, "Start", "start")
            end = _int(row, "End", "end")
            total = _first(
                row,
                "Copy_Number",
                "total_cn",
                "Segment_Copy_Number",
            )
            minor = _first(row, "Minor_Copy_Number", "minor_cn")
            if not all((sample, chrom, start, end)):
                continue
            try:
                total_cn = float(total)
            except ValueError:
                continue
            try:
                minor_cn = float(minor) if minor else max(0.0, total_cn / 2.0)
            except ValueError:
                continue
            try:
                segments.append(
                    CopyNumberSegment(
                        sample=sample,
                        chrom=chrom,
                        start=start,
                        end=end,
                        total_cn=total_cn,
                        minor_cn=minor_cn,
                    )
                )
            except ValidationError:
                continue
    return segments


def _write_tsv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def ingest_gdc(
    project: str,
    cache_dir: Path | str = "data/gdc",
    output_dir: Path | str = "data/processed",
    force: bool = False,
) -> GdcIngestResult:
    """Download open GDC MAF, segment, and clinical files and build canonical tables."""

    project = project.strip().upper()
    cache = Path(cache_dir)
    output = Path(output_dir)
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)

    maf_files = search_open_files(project, "Simple Nucleotide Variation")
    maf_files = [
        file
        for file in maf_files
        if file.file_name.lower().endswith(".maf") or file.file_name.lower().endswith(".maf.gz")
    ]
    segment_files = search_open_files(project, "Copy Number Variation")
    clinical_files = search_open_files(project, "Clinical")

    downloads: list[tuple[GdcFile, Path]] = []
    for file in maf_files + segment_files + clinical_files:
        downloads.append((file, download_open_file(file, cache, force=force)))

    mutations: list[MutationRecord] = []
    for file in maf_files:
        path = next(path for selected, path in downloads if selected.file_id == file.file_id)
        mutations.extend(parse_maf(path))

    segments: list[CopyNumberSegment] = []
    for file in segment_files:
        path = next(path for selected, path in downloads if selected.file_id == file.file_id)
        if (
            path.name.lower().endswith(".seg")
            or path.name.lower().endswith(".txt")
            or ".seg." in path.name.lower()
        ):
            segments.extend(parse_segments(path))

    mutation_table = MutationTable(records=mutations)
    segment_table = CopyNumberTable(records=segments)
    _write_tsv(
        output / "mutations.tsv",
        [record.model_dump(mode="python") for record in mutation_table.records],
    )
    _write_tsv(
        output / "copy_number_segments.tsv",
        [record.model_dump(mode="python") for record in segment_table.records],
    )

    manifest: dict[str, Any] = {
        "project": project,
        "source": "GDC",
        "access": "open",
        "files": [],
        "tables": {
            "mutations": str(output / "mutations.tsv"),
            "copy_number_segments": str(output / "copy_number_segments.tsv"),
        },
    }
    for file, path in downloads:
        manifest["files"].append(
            {
                **asdict(file),
                "cached_path": str(path),
                "md5": _md5(path),
                "sha256": _sha256(path),
                "size": path.stat().st_size,
            }
        )
    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    return GdcIngestResult(
        project=project,
        cache_dir=cache,
        output_dir=output,
        manifest=manifest_path,
        maf_count=len(maf_files),
        segment_count=len(segment_files),
        clinical_count=len(clinical_files),
        mutation_count=len(mutation_table.records),
        segment_record_count=len(segment_table.records),
    )
