"""GDC open-tier data ingestion."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import urllib.parse
import urllib.request
from collections.abc import Iterable
from dataclasses import dataclass, field
from io import StringIO, TextIOWrapper
from pathlib import Path
from typing import Any, Literal, cast

from pydantic import ValidationError

from clonetrace.ingest.cohort import build_cohort_table
from clonetrace.schemas import (
    ClinicalRecord,
    ClinicalTable,
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
    data_format: str = ""
    access: str = "open"
    sample_by_aliquot: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class GdcIngestResult:
    project: str
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


def _open_filter(project: str, data_type: str, data_format: str | None) -> dict[str, Any]:
    filters: list[dict[str, Any]] = [
        {
            "op": "in",
            "content": {"field": "cases.project.project_id", "value": [project]},
        },
        {"op": "in", "content": {"field": "access", "value": ["open"]}},
        {"op": "in", "content": {"field": "data_type", "value": [data_type]}},
    ]
    if data_format:
        filters.append({"op": "in", "content": {"field": "data_format", "value": [data_format]}})
    return {
        "op": "and",
        "content": filters,
    }


def search_open_files(
    project: str,
    data_type: str,
    data_format: str | None = None,
    include_sample_map: bool = False,
    page_size: int = 100,
) -> list[GdcFile]:
    """List only open-tier files of the requested GDC data type."""

    files: list[GdcFile] = []
    offset = 0
    while True:
        query: dict[str, str | int] = {
            "filters": json.dumps(_open_filter(project, data_type, data_format)),
            "fields": (
                "file_id,file_name,md5sum,file_size,data_category,data_type,data_format,access,"
                "cases.submitter_id,cases.samples.submitter_id,"
                "cases.samples.portions.analytes.aliquots.aliquot_id"
            ),
            "from": offset,
            "size": page_size,
            "sort": "file_name:asc",
        }
        if include_sample_map:
            query["expand"] = "cases.samples.portions.analytes.aliquots"
        params = urllib.parse.urlencode(query)
        payload = _get_json(f"{GDC_API}/files?{params}")
        data = payload["data"]
        hits = data["hits"]
        for hit in hits:
            files.append(
                GdcFile(
                    file_id=str(hit["file_id"]),
                    file_name=str(hit["file_name"]),
                    md5sum=str(hit.get("md5sum") or ""),
                    size=int(hit.get("file_size") or 0),
                    data_category=str(hit.get("data_category") or ""),
                    data_type=str(hit.get("data_type") or ""),
                    data_format=str(hit.get("data_format") or ""),
                    access=str(hit.get("access") or ""),
                    sample_by_aliquot=_aliquot_map(hit) if include_sample_map else {},
                )
            )

        offset += len(hits)
        pagination = data["pagination"]
        if not hits or offset >= int(pagination["total"]):
            break
    return files


def fetch_open_file(file: GdcFile) -> bytes:
    """Fetch one open GDC data file through the API without writing it to disk."""

    if file.access != "open":
        raise ValueError(f"refusing non-open GDC file: {file.file_id}")
    request = urllib.request.Request(f"{GDC_DATA}/{file.file_id}")
    with urllib.request.urlopen(request, timeout=300) as response:
        content = cast(bytes, response.read())
    if file.md5sum:
        actual_md5 = hashlib.md5(content).hexdigest()
        if actual_md5.lower() != file.md5sum.lower():
            raise RuntimeError(
                f"MD5 mismatch for {file.file_name}: expected {file.md5sum}, got {actual_md5}"
            )
    if file.size and len(content) != file.size:
        raise RuntimeError(f"size mismatch for {file.file_name}")
    return content


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
        return _parse_maf_reader(_maf_reader(handle))


def parse_maf_content(content: bytes, file_name: str) -> list[MutationRecord]:
    """Parse MAF bytes returned by the GDC API without creating a local file."""

    if file_name.lower().endswith(".gz"):
        content = gzip.decompress(content)
    text = content.decode("utf-8-sig")
    return _parse_maf_reader(_maf_reader(StringIO(text)))


def _maf_reader(handle: Iterable[str]) -> csv.DictReader[str]:
    data_lines = (line for line in handle if not line.lstrip().startswith("#"))
    return csv.DictReader(data_lines, delimiter="\t")


def _parse_maf_reader(reader: csv.DictReader[str]) -> list[MutationRecord]:
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
        return _parse_segment_reader(csv.DictReader(handle, delimiter="\t"), {})


def parse_segment_content(
    content: bytes,
    file_name: str,
    sample_by_aliquot: dict[str, str],
) -> list[CopyNumberSegment]:
    """Parse allele-specific segments returned by the GDC API in memory."""

    if file_name.lower().endswith(".gz"):
        content = gzip.decompress(content)
    text = content.decode("utf-8-sig")
    return _parse_segment_reader(
        csv.DictReader(StringIO(text), delimiter="\t"),
        sample_by_aliquot,
    )


def _parse_segment_reader(
    reader: csv.DictReader[str],
    sample_by_aliquot: dict[str, str],
) -> list[CopyNumberSegment]:
    segments: list[CopyNumberSegment] = []
    for row in reader:
        aliquot = _first(row, "GDC_Aliquot", "GDC_Aliquot_ID", "Sample", "sample")
        sample = sample_by_aliquot.get(aliquot, aliquot)
        chrom = _first(row, "Chromosome", "chrom")
        start = _int(row, "Start", "start")
        end = _int(row, "End", "end")
        total = _first(row, "Copy_Number", "total_cn")
        minor = _first(row, "Minor_Copy_Number", "minor_cn")
        if not all((sample, chrom, start, end, total, minor)):
            continue
        try:
            total_cn = float(total)
            minor_cn = float(minor)
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


def _aliquot_map(hit: dict[str, Any]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    cases = hit.get("cases", [])
    for case in cases if isinstance(cases, list) else []:
        samples = case.get("samples", []) if isinstance(case, dict) else []
        for sample in samples if isinstance(samples, list) else []:
            sample_id = str(sample.get("submitter_id", "")) if isinstance(sample, dict) else ""
            portions = sample.get("portions", []) if isinstance(sample, dict) else []
            for portion in portions if isinstance(portions, list) else []:
                analytes = portion.get("analytes", []) if isinstance(portion, dict) else []
                for analyte in analytes if isinstance(analytes, list) else []:
                    aliquots = analyte.get("aliquots", []) if isinstance(analyte, dict) else []
                    for aliquot in aliquots if isinstance(aliquots, list) else []:
                        if isinstance(aliquot, dict) and aliquot.get("aliquot_id") and sample_id:
                            mapping[str(aliquot["aliquot_id"])] = sample_id
    return mapping


def _source_file_manifest(file: GdcFile, content: bytes) -> dict[str, Any]:
    return {
        "file_id": file.file_id,
        "file_name": file.file_name,
        "data_type": file.data_type,
        "data_format": file.data_format,
        "access": file.access,
        "md5sum": file.md5sum,
        "sha256": hashlib.sha256(content).hexdigest(),
        "size": len(content),
    }


def fetch_open_cases(project: str, page_size: int = 1000) -> list[dict[str, Any]]:
    """Fetch public clinical fields from the GDC Cases API."""

    fields = (
        "submitter_id,samples.submitter_id,demographic.vital_status,"
        "demographic.days_to_death,diagnoses.days_to_last_follow_up,"
        "diagnoses.treatments.treatment_type"
    )
    filters = {
        "op": "in",
        "content": {"field": "project.project_id", "value": [project]},
    }
    cases: list[dict[str, Any]] = []
    offset = 0
    while True:
        payload = _post_json(
            f"{GDC_API}/cases",
            {
                "filters": filters,
                "fields": fields,
                "from": offset,
                "size": page_size,
                "sort": "submitter_id:asc",
            },
        )
        data = payload["data"]
        hits = data["hits"]
        cases.extend(hit for hit in hits if isinstance(hit, dict))
        offset += len(hits)
        if not hits or offset >= int(data["pagination"]["total"]):
            break
    return cases


def parse_clinical_cases(cases: list[dict[str, Any]]) -> list[ClinicalRecord]:
    records: list[ClinicalRecord] = []
    for case in cases:
        demographic = case.get("demographic") or {}
        vital_status = str(demographic.get("vital_status") or "").lower()
        status: Literal["event", "censored"] | None = (
            "event" if vital_status == "dead" else "censored" if vital_status == "alive" else None
        )
        survival_time = demographic.get("days_to_death") if status == "event" else None
        diagnoses = case.get("diagnoses") or []
        therapies: list[str] = []
        if survival_time is None:
            for diagnosis in diagnoses:
                if survival_time is None:
                    survival_time = diagnosis.get("days_to_last_follow_up")
        for diagnosis in diagnoses:
            for treatment in diagnosis.get("treatments") or []:
                treatment_type = treatment.get("treatment_type")
                if treatment_type:
                    therapies.append(str(treatment_type))

        sample_ids = [
            str(sample["submitter_id"])
            for sample in case.get("samples") or []
            if isinstance(sample, dict) and sample.get("submitter_id")
        ]
        if not sample_ids and case.get("submitter_id"):
            sample_ids = [str(case["submitter_id"])]
        for sample_id in sample_ids:
            try:
                records.append(
                    ClinicalRecord(
                        sample=sample_id,
                        survival_time=survival_time,
                        survival_status=status,
                        therapy=therapies,
                    )
                )
            except ValidationError:
                continue
    return records


def _write_tsv(
    path: Path,
    rows: Iterable[dict[str, Any]],
    fieldnames: list[str],
) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def ingest_gdc(
    project: str,
    output_dir: Path | str = "data/processed",
) -> GdcIngestResult:
    """Use GDC APIs to build canonical tables without storing source files locally."""

    project = project.strip().upper()
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    maf_files = search_open_files(project, "Masked Somatic Mutation", "MAF")
    segment_files = search_open_files(
        project,
        "Allele-specific Copy Number Segment",
        "TXT",
        include_sample_map=True,
    )
    mutations: list[MutationRecord] = []
    source_files: list[dict[str, Any]] = []
    for file in maf_files:
        content = fetch_open_file(file)
        source_files.append(_source_file_manifest(file, content))
        mutations.extend(parse_maf_content(content, file.file_name))

    segments: list[CopyNumberSegment] = []
    for file in segment_files:
        content = fetch_open_file(file)
        source_files.append(_source_file_manifest(file, content))
        segments.extend(parse_segment_content(content, file.file_name, file.sample_by_aliquot))

    case_rows = fetch_open_cases(project)
    clinical = parse_clinical_cases(case_rows)

    mutation_table = MutationTable(records=mutations)
    segment_table = CopyNumberTable(records=segments)
    clinical_table = ClinicalTable(records=clinical)
    _write_tsv(
        output / "mutations.tsv",
        [record.model_dump(mode="python") for record in mutation_table.records],
        list(MutationRecord.model_fields),
    )
    _write_tsv(
        output / "copy_number_segments.tsv",
        [record.model_dump(mode="python") for record in segment_table.records],
        list(CopyNumberSegment.model_fields),
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
        ["sample", "survival_time", "survival_status", "therapy"],
    )
    cohort_path = build_cohort_table(output, output, source="gdc")

    manifest: dict[str, Any] = {
        "project": project,
        "source": "GDC",
        "access": "open",
        "raw_files_persisted": False,
        "clinical_api": {
            "endpoint": f"{GDC_API}/cases",
            "record_count": len(clinical_table.records),
            "sha256": hashlib.sha256(
                json.dumps(case_rows, sort_keys=True).encode("utf-8")
            ).hexdigest(),
            "time_unit": "days",
        },
        "files": source_files,
        "tables": {
            "mutations": str(output / "mutations.tsv"),
            "copy_number_segments": str(output / "copy_number_segments.tsv"),
            "clinical": str(output / "clinical.tsv"),
            "cohort": str(cohort_path),
        },
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    return GdcIngestResult(
        project=project,
        output_dir=output,
        manifest=manifest_path,
        maf_count=len(maf_files),
        segment_count=len(segment_files),
        clinical_count=len(clinical_table.records),
        mutation_count=len(mutation_table.records),
        segment_record_count=len(segment_table.records),
    )
