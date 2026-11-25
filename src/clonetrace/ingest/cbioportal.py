"""cBioPortal public API ingestion."""

from __future__ import annotations

import hashlib
import json
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, cast

from pydantic import ValidationError

from clonetrace.schemas import (
    ClinicalRecord,
    ClinicalTable,
    MutationRecord,
    MutationTable,
)

CBIOPORTAL_API = "https://www.cbioportal.org/api"


@dataclass(frozen=True)
class CbioPortalIngestResult:
    study: str
    cache_dir: Path
    output_dir: Path
    manifest: Path
    mutation_count: int
    clinical_count: int
    cna_count: int


def _get_json(url: str) -> dict[str, Any] | list[Any]:
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return cast(dict[str, Any] | list[Any], json.loads(response.read()))


def _post_json(url: str, body: dict[str, Any]) -> dict[str, Any] | list[Any]:
    request = urllib.request.Request(
        url,
        data=json.dumps(body).encode("utf-8"),
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=300) as response:
        return cast(dict[str, Any] | list[Any], json.loads(response.read()))


def _unwrap(payload: dict[str, Any] | list[Any]) -> list[Any]:
    if isinstance(payload, list):
        return payload
    value = payload.get("value")
    if isinstance(value, list):
        return value
    return []


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _cache_json(
    cache_dir: Path,
    name: str,
    payload: dict[str, Any] | list[Any],
    force: bool = False,
) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{name}.json"
    if path.exists() and not force:
        return path
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return path


def get_study(study: str) -> dict[str, Any]:
    payload = _get_json(f"{CBIOPORTAL_API}/studies/{urllib.parse.quote(study)}")
    if not isinstance(payload, dict):
        raise RuntimeError(f"unexpected cBioPortal study response for {study}")
    return payload


def get_molecular_profiles(study: str) -> list[Any]:
    payload = _get_json(f"{CBIOPORTAL_API}/studies/{urllib.parse.quote(study)}/molecular-profiles")
    return _unwrap(payload)


def get_mutations(study: str, molecular_profile: str) -> list[Any]:
    payload = _post_json(
        f"{CBIOPORTAL_API}/mutations/fetch",
        {
            "molecularProfileId": molecular_profile,
            "sampleListId": f"{study}_all",
            "projection": "DETAILED",
        },
    )
    return _unwrap(payload)


def get_discrete_cna(study: str, molecular_profile: str) -> list[Any]:
    payload = _post_json(
        f"{CBIOPORTAL_API}/discrete-copy-number/fetch",
        {
            "molecularProfileIds": [molecular_profile],
            "sampleIds": [],
            "discreteCopyNumberEventType": "ALL",
            "projection": "DETAILED",
        },
    )
    return _unwrap(payload)


def get_clinical_data(study: str, page_size: int = 1000) -> list[Any]:
    rows: list[Any] = []
    page = 0
    while True:
        params = urllib.parse.urlencode(
            {
                "direction": "ASC",
                "pageNumber": page,
                "pageSize": page_size,
                "projection": "DETAILED",
            }
        )
        payload = _get_json(
            f"{CBIOPORTAL_API}/studies/{urllib.parse.quote(study)}/clinical-data?{params}"
        )
        batch = _unwrap(payload)
        if not batch:
            break
        rows.extend(batch)
        if len(batch) < page_size:
            break
        page += 1
    return rows


def _first(row: dict[str, Any], *names: str) -> Any:
    for name in names:
        value = row.get(name)
        if value is not None:
            return value
    return None


def parse_mutations(rows: list[Any]) -> list[MutationRecord]:
    mutations: list[MutationRecord] = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        gene = row.get("gene")
        gene_symbol = gene.get("hugoGeneSymbol") if isinstance(gene, dict) else None
        chromosome = _first(row, "chromosome", "chr")
        start = _first(row, "startPosition", "start")
        ref = _first(row, "referenceAllele", "ref")
        alt = _first(row, "variantAllele", "tumorSeqAllele", "alt")
        sample = _first(row, "sampleId", "sample")
        ref_count = _first(row, "normalRefCount", "t_ref_count", "refCount")
        alt_count = _first(row, "tumorAltCount", "t_alt_count", "altCount")
        try:
            if not all((sample, gene_symbol, chromosome, start, ref, alt)):
                continue
            mutations.append(
                MutationRecord(
                    sample=str(sample),
                    gene=str(gene_symbol),
                    chrom=str(chromosome),
                    pos=int(start),
                    ref=str(ref),
                    alt=str(alt),
                    ref_count=int(ref_count or 0),
                    alt_count=int(alt_count or 0),
                )
            )
        except (TypeError, ValueError, ValidationError):
            continue
    return mutations


def parse_clinical(rows: list[Any]) -> list[ClinicalRecord]:
    by_sample: dict[str, dict[str, str]] = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        sample = _first(row, "sampleId", "patientId")
        attribute = _first(row, "clinicalAttributeId")
        value = _first(row, "value")
        if sample and attribute and value is not None:
            by_sample.setdefault(str(sample), {})[str(attribute)] = str(value)

    records: list[ClinicalRecord] = []
    for sample, attributes in sorted(by_sample.items()):
        os_status = attributes.get("OS_STATUS")
        status: bool | Literal["event", "censored"] | int | None = None
        if os_status:
            status = (
                "event"
                if os_status.lower().startswith("1") or "deceased" in os_status.lower()
                else "censored"
            )
        try:
            os_months = float(attributes["OS_MONTHS"]) if "OS_MONTHS" in attributes else None
        except ValueError:
            os_months = None

        therapy = [
            name
            for name in (
                "CHEMOTHERAPY",
                "HORMONE_THERAPY",
                "RADIO_THERAPY",
            )
            if attributes.get(name, "").strip().lower() not in {"", "no", "none", "not recorded"}
        ]
        therapy_labels = [
            attributes[name] if attributes.get(name) else name.lower().replace("_", " ")
            for name in therapy
        ]
        try:
            records.append(
                ClinicalRecord(
                    sample=sample,
                    survival_time=os_months,
                    survival_status=status,
                    therapy=therapy_labels,
                )
            )
        except ValidationError:
            continue
    return records


def parse_discrete_cna(rows: list[Any]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        gene = row.get("gene")
        records.append(
            {
                "sample": str(row.get("sampleId", "")),
                "gene": str(gene.get("hugoGeneSymbol", "") if isinstance(gene, dict) else ""),
                "entrez_gene_id": row.get("entrezGeneId"),
                "cna_value": row.get("value"),
                "molecular_profile_id": row.get("molecularProfileId"),
            }
        )
    return [record for record in records if record["sample"] and record["gene"]]


def _write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        if not rows:
            path.write_text("", encoding="utf-8")
            return
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def ingest_cbioportal(
    study: str,
    cache_dir: Path | str = "data/cbioportal",
    output_dir: Path | str = "data/processed",
    force: bool = False,
) -> CbioPortalIngestResult:
    """Download and harmonize mutations, CNA calls, and clinical data for a public study."""

    study = study.strip()
    cache = Path(cache_dir) / study
    output = Path(output_dir)
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)

    study_metadata = get_study(study)
    profiles = get_molecular_profiles(study)
    mutation_profile = next(
        (
            profile["molecularProfileId"]
            for profile in profiles
            if isinstance(profile, dict)
            and profile.get("molecularAlterationType") == "MUTATION_EXTENDED"
        ),
        f"{study}_mutations",
    )
    cna_profile = next(
        (
            profile["molecularProfileId"]
            for profile in profiles
            if isinstance(profile, dict)
            and profile.get("molecularAlterationType") == "COPY_NUMBER_ALTERATION"
        ),
        f"{study}_cna",
    )

    mutation_rows = get_mutations(study, mutation_profile)
    cna_rows = get_discrete_cna(study, cna_profile)
    clinical_rows = get_clinical_data(study)

    study_cache = _cache_json(cache, "study", study_metadata, force)
    profiles_cache = _cache_json(cache, "profiles", profiles, force)
    mutations_cache = _cache_json(cache, "mutations", mutation_rows, force)
    cna_cache = _cache_json(cache, "cna", cna_rows, force)
    clinical_cache = _cache_json(cache, "clinical", clinical_rows, force)

    mutations = parse_mutations(mutation_rows)
    clinical = parse_clinical(clinical_rows)
    cna = parse_discrete_cna(cna_rows)

    mutation_table = MutationTable(records=mutations)
    clinical_table = ClinicalTable(records=clinical)

    _write_tsv(
        output / "mutations.tsv",
        list(MutationRecord.model_fields),
        [record.model_dump(mode="python") for record in mutation_table.records],
    )
    _write_tsv(
        output / "clinical.tsv",
        ["sample", "survival_time", "survival_status", "therapy"],
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
    _write_tsv(
        output / "copy_number_calls.tsv",
        ["sample", "gene", "entrez_gene_id", "cna_value", "molecular_profile_id"],
        cna,
    )

    manifest: dict[str, Any] = {
        "study": study,
        "source": "cBioPortal",
        "public_api": f"{CBIOPORTAL_API}",
        "study_metadata": str(study_cache),
        "profiles": str(profiles_cache),
        "mutations": str(mutations_cache),
        "cna": str(cna_cache),
        "clinical": str(clinical_cache),
        "files": [],
    }
    for path in (study_cache, profiles_cache, mutations_cache, cna_cache, clinical_cache):
        manifest["files"].append(
            {
                "path": str(path),
                "size": path.stat().st_size,
                "sha256": _sha256(path),
            }
        )
    manifest_path = output / "cbioportal_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")

    return CbioPortalIngestResult(
        study=study,
        cache_dir=cache,
        output_dir=output,
        manifest=manifest_path,
        mutation_count=len(mutation_table.records),
        clinical_count=len(clinical_table.records),
        cna_count=len(cna),
    )
