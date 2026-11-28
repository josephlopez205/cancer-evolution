"""Loading externally estimated sample purity and ploidy with provenance."""

import csv
from pathlib import Path

from pydantic import ValidationError

from clonetrace.schemas import PurityPloidyEstimate

_REQUIRED_COLUMNS = {"sample", "purity", "ploidy", "source"}
_ALLOWED_COLUMNS = _REQUIRED_COLUMNS | {"method"}


def _parse_optional_float(value: str | None, *, column: str, line: int) -> float | None:
    if value is None or not value.strip():
        return None
    try:
        return float(value)
    except ValueError as error:
        raise ValueError(f"line {line}: {column} must be numeric or blank") from error


def read_purity_ploidy_estimates(path: str | Path) -> list[PurityPloidyEstimate]:
    """Read a strict TSV with sample, purity, ploidy, and source columns.

    Purity or ploidy may be blank when the other estimate is present. The
    source column records provenance, such as a tool name or source file.
    """

    estimates: list[PurityPloidyEstimate] = []
    seen_samples: set[str] = set()
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("purity/ploidy TSV must include a header")
        headers = reader.fieldnames
        if len(headers) != len(set(headers)):
            raise ValueError("purity/ploidy TSV contains duplicate column names")
        missing = _REQUIRED_COLUMNS - set(headers)
        unknown = set(headers) - _ALLOWED_COLUMNS
        if missing:
            raise ValueError(f"purity/ploidy TSV is missing columns: {', '.join(sorted(missing))}")
        if unknown:
            columns = ", ".join(sorted(unknown))
            raise ValueError(f"purity/ploidy TSV has unsupported columns: {columns}")

        for line, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"line {line}: row has more values than the header")
            if not row or all(not value for value in row.values() if value is not None):
                continue
            sample = (row.get("sample") or "").strip()
            if sample in seen_samples:
                raise ValueError(
                    f"line {line}: duplicate purity/ploidy estimate for sample {sample!r}"
                )

            estimate_data = {
                "sample": sample,
                "purity": _parse_optional_float(row.get("purity"), column="purity", line=line),
                "ploidy": _parse_optional_float(row.get("ploidy"), column="ploidy", line=line),
                "source": row.get("source") or "",
                "method": row.get("method") or None,
            }
            try:
                estimate = PurityPloidyEstimate.model_validate(estimate_data)
            except ValidationError as error:
                raise ValueError(f"line {line}: invalid purity/ploidy estimate: {error}") from error
            estimates.append(estimate)
            seen_samples.add(estimate.sample)

    return estimates
