import csv
import hashlib
import json
import urllib.parse

import pytest

from clonetrace.ingest.gdc import (
    GdcFile,
    fetch_open_file,
    ingest_gdc,
    parse_clinical_cases,
    parse_maf,
    parse_maf_content,
    parse_segment_content,
    parse_segments,
    search_open_files,
)
from clonetrace.schemas import CopyNumberSegment, MutationRecord


def test_search_open_files_only_accepts_api_hits(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []

    def fake_get(url: str) -> dict[str, object]:
        calls.append(url)
        return {
            "data": {
                "hits": [
                    {
                        "file_id": "maf-1",
                        "file_name": "sample.maf.gz",
                        "md5sum": "abc",
                        "file_size": 10,
                        "data_category": "Simple Nucleotide Variation",
                        "data_type": "Masked Somatic Mutation",
                        "data_format": "MAF",
                        "access": "open",
                    }
                ],
                "pagination": {"total": 1},
            }
        }

    monkeypatch.setattr("clonetrace.ingest.gdc._get_json", fake_get)

    files = search_open_files("TCGA-BRCA", "Masked Somatic Mutation", "MAF")

    assert len(files) == 1
    assert files[0].file_id == "maf-1"
    assert files[0].size == 10
    assert len(calls) == 1
    query = urllib.parse.parse_qs(urllib.parse.urlparse(calls[0]).query)
    filters = json.loads(query["filters"][0])
    filter_contents = filters["content"]
    assert {
        "op": "in",
        "content": {"field": "cases.project.project_id", "value": ["TCGA-BRCA"]},
    } in filter_contents
    assert {
        "op": "in",
        "content": {"field": "access", "value": ["open"]},
    } in filter_contents
    assert {
        "op": "in",
        "content": {"field": "data_type", "value": ["Masked Somatic Mutation"]},
    } in filter_contents
    assert {
        "op": "in",
        "content": {"field": "data_format", "value": ["MAF"]},
    } in filter_contents


def test_parse_maf(tmp_path) -> None:
    path = tmp_path / "sample.maf"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Tumor_Sample_Barcode",
                "Hugo_Symbol",
                "Chromosome",
                "Start_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
                "t_ref_count",
                "t_alt_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "Tumor_Sample_Barcode": "TCGA-01",
                "Hugo_Symbol": "TP53",
                "Chromosome": "chr17",
                "Start_Position": "7579472",
                "Reference_Allele": "C",
                "Tumor_Seq_Allele2": "T",
                "t_ref_count": "40",
                "t_alt_count": "10",
            }
        )

    mutations = parse_maf(path)

    assert len(mutations) == 1
    assert isinstance(mutations[0], MutationRecord)
    assert mutations[0].chrom == "17"
    assert mutations[0].ref_count == 40


def test_parse_segments(tmp_path) -> None:
    path = tmp_path / "sample.seg.txt"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "GDC_Aliquot_ID",
                "Chromosome",
                "Start",
                "End",
                "Copy_Number",
                "Minor_Copy_Number",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "GDC_Aliquot_ID": "aliquot-1",
                "Chromosome": "chr1",
                "Start": "100",
                "End": "500",
                "Copy_Number": "3",
                "Minor_Copy_Number": "1",
            }
        )

    segments = parse_segments(path)

    assert len(segments) == 1
    assert isinstance(segments[0], CopyNumberSegment)
    assert segments[0].total_cn == 3
    assert segments[0].minor_cn == 1


def test_fetch_open_file_verifies_checksums_without_persisting(monkeypatch) -> None:
    content = b"sample"
    md5 = hashlib.md5(content).hexdigest()
    file = GdcFile(
        file_id="file-1",
        file_name="sample.maf",
        md5sum=md5,
        size=len(content),
        data_category="Simple Nucleotide Variation",
        data_type="Masked Somatic Mutation",
    )

    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return None

        def read(self):
            return content

    monkeypatch.setattr("urllib.request.urlopen", lambda *args, **kwargs: Response())

    result = fetch_open_file(file)

    assert result == content


def test_fetch_open_file_rejects_non_open_metadata(monkeypatch) -> None:
    file = GdcFile(
        file_id="file-1",
        file_name="sample.maf",
        md5sum="",
        size=0,
        data_category="Simple Nucleotide Variation",
        data_type="Masked Somatic Mutation",
        access="controlled",
    )
    with pytest.raises(ValueError, match="non-open"):
        fetch_open_file(file)


def test_parse_api_content_in_memory() -> None:
    maf = (
        b"#version gdc-1.0.0\n"
        b"#tumor.aliquot fake-uuid\n"
        b"Hugo_Symbol\tTumor_Sample_Barcode\tChromosome\tStart_Position\t"
        b"Reference_Allele\tTumor_Seq_Allele2\tt_ref_count\tt_alt_count\n"
        b"TP53\tTCGA-01-0001-01A\t17\t100\tC\tT\t30\t10\n"
    )
    mutations = parse_maf_content(maf, "sample.maf")
    assert mutations[0].sample == "TCGA-01-0001-01A"

    segments = (
        b"GDC_Aliquot\tChromosome\tStart\tEnd\tCopy_Number\tMajor_Copy_Number\tMinor_Copy_Number\n"
        b"aliquot-1\tchr17\t100\t500\t3\t2\t1\n"
    )
    parsed_segments = parse_segment_content(
        segments,
        "sample.seg.txt",
        {"aliquot-1": "TCGA-01-0001-01A"},
    )
    assert parsed_segments[0].sample == "TCGA-01-0001-01A"
    assert parsed_segments[0].minor_cn == 1


def test_parse_clinical_cases_maps_survival_and_therapy() -> None:
    records = parse_clinical_cases(
        [
            {
                "submitter_id": "TCGA-01-0001",
                "samples": [{"submitter_id": "TCGA-01-0001-01A"}],
                "demographic": {"vital_status": "Dead", "days_to_death": 100},
                "diagnoses": [
                    {
                        "days_to_last_follow_up": None,
                        "treatments": [{"treatment_type": "Chemotherapy"}],
                    }
                ],
            }
        ]
    )
    assert records[0].sample == "TCGA-01-0001-01A"
    assert records[0].survival_time == 100
    assert records[0].survival_status == "event"
    assert records[0].therapy == ["Chemotherapy"]


def test_ingest_gdc_writes_manifest_and_tables(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cache = tmp_path / "cache"
    output = tmp_path / "processed"
    maf = tmp_path / "test.maf"
    seg = tmp_path / "test.seg.txt"

    maf.write_text(
        "Tumor_Sample_Barcode\tHugo_Symbol\tChromosome\tStart_Position\t"
        "Reference_Allele\tTumor_Seq_Allele2\tt_ref_count\tt_alt_count\n"
        "TCGA-01\tTP53\t17\t100\tC\tT\t30\t10\n",
        encoding="utf-8",
    )
    seg.write_text(
        "GDC_Aliquot_ID\tChromosome\tStart\tEnd\tCopy_Number\tMinor_Copy_Number\n"
        "aliquot-1\t1\t100\t500\t3\t1\n",
        encoding="utf-8",
    )

    maf_file = GdcFile(
        "maf-1",
        "test.maf",
        "",
        maf.stat().st_size,
        "Simple Nucleotide Variation",
        "Masked Somatic Mutation",
    )
    seg_file = GdcFile(
        "seg-1",
        "test.seg.txt",
        "",
        seg.stat().st_size,
        "Copy Number Variation",
        "Copy Number Segment",
    )

    def fake_search(
        project: str,
        data_type: str,
        data_format: str | None = None,
        include_sample_map: bool = False,
        page_size: int = 100,
    ) -> list[GdcFile]:
        if data_type == "Masked Somatic Mutation":
            return [maf_file]
        if data_type == "Allele-specific Copy Number Segment":
            return [seg_file]
        return []

    def fake_fetch(file: GdcFile) -> bytes:
        return maf.read_bytes() if file.file_id == "maf-1" else seg.read_bytes()

    monkeypatch.setattr("clonetrace.ingest.gdc.search_open_files", fake_search)
    monkeypatch.setattr("clonetrace.ingest.gdc.fetch_open_file", fake_fetch)
    monkeypatch.setattr(
        "clonetrace.ingest.gdc.fetch_open_cases",
        lambda project: [
            {
                "submitter_id": "TCGA-01",
                "samples": [{"submitter_id": "TCGA-01-0001-01A"}],
                "demographic": {"vital_status": "Alive"},
                "diagnoses": [],
            }
        ],
    )

    seg_file = GdcFile(
        "seg-1",
        "test.seg.txt",
        "",
        seg.stat().st_size,
        "Copy Number Variation",
        "Allele-specific Copy Number Segment",
        "TXT",
        "open",
        {"aliquot-1": "TCGA-01-0001-01A"},
    )

    result = ingest_gdc("TCGA-BRCA", output_dir=output)

    assert result.mutation_count == 1
    assert result.segment_record_count == 1
    assert result.clinical_count == 1
    assert result.manifest.exists()
    assert (output / "mutations.tsv").exists()
    assert (output / "copy_number_segments.tsv").exists()
    assert (output / "clinical.tsv").exists()
    assert (output / "cohort.tsv").exists()
    assert not cache.exists()
    assert {path.name for path in output.iterdir()} == {
        "mutations.tsv",
        "copy_number_segments.tsv",
        "clinical.tsv",
        "cohort.tsv",
        "manifest.json",
    }
    manifest = json.loads(result.manifest.read_text(encoding="utf-8"))
    assert manifest["raw_files_persisted"] is False
    assert all(item["access"] == "open" for item in manifest["files"])
