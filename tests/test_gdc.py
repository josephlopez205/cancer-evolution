import csv
import hashlib

import pytest

from clonetrace.ingest.gdc import (
    GdcFile,
    download_open_file,
    ingest_gdc,
    parse_maf,
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
                        "size": 10,
                        "data_category": "Simple Nucleotide Variation",
                        "data_type": "Masked Somatic Mutation",
                    }
                ],
                "pagination": {"total": 1},
            }
        }

    monkeypatch.setattr("clonetrace.ingest.gdc._get_json", fake_get)

    files = search_open_files("TCGA-BRCA", "Simple Nucleotide Variation")

    assert len(files) == 1
    assert files[0].file_id == "maf-1"
    assert len(calls) == 1
    assert "TCGA-BRCA" in calls[0]
    assert "%22open%22" in calls[0]


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


def test_download_open_file_uses_cache_and_md5(tmp_path, monkeypatch) -> None:
    cache = tmp_path / "cache"
    cached = cache / "file-1" / "sample.maf"
    cached.parent.mkdir(parents=True)
    cached.write_text("sample", encoding="utf-8")
    md5 = hashlib.md5(cached.read_bytes()).hexdigest()
    file = GdcFile(
        file_id="file-1",
        file_name="sample.maf",
        md5sum=md5,
        size=cached.stat().st_size,
        data_category="Simple Nucleotide Variation",
        data_type="Masked Somatic Mutation",
    )

    def unexpected_download(*args, **kwargs) -> None:
        raise AssertionError("cached file should not be downloaded")

    monkeypatch.setattr("urllib.request.urlopen", unexpected_download)

    result = download_open_file(file, cache)

    assert result == cached
    assert result.read_text(encoding="utf-8") == "sample"


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

    def fake_search(project: str, category: str, page_size: int = 100) -> list[GdcFile]:
        if category == "Simple Nucleotide Variation":
            return [maf_file]
        if category == "Copy Number Variation":
            return [seg_file]
        return []

    def fake_download(file: GdcFile, cache_dir, force: bool = False):
        return maf if file.file_id == "maf-1" else seg

    monkeypatch.setattr("clonetrace.ingest.gdc.search_open_files", fake_search)
    monkeypatch.setattr("clonetrace.ingest.gdc.download_open_file", fake_download)

    result = ingest_gdc("TCGA-BRCA", cache_dir=cache, output_dir=output)

    assert result.mutation_count == 1
    assert result.segment_record_count == 1
    assert result.manifest.exists()
    assert (output / "mutations.tsv").exists()
    assert (output / "copy_number_segments.tsv").exists()
