import json

import pytest

from clonetrace.ingest.geo import ingest_geo, parse_clinical, validate_geo_accession


def test_validate_geo_accession() -> None:
    assert validate_geo_accession(" gse12345 ") == "GSE12345"

    with pytest.raises(ValueError):
        validate_geo_accession("SRP123")


def test_parse_clinical(tmp_path) -> None:
    path = tmp_path / "clinical.tsv"
    path.write_text(
        "sample\tOS_MONTHS\tOS_STATUS\ttherapy\n"
        "region-1\t36.5\t1:DECEASED\ttamoxifen,radiotherapy\n",
        encoding="utf-8",
    )

    records = parse_clinical(path)

    assert len(records) == 1
    assert records[0].sample == "region-1"
    assert records[0].survival_time == 36.5
    assert records[0].survival_status == "event"
    assert records[0].therapy == ["tamoxifen", "radiotherapy"]


def test_ingest_geo_requires_license_url(tmp_path) -> None:
    with pytest.raises(ValueError):
        ingest_geo(
            accession="GSE12345",
            processed_dir=tmp_path,
            license_url="not-a-url",
        )


def test_ingest_geo_loads_processed_files_and_excludes_raw_reads(tmp_path) -> None:
    source = tmp_path / "processed"
    source.mkdir()
    (source / "mutations.maf").write_text(
        "Tumor_Sample_Barcode\tHugo_Symbol\tChromosome\tStart_Position\t"
        "Reference_Allele\tTumor_Seq_Allele2\tt_ref_count\tt_alt_count\n"
        "region-1\tTP53\t17\t100\tC\tT\t30\t10\n",
        encoding="utf-8",
    )
    (source / "copy_number.seg").write_text(
        "GDC_Aliquot_ID\tChromosome\tStart\tEnd\tCopy_Number\tMinor_Copy_Number\n"
        "region-1\t17\t100\t500\t3\t1\n",
        encoding="utf-8",
    )
    (source / "clinical.tsv").write_text(
        "sample\tsurvival_time\tsurvival_status\ttherapy\n"
        "region-1\t36.5\tevent\ttamoxifen,radiotherapy\n",
        encoding="utf-8",
    )
    (source / "sample.bam").write_bytes(b"raw-read-placeholder")

    cache = tmp_path / "cache"
    output = tmp_path / "output"
    result = ingest_geo(
        accession="GSE12345",
        processed_dir=source,
        license_url="https://example.org/license",
        cache_dir=cache,
        output_dir=output,
    )

    assert result.mutation_count == 1
    assert result.segment_count == 1
    assert result.clinical_count == 1
    assert not (cache / "GSE12345" / "sample.bam").exists()
    assert (output / "mutations.tsv").exists()
    assert (output / "copy_number_segments.tsv").exists()
    assert (output / "clinical.tsv").exists()
    manifest = json.loads(result.manifest.read_text(encoding="utf-8"))
    assert manifest["raw_reads_excluded"] is True
    assert manifest["license_url"] == "https://example.org/license"
