from pathlib import Path

from clonetrace.ingest.cohort import build_cohort_table


def test_build_cohort_table_merges_canonical_tables(tmp_path: Path) -> None:
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()

    (input_dir / "mutations.tsv").write_text(
        "sample\tgene\tchrom\tpos\tref\talt\tref_count\talt_count\n"
        "sample-1\tTP53\t17\t100\tC\tT\t30\t10\n"
        "sample-1\tKRAS\t12\t200\tG\tA\t20\t5\n",
        encoding="utf-8",
    )
    (input_dir / "copy_number_segments.tsv").write_text(
        "sample\tchrom\tstart\tend\ttotal_cn\tminor_cn\n"
        "sample-1\t17\t100\t500\t3\t1\n"
        "sample-2\t1\t20\t100\t2\t1\n",
        encoding="utf-8",
    )
    (input_dir / "clinical.tsv").write_text(
        "sample\tsurvival_time\tsurvival_status\ttherapy\n"
        "sample-1\t36.5\tevent\ttamoxifen,radiotherapy\n",
        encoding="utf-8",
    )

    cohort_path = build_cohort_table(input_dir, output_dir, "cbioportal")
    lines = cohort_path.read_text(encoding="utf-8").splitlines()

    assert lines[0].split("\t") == [
        "sample",
        "source",
        "mutation_count",
        "segment_count",
        "survival_time",
        "survival_status",
        "therapy",
    ]
    assert "sample-1\tcbioportal\t2\t1\t36.5\tevent\ttamoxifen,radiotherapy" in lines
    assert "sample-2\tcbioportal\t0\t1\t\t\t" in lines
