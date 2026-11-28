import sys

import pytest

from clonetrace.cluster import (
    PyCloneObservation,
    PyCloneVIError,
    PyCloneVINotFoundError,
    run_pyclonevi,
)
from clonetrace.schemas import MutationRecord


def _observations() -> list[PyCloneObservation]:
    return [
        PyCloneObservation(
            mutation=MutationRecord(
                sample="sample-1",
                gene=f"GENE{index}",
                chrom="1",
                pos=index + 1,
                ref="A",
                alt="C",
                ref_count=70,
                alt_count=30,
            ),
            major_cn=2,
            minor_cn=1,
            normal_cn=2,
            tumour_content=0.75,
        )
        for index in range(2)
    ]


def _write_fake_pyclonevi(path, *, omit_last_result: bool = False) -> None:
    limit = "rows[:-1]" if omit_last_result else "rows"
    path.write_text(
        "\n".join(
            [
                "import csv",
                "import pathlib",
                "import sys",
                "args = sys.argv[1:]",
                "def option(name):",
                "    return args[args.index(name) + 1]",
                "if args[0] == 'fit':",
                "    pathlib.Path(option('--out-file')).write_text(option('--in-file'))",
                "else:",
                "    input_path = pathlib.Path(option('--in-file')).read_text()",
                "    with open(input_path, newline='') as source:",
                "        rows = list(csv.DictReader(source, delimiter='\\t'))",
                "    with open(option('--out-file'), 'w', newline='') as output:",
                "        fields = ['mutation_id', 'sample_id', 'cluster_id', "
                "'cellular_prevalence', 'cellular_prevalence_std', "
                "'cluster_assignment_prob']",
                "        writer = csv.DictWriter(output, fieldnames=fields, delimiter='\\t')",
                "        writer.writeheader()",
                f"        for row in reversed({limit}):",
                "            index = int(row['mutation_id'].rsplit('-', 1)[1])",
                "            writer.writerow({",
                "                'mutation_id': row['mutation_id'],",
                "                'sample_id': row['sample_id'],",
                "                'cluster_id': index,",
                "                'cellular_prevalence': 0.6 + index * 0.1,",
                "                'cellular_prevalence_std': 0.05,",
                "                'cluster_assignment_prob': 0.9,",
                "            })",
            ]
        ),
        encoding="utf-8",
    )


def test_pyclonevi_writes_input_runs_both_stages_and_normalizes_output(tmp_path) -> None:
    executable_script = tmp_path / "fake_pyclonevi.py"
    _write_fake_pyclonevi(executable_script)
    output_dir = tmp_path / "pyclonevi-run"

    result = run_pyclonevi(
        _observations(),
        output_dir,
        executable=[sys.executable, str(executable_script)],
        num_clusters=4,
        density="beta-binomial",
        seed=17,
    )

    assert result.sample == "sample-1"
    assert result.method == "pyclonevi"
    assert result.cluster_count == 2
    assert [item.cluster_id for item in result.assignments] == [0, 1]
    assert result.assignments[0].mutation.gene == "GENE0"
    assert result.assignments[1].cellular_prevalence == pytest.approx(0.7)
    assert result.assignments[0].cluster_probabilities is None
    assert [cluster.mutation_count for cluster in result.clusters] == [1, 1]
    assert result.input_path.is_file()
    assert result.fit_path.is_file()
    assert result.results_path.is_file()
    assert "tumour_content" in result.input_path.read_text(encoding="utf-8").splitlines()[0]


def test_pyclonevi_rejects_missing_result_rows(tmp_path) -> None:
    executable_script = tmp_path / "fake_pyclonevi.py"
    _write_fake_pyclonevi(executable_script, omit_last_result=True)

    with pytest.raises(PyCloneVIError, match="omitted results"):
        run_pyclonevi(
            _observations(),
            tmp_path / "run",
            executable=[sys.executable, str(executable_script)],
        )


def test_pyclonevi_surfaces_nonzero_process_exit(tmp_path) -> None:
    executable_script = tmp_path / "failing_pyclonevi.py"
    executable_script.write_text(
        "import sys\nprint('model failed', file=sys.stderr)\nsys.exit(3)\n",
        encoding="utf-8",
    )

    with pytest.raises(PyCloneVIError, match="model failed"):
        run_pyclonevi(
            _observations(),
            tmp_path / "run",
            executable=[sys.executable, str(executable_script)],
        )


def test_pyclonevi_reports_unavailable_executable(tmp_path) -> None:
    with pytest.raises(PyCloneVINotFoundError, match="Install PyClone-VI"):
        run_pyclonevi(
            _observations(),
            tmp_path / "run",
            executable="pyclone-vi-definitely-not-installed",
        )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("major_cn", -1, "major_cn >= minor_cn"),
        ("normal_cn", 0, "normal_cn must be greater"),
        ("tumour_content", 1.1, "tumour_content must be between"),
    ],
)
def test_pyclone_observation_rejects_invalid_inputs(field, value, message: str) -> None:
    arguments = {
        "mutation": _observations()[0].mutation,
        "major_cn": 2,
        "minor_cn": 1,
        "normal_cn": 2,
        "tumour_content": 0.75,
    }
    arguments[field] = value

    with pytest.raises(ValueError, match=message):
        PyCloneObservation(**arguments)
