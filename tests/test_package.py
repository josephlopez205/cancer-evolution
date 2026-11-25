from clonetrace import __version__
from clonetrace.schemas import Clinical, CopyNumberSegment, Mutation


def test_package_version() -> None:
    assert __version__ == "0.1.0"


def test_canonical_schema_module_imports() -> None:
    assert Mutation.__name__ == "MutationRecord"
    assert CopyNumberSegment.__name__ == "CopyNumberSegment"
    assert Clinical.__name__ == "ClinicalRecord"
