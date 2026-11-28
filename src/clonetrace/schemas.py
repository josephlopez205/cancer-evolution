"""Canonical data schemas for CloneTrace inputs and outputs."""

from typing import Literal

from pydantic import AliasChoices, BaseModel, ConfigDict, Field, field_validator, model_validator


class CanonicalModel(BaseModel):
    """Shared strict configuration for canonical CloneTrace records."""

    model_config = ConfigDict(
        extra="forbid",
        str_strip_whitespace=True,
        validate_assignment=True,
    )


class MutationRecord(CanonicalModel):
    """A somatic mutation observation in one sample."""

    sample: str = Field(validation_alias=AliasChoices("sample", "sample_id"))
    gene: str
    chrom: str
    pos: int = Field(ge=1)
    ref: str
    alt: str
    ref_count: int = Field(ge=0)
    alt_count: int = Field(ge=0)

    @field_validator("sample", "gene")
    @classmethod
    def validate_nonempty(cls, value: str) -> str:
        if not value:
            raise ValueError("must not be empty")
        return value

    @field_validator("chrom")
    @classmethod
    def normalize_chromosome(cls, value: str) -> str:
        chromosome = value.strip()
        if chromosome.lower().startswith("chr"):
            chromosome = chromosome[3:]
        chromosome = chromosome.upper() if chromosome.lower() in {"x", "y", "m"} else chromosome
        if not chromosome:
            raise ValueError("chromosome must not be empty")
        return chromosome

    @field_validator("ref", "alt")
    @classmethod
    def normalize_allele(cls, value: str) -> str:
        if not value:
            raise ValueError("allele must not be empty")
        return value.upper()


class CopyNumberSegment(CanonicalModel):
    """A copy-number segment in one sample."""

    sample: str = Field(validation_alias=AliasChoices("sample", "sample_id"))
    chrom: str
    start: int = Field(ge=0)
    end: int = Field(ge=1)
    total_cn: float = Field(ge=0, allow_inf_nan=False)
    minor_cn: float = Field(ge=0, allow_inf_nan=False)

    @field_validator("sample")
    @classmethod
    def validate_sample(cls, value: str) -> str:
        if not value:
            raise ValueError("sample must not be empty")
        return value

    @field_validator("chrom")
    @classmethod
    def normalize_chromosome(cls, value: str) -> str:
        chromosome = value.strip()
        if chromosome.lower().startswith("chr"):
            chromosome = chromosome[3:]
        chromosome = chromosome.upper() if chromosome.lower() in {"x", "y", "m"} else chromosome
        if not chromosome:
            raise ValueError("chromosome must not be empty")
        return chromosome

    @model_validator(mode="after")
    def validate_segment(self) -> "CopyNumberSegment":
        if self.end <= self.start:
            raise ValueError("end must be greater than start")
        if self.minor_cn > self.total_cn:
            raise ValueError("minor_cn must not exceed total_cn")
        return self


class PurityPloidyEstimate(CanonicalModel):
    """Externally estimated sample purity and/or tumor ploidy with provenance."""

    sample: str = Field(validation_alias=AliasChoices("sample", "sample_id"))
    purity: float | None = Field(default=None, ge=0, le=1, allow_inf_nan=False)
    ploidy: float | None = Field(default=None, gt=0, allow_inf_nan=False)
    source: str
    method: str | None = None

    @field_validator("sample", "source")
    @classmethod
    def validate_nonempty(cls, value: str) -> str:
        if not value:
            raise ValueError("must not be empty")
        return value

    @model_validator(mode="after")
    def require_estimate(self) -> "PurityPloidyEstimate":
        if self.purity is None and self.ploidy is None:
            raise ValueError("at least one of purity or ploidy must be provided")
        return self


class Survival(CanonicalModel):
    """Survival outcome for a clinical record."""

    time: float | None = Field(default=None, ge=0, allow_inf_nan=False)
    status: bool | Literal["event", "censored"] | int | None = None

    @field_validator("status")
    @classmethod
    def normalize_status(
        cls,
        value: bool | Literal["event", "censored"] | int | None,
    ) -> bool | Literal["event", "censored"] | int | None:
        if isinstance(value, int) and value not in {0, 1}:
            raise ValueError("survival status integer must be 0 or 1")
        return value


class ClinicalRecord(CanonicalModel):
    """Clinical and treatment annotations for one sample."""

    sample: str = Field(validation_alias=AliasChoices("sample", "sample_id"))
    survival: Survival | dict[str, object] | None = None
    survival_time: float | None = Field(default=None, ge=0, allow_inf_nan=False)
    survival_status: bool | Literal["event", "censored"] | int | None = None
    therapy: list[str] = Field(default_factory=list)

    @field_validator("sample")
    @classmethod
    def validate_sample(cls, value: str) -> str:
        if not value:
            raise ValueError("sample must not be empty")
        return value

    @field_validator("survival_status")
    @classmethod
    def normalize_status(
        cls,
        value: bool | Literal["event", "censored"] | int | None,
    ) -> bool | Literal["event", "censored"] | int | None:
        if isinstance(value, int) and value not in {0, 1}:
            raise ValueError("survival status integer must be 0 or 1")
        return value

    @field_validator("therapy")
    @classmethod
    def normalize_therapy(cls, value: list[str]) -> list[str]:
        therapies = [therapy.strip() for therapy in value if therapy.strip()]
        return list(dict.fromkeys(therapies))

    @model_validator(mode="after")
    def validate_survival(self) -> "ClinicalRecord":
        nested = self.survival
        if isinstance(nested, dict):
            nested = Survival.model_validate(nested)
            object.__setattr__(self, "survival", nested)
        if nested is not None:
            if self.survival_time is None:
                object.__setattr__(self, "survival_time", nested.time)
            if self.survival_status is None:
                object.__setattr__(self, "survival_status", nested.status)
        return self


class MutationTable(CanonicalModel):
    """Canonical collection of mutation records."""

    records: list[MutationRecord]


class CopyNumberTable(CanonicalModel):
    """Canonical collection of copy-number segments."""

    records: list[CopyNumberSegment]


class ClinicalTable(CanonicalModel):
    """Canonical collection of clinical records."""

    records: list[ClinicalRecord]


class CohortRecord(CanonicalModel):
    """One sample summary in a harmonized cohort table."""

    sample: str = Field(validation_alias=AliasChoices("sample", "sample_id"))
    source: str
    mutation_count: int = Field(default=0, ge=0)
    segment_count: int = Field(default=0, ge=0)
    survival_time: float | None = Field(default=None, ge=0, allow_inf_nan=False)
    survival_status: bool | Literal["event", "censored"] | int | None = None
    therapy: list[str] = Field(default_factory=list)

    @field_validator("sample", "source")
    @classmethod
    def validate_nonempty(cls, value: str) -> str:
        if not value:
            raise ValueError("must not be empty")
        return value

    @field_validator("therapy")
    @classmethod
    def normalize_therapy(cls, value: list[str]) -> list[str]:
        return list(dict.fromkeys(item.strip() for item in value if item.strip()))


class CohortTable(CanonicalModel):
    """Canonical collection of cohort sample summaries."""

    records: list[CohortRecord]


# Concise aliases for callers that prefer domain names over table-oriented names.
Mutation = MutationRecord
Segment = CopyNumberSegment
Clinical = ClinicalRecord
