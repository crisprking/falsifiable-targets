"""
Claim YAML schema validation.

Researchers writing their own claim files need fast feedback on whether
their YAML is well-formed BEFORE they run a 30-second live audit. This
module provides both a programmatic API and a CLI:

    ft-validate claims/my_claim.yaml
    python validate_claim.py claims/my_claim.yaml

Validation strategy:
  - Tries pydantic v2 first (if installed via `pip install ...[validate]`)
  - Falls back to manual schema check if pydantic missing
  - Both paths produce the same human-readable error report

The schema below is the canonical definition of a claim YAML. Any change
here must be reflected in docs/CLAIM_SCHEMA.md.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import yaml

# -- Schema spec (used by both pydantic and manual fallback) ---------------

VALID_CLAIM_TYPES = {
    "validated_mechanism",
    "novel_target",
    "chemistry_series",
    "extraordinary_claim",
}

VALID_VERDICTS = {
    "SURVIVED",
    "FALSIFIED_WITH_CAVEATS",
    "FALSIFIED",
    "INSUFFICIENT_DATA",
}

REQUIRED_CLAIM_FIELDS = {"target_symbol", "indication", "mechanism", "claim_type"}
OPTIONAL_CLAIM_FIELDS = {
    "uniprot_id",
    "ensembl_id",
    "chembl_id",
    "notes",
    "source_url",       # provenance URL (FDA announcement, paper, etc.)
    "citation_doi",     # DOI of the primary claim source
    "citation_pmid",    # PubMed ID
    "author",           # who authored this claim file (not the science)
    "date",             # YYYY-MM-DD when the claim was filed
}

# Sections that the fixture dict may supply. Adapters serve the same names.
VALID_FIXTURE_SECTIONS = {
    "orthology", "chemistry", "genetics", "expression",
    "reproducibility", "selectivity", "structure",
}


class ValidationError(Exception):
    """Raised when a claim YAML fails schema validation."""

    def __init__(self, errors: list[str]):
        self.errors = errors
        msg = "Claim validation failed:\n" + "\n".join(f"  - {e}" for e in errors)
        super().__init__(msg)


# -- Pydantic path (preferred when available) ------------------------------

def _validate_with_pydantic(data: dict[str, Any]) -> list[str]:
    """Returns a list of error messages, empty if valid."""
    try:
        from pydantic import BaseModel, Field, field_validator
        from pydantic import ValidationError as PydanticVE
    except ImportError:
        return _validate_manually(data)

    class ClaimModel(BaseModel):
        target_symbol: str = Field(min_length=1, max_length=64)
        indication: str = Field(min_length=1, max_length=256)
        mechanism: str = Field(min_length=1, max_length=2048)
        claim_type: str
        uniprot_id: str | None = Field(default=None, max_length=16)
        ensembl_id: str | None = Field(default=None, max_length=32)
        chembl_id: str | None = Field(default=None, max_length=32)
        notes: str | None = None
        source_url: str | None = None
        citation_doi: str | None = None
        citation_pmid: str | None = None
        author: str | None = None
        date: str | None = None

        @field_validator("claim_type")
        @classmethod
        def claim_type_must_be_valid(cls, v):
            if v not in VALID_CLAIM_TYPES:
                raise ValueError(
                    f"claim_type {v!r} not in {sorted(VALID_CLAIM_TYPES)}"
                )
            return v

        @field_validator("uniprot_id")
        @classmethod
        def uniprot_id_format(cls, v):
            if v is None:
                return v
            # UniProt accessions: 6 or 10 chars, [OPQ][0-9][A-Z0-9]{3}[0-9] etc.
            # Loose check - just alnum and reasonable length
            if not v.replace("-", "").isalnum():
                raise ValueError(f"uniprot_id {v!r} has non-alphanumeric chars")
            if len(v) < 6 or len(v) > 12:
                raise ValueError(f"uniprot_id {v!r} length {len(v)} not in [6, 12]")
            return v

    class SpecModel(BaseModel):
        claim: ClaimModel
        fixture: dict[str, dict] = Field(default_factory=dict)
        # Optional sentinel-only fields - tolerated but not required
        expected_verdict: str | None = None
        id: str | None = None

        @field_validator("fixture")
        @classmethod
        def fixture_sections_known(cls, v):
            unknown = set(v.keys()) - VALID_FIXTURE_SECTIONS
            if unknown:
                raise ValueError(
                    f"fixture has unknown sections {sorted(unknown)}. "
                    f"Valid: {sorted(VALID_FIXTURE_SECTIONS)}"
                )
            return v

        @field_validator("expected_verdict")
        @classmethod
        def verdict_valid(cls, v):
            if v is not None and v not in VALID_VERDICTS:
                raise ValueError(f"expected_verdict {v!r} not in {sorted(VALID_VERDICTS)}")
            return v

    try:
        SpecModel.model_validate(data)
        return []
    except PydanticVE as e:
        return [
            f"{'.'.join(str(p) for p in err['loc'])}: {err['msg']}"
            for err in e.errors()
        ]


# -- Manual path (no pydantic) ---------------------------------------------

def _validate_manually(data: dict[str, Any]) -> list[str]:
    errors: list[str] = []

    if not isinstance(data, dict):
        return [f"top-level must be a dict, got {type(data).__name__}"]

    if "claim" not in data:
        return ["missing required top-level field: 'claim'"]

    claim = data["claim"]
    if not isinstance(claim, dict):
        errors.append(f"claim must be a dict, got {type(claim).__name__}")
        return errors

    # Required fields
    for f in REQUIRED_CLAIM_FIELDS:
        if f not in claim:
            errors.append(f"claim.{f}: required field missing")
        elif not isinstance(claim[f], str) or not claim[f].strip():
            errors.append(f"claim.{f}: must be a non-empty string")

    # Claim type membership
    if "claim_type" in claim and claim["claim_type"] not in VALID_CLAIM_TYPES:
        errors.append(
            f"claim.claim_type: {claim['claim_type']!r} not in "
            f"{sorted(VALID_CLAIM_TYPES)}"
        )

    # Optional uniprot_id format
    upi = claim.get("uniprot_id")
    if upi is not None:
        if not isinstance(upi, str):
            errors.append(f"claim.uniprot_id: must be a string, got {type(upi).__name__}")
        elif not (6 <= len(upi) <= 12):
            errors.append(f"claim.uniprot_id: length {len(upi)} not in [6, 12]")

    # Unknown claim fields - warn (don't error - allow user notes)
    known = REQUIRED_CLAIM_FIELDS | OPTIONAL_CLAIM_FIELDS
    extra = set(claim.keys()) - known
    if extra:
        errors.append(
            f"claim has unknown field(s) {sorted(extra)}. "
            f"Known: {sorted(known)}. (If intentional metadata, prefix with 'x_'.)"
        )

    # Fixture sections
    fixture = data.get("fixture", {})
    if not isinstance(fixture, dict):
        errors.append(f"fixture must be a dict, got {type(fixture).__name__}")
    else:
        unknown = set(fixture.keys()) - VALID_FIXTURE_SECTIONS
        if unknown:
            errors.append(
                f"fixture has unknown section(s) {sorted(unknown)}. "
                f"Valid: {sorted(VALID_FIXTURE_SECTIONS)}"
            )
        for section_name, section_val in fixture.items():
            if not isinstance(section_val, dict):
                errors.append(
                    f"fixture.{section_name}: must be a dict, "
                    f"got {type(section_val).__name__}"
                )

    # Optional sentinel fields
    if "expected_verdict" in data and data["expected_verdict"] not in VALID_VERDICTS:
        errors.append(
            f"expected_verdict: {data['expected_verdict']!r} not in "
            f"{sorted(VALID_VERDICTS)}"
        )

    return errors


# -- Public API ------------------------------------------------------------

def validate_claim_dict(data: dict[str, Any]) -> list[str]:
    """Validate a parsed claim dict. Returns list of error strings (empty if valid)."""
    try:
        import pydantic  # noqa: F401
        return _validate_with_pydantic(data)
    except ImportError:
        return _validate_manually(data)


def validate_claim_file(path: str | Path) -> list[str]:
    """Validate a claim YAML file by path. Returns list of error strings."""
    p = Path(path)
    if not p.exists():
        return [f"file not found: {p}"]
    try:
        with open(p) as f:
            data = yaml.safe_load(f)
    except yaml.YAMLError as e:
        return [f"YAML parse error: {e}"]
    if data is None:
        return ["file is empty or contains only comments"]
    return validate_claim_dict(data)


# -- CLI -------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate a claim YAML against the falsifiable-targets schema.",
        epilog="Returns exit code 0 if valid, 1 if invalid.",
    )
    parser.add_argument("claim_path", help="Path to claim YAML file")
    parser.add_argument(
        "--quiet", action="store_true",
        help="No output on success (exit code only)",
    )
    args = parser.parse_args()

    errors = validate_claim_file(args.claim_path)

    if errors:
        print(f"FAIL: {args.claim_path}", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        return 1

    if not args.quiet:
        print(f"OK: {args.claim_path} validates against the claim schema.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
