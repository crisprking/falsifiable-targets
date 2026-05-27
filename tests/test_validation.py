"""
Schema-validation tests.

Every shipped claim YAML must validate. If a contributor adds a claim
file that doesn't conform, this catches it in CI.

Also tests the validate_claim module itself, since it's a primary
researcher-facing entry point (the `ft-validate` CLI).
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from validate_claim import (
    VALID_CLAIM_TYPES,
    VALID_FIXTURE_SECTIONS,
    validate_claim_dict,
    validate_claim_file,
)

# -- Discover claim files automatically. Skip non-claim YAMLs like paralog_map. -

def _claim_files() -> list[Path]:
    files = []
    for p in sorted((ROOT / "claims").glob("*.yaml")):
        data = yaml.safe_load(p.read_text())
        if isinstance(data, dict) and "claim" in data:
            files.append(p)
    return files


CLAIM_FILES = _claim_files()


@pytest.mark.parametrize("path", CLAIM_FILES, ids=lambda p: p.name)
def test_shipped_claim_validates(path: Path):
    """Every shipped claim file must pass schema validation."""
    errors = validate_claim_file(path)
    assert not errors, (
        f"{path.name} failed schema validation:\n" +
        "\n".join(f"  - {e}" for e in errors)
    )


# -- Sentinel files use the same schema with extra fields ------------------

def test_all_sentinels_validate():
    """Every sentinel in v1_sentinels.yaml validates as a claim."""
    with open(ROOT / "sentinels" / "v1_sentinels.yaml") as f:
        suite = yaml.safe_load(f)

    bad = []
    for s in suite["sentinels"]:
        errors = validate_claim_dict(s)
        if errors:
            bad.append((s.get("id", "<unknown>"), errors))

    if bad:
        msg = ["Sentinels failed schema validation:"]
        for sid, errs in bad:
            msg.append(f"  {sid}:")
            for e in errs:
                msg.append(f"    - {e}")
        raise AssertionError("\n".join(msg))


# -- Validator catches actual problems -------------------------------------

class TestValidatorRejects:
    """The validator must reject malformed input. Each case is a known bug class."""

    def test_missing_claim_block(self):
        errors = validate_claim_dict({"fixture": {}})
        # Both pydantic and manual paths produce different exact wording;
        # check the error mentions 'claim'
        assert errors, "expected errors for missing 'claim' block"
        assert any("claim" in e.lower() for e in errors)

    def test_missing_claim_type(self):
        errors = validate_claim_dict({
            "claim": {
                "target_symbol": "TEST",
                "indication": "test",
                "mechanism": "test",
            }
        })
        assert any("claim_type" in e and "required" in e.lower() for e in errors)

    def test_invalid_claim_type(self):
        errors = validate_claim_dict({
            "claim": {
                "target_symbol": "TEST",
                "indication": "test",
                "mechanism": "test",
                "claim_type": "not_a_real_type",
            }
        })
        assert any("claim_type" in e and "not_a_real_type" in e for e in errors)

    def test_empty_target_symbol(self):
        errors = validate_claim_dict({
            "claim": {
                "target_symbol": "",
                "indication": "test",
                "mechanism": "test",
                "claim_type": "novel_target",
            }
        })
        assert errors  # at least one error

    def test_unknown_fixture_section(self):
        errors = validate_claim_dict({
            "claim": {
                "target_symbol": "TEST",
                "indication": "test",
                "mechanism": "test",
                "claim_type": "novel_target",
            },
            "fixture": {
                "made_up_section": {"foo": "bar"},
            },
        })
        assert any("made_up_section" in e for e in errors)

    def test_uniprot_id_too_short(self):
        errors = validate_claim_dict({
            "claim": {
                "target_symbol": "TEST",
                "indication": "test",
                "mechanism": "test",
                "claim_type": "novel_target",
                "uniprot_id": "X",
            }
        })
        assert any("uniprot_id" in e for e in errors)


class TestValidatorAccepts:
    """The validator must accept legitimate inputs."""

    def test_minimal_valid_claim(self):
        errors = validate_claim_dict({
            "claim": {
                "target_symbol": "TEST",
                "indication": "test indication",
                "mechanism": "test mechanism",
                "claim_type": "novel_target",
            }
        })
        assert errors == []

    def test_all_optional_fields(self):
        errors = validate_claim_dict({
            "claim": {
                "target_symbol": "TEST",
                "indication": "test",
                "mechanism": "test mechanism",
                "claim_type": "validated_mechanism",
                "uniprot_id": "P12345",
                "ensembl_id": "ENSG00000123456",
                "chembl_id": "CHEMBL1234",
                "source_url": "https://example.com",
                "citation_doi": "10.1234/test",
                "citation_pmid": "12345",
                "author": "test author",
                "date": "2026-05-26",
                "notes": "test notes",
            }
        })
        assert errors == []


# -- Constants sanity ------------------------------------------------------

def test_claim_types_match_smoke_test_enum():
    """validate_claim's VALID_CLAIM_TYPES must match smoke_test.ClaimType."""
    from smoke_test import ClaimType
    enum_values = {ct.value for ct in ClaimType}
    assert VALID_CLAIM_TYPES == enum_values, (
        f"drift between validate_claim.VALID_CLAIM_TYPES ({VALID_CLAIM_TYPES}) "
        f"and smoke_test.ClaimType ({enum_values})"
    )


def test_fixture_sections_match_adapter_protocol():
    """validate_claim's VALID_FIXTURE_SECTIONS must match adapters.SECTIONS."""
    from adapters import SECTIONS
    assert VALID_FIXTURE_SECTIONS == set(SECTIONS), (
        "drift between validate_claim.VALID_FIXTURE_SECTIONS and adapters.SECTIONS"
    )
