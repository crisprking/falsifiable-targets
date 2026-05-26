"""TYK2 external-audit pin: ensures the first external audit returns
SURVIVED under v1.1.0 with the documented fixture inputs. Hermetic -
runs fixture-only so it does not require network."""
from pathlib import Path
import sys

import pytest
import yaml

ROOT = Path("/kaggle/working/falsifiable-targets")
sys.path.insert(0, str(ROOT))

from smoke_test import TargetClaim, ClaimType, Verdict, run_audit, RuleStatus


def test_tyk2_audit_survives_under_v1_1_0_fixture_only():
    """
    TYK2 / psoriasis: a validated_mechanism claim with full
    deucravacitinib-quality evidence. Should SURVIVE.

    This test is the regression lock for the first external audit.
    If a future rule change (e.g. R6 scope expansion in v1.3.0+)
    flips this verdict, the change must be deliberate - which means
    updating this test with an explanation.
    """
    with open(ROOT / "claims" / "tyk2_psoriasis.yaml") as f:
        spec = yaml.safe_load(f)
    c = spec["claim"]
    claim = TargetClaim(
        target_symbol=c["target_symbol"],
        indication=c["indication"],
        mechanism=c["mechanism"],
        claim_type=ClaimType(c["claim_type"]),
        uniprot_id=c.get("uniprot_id"),
    )
    (verdict, score, _, substantive, _), results = run_audit(claim, spec["fixture"])
    assert verdict == Verdict.SURVIVED, (
        f"TYK2 expected SURVIVED under fixture-only v1.1.0, got {verdict.value}. "
        f"Substantive caveats: {[(c.rule_id, c.text[:50]) for c in substantive]}. "
        f"Per-rule: {[(r.rule_id, r.status.value) for r in results]}"
    )
    assert score == 0.0
    assert substantive == [], f"unexpected substantive caveats: {substantive}"


def test_tyk2_r6_not_applicable_under_v1_1_0():
    """
    Pin: R6 does NOT apply to validated_mechanism claims in v1.1.0.
    This is the documented scope limitation surfaced by the TYK2 audit.
    When/if v1.3.0 expands R6 to validated_mechanism, this test will
    fail loudly and force a deliberate update with rationale.
    """
    with open(ROOT / "claims" / "tyk2_psoriasis.yaml") as f:
        spec = yaml.safe_load(f)
    c = spec["claim"]
    claim = TargetClaim(
        target_symbol=c["target_symbol"],
        indication=c["indication"],
        mechanism=c["mechanism"],
        claim_type=ClaimType(c["claim_type"]),
        uniprot_id=c.get("uniprot_id"),
    )
    _, results = run_audit(claim, spec["fixture"])
    r6 = next(r for r in results if r.rule_id == "R6_chemistry_class_collapse")
    assert r6.status == RuleStatus.NOT_APPLICABLE, (
        f"R6 expected NOT_APPLICABLE for validated_mechanism claim under "
        f"v1.1.0, got {r6.status.value}. This may indicate an undocumented "
        f"R6 scope change."
    )
