"""Pytest sentinel suite - runs every sentinel as a parametrized test."""
from pathlib import Path
import sys
import yaml
import pytest

ROOT = Path("/kaggle/working/falsifiable-targets")
sys.path.insert(0, str(ROOT))

from smoke_test import ClaimType, TargetClaim, Verdict, run_audit, RULES


# Locked v1.1.0 ruleset SHA. Computed from the canonical ASCII rule descriptions.
RULESET_SHA_V1_1_0 = "2f9aab7d0ebc209f62c16eb35be31bc5b65fa2eb09adc02bea5bff5176269b32"


def _load_sentinels():
    with open(ROOT / "sentinels" / "v1_sentinels.yaml") as f:
        return yaml.safe_load(f)["sentinels"]


@pytest.mark.parametrize("sentinel", _load_sentinels(), ids=lambda s: s["id"])
def test_sentinel(sentinel):
    c = sentinel["claim"]
    claim = TargetClaim(
        target_symbol=c["target_symbol"],
        indication=c["indication"],
        mechanism=c["mechanism"],
        claim_type=ClaimType(c["claim_type"]),
        uniprot_id=c.get("uniprot_id"),
    )
    (verdict, _, _, _, _), _ = run_audit(claim, sentinel["fixture"])
    expected = Verdict(sentinel["expected_verdict"])
    sid = sentinel["id"]
    assert verdict == expected, (
        f"{sid}: expected {expected.value}, got {verdict.value}"
    )


def test_ruleset_sha_stable():
    """Locks the v1.1.0 ruleset hash. Any rule change must bump the version."""
    import hashlib
    import json
    descriptor = {
        "version": "1.1.0",
        "rules": [
            {"rule_id": r.rule_id, "version": r.version, "description": r.description}
            for r in sorted(RULES, key=lambda x: x.rule_id)
        ],
    }
    canon = json.dumps(descriptor, sort_keys=True, separators=(",", ":"))
    sha = hashlib.sha256(canon.encode()).hexdigest()
    assert sha == RULESET_SHA_V1_1_0, (
        f"Ruleset SHA drifted: {sha}. "
        "If this change is intentional, bump ruleset version."
    )


def test_ipi1_inaugural_audit_two_substantive_caveats():
    """
    Pins the Ipi1 inaugural audit result under v1.1.0:
    FALSIFIED_WITH_CAVEATS with exactly 2 substantive caveats
    (R1 orthology disagreement + R7 selectivity gap).

    This test catches the specific R7 regression that the
    sentinel-verdict test cannot catch: under v1.0.0 R7 logic,
    Ipi1 would have 1 caveat not 2, but the verdict would still
    be FALSIFIED_WITH_CAVEATS.
    """
    import yaml
    with open(ROOT / "claims" / "ipi1_madurella.yaml") as f:
        spec = yaml.safe_load(f)
    c = spec["claim"]
    claim = TargetClaim(
        target_symbol=c["target_symbol"],
        indication=c["indication"],
        mechanism=c["mechanism"],
        claim_type=ClaimType(c["claim_type"]),
        uniprot_id=c.get("uniprot_id"),
    )
    (verdict, score, _, substantive, _), _ = run_audit(claim, spec["fixture"])
    assert verdict == Verdict.FALSIFIED_WITH_CAVEATS, f"got {verdict.value}"
    assert score == 0.0, f"got score {score}"
    rule_ids = sorted(c.rule_id for c in substantive)
    assert rule_ids == ["R1_orthology", "R7_selectivity_counterscreen"], (
        f"expected substantive caveats on R1 and R7, got {rule_ids}"
    )


def test_novel_target_gap_has_r7_substantive_caveat():
    """
    Pins R7 v1.1.0 behavior on the dedicated sentinel:
    a novel target with majority orthology and no selectivity data
    must produce exactly one substantive caveat, from R7.
    """
    sentinels = _load_sentinels()
    gap = next(s for s in sentinels if s["id"] == "NOVEL_TARGET_SELECTIVITY_GAP")
    c = gap["claim"]
    claim = TargetClaim(
        target_symbol=c["target_symbol"],
        indication=c["indication"],
        mechanism=c["mechanism"],
        claim_type=ClaimType(c["claim_type"]),
        uniprot_id=c.get("uniprot_id"),
    )
    (verdict, _, _, substantive, _), _ = run_audit(claim, gap["fixture"])
    assert verdict == Verdict.FALSIFIED_WITH_CAVEATS
    rule_ids = [c.rule_id for c in substantive]
    assert rule_ids == ["R7_selectivity_counterscreen"], (
        f"expected exactly one substantive caveat from R7, got {rule_ids}"
    )
