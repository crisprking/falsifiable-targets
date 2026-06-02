"""TYK2 external-audit pins for v1.2.0.

Under v1.1.0 R6 did NOT apply to validated_mechanism claims, so the
v1.1.0 TYK2 audit returned SURVIVED with R6=NOT_APPLICABLE. Under
v1.2.0 R6 DOES apply to validated_mechanism (scope expansion), but
the new heuristic substantive-caveat path requires live paralog
chemistry data. In fixture-only mode without paralog counts, R6
falls through to PASSED without caveat - so TYK2 fixture-only still
returns SURVIVED. The live-mode TYK2 audit (with paralog map loaded)
is expected to fire the R6 substantive caveat; that's exercised in
the heuristic-path test below using a synthetic fixture.
"""
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from smoke_test import ClaimType, TargetClaim, Verdict, run_audit


def test_tyk2_audit_survives_under_v1_2_0_fixture_only():
    """TYK2 fixture-only mode under v1.2.0 returns SURVIVED.
    The R6 heuristic requires chembl_paralog_compound_counts which
    only the live ChEMBL adapter populates."""
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
        f"TYK2 expected SURVIVED under fixture-only v1.2.0, got {verdict.value}. "
        f"Substantive caveats: {[(c.rule_id, c.text[:50]) for c in substantive]}. "
        f"Per-rule: {[(r.rule_id, r.status.value) for r in results]}"
    )
    assert score == 0.0
    assert substantive == [], f"unexpected substantive caveats: {substantive}"


def test_tyk2_r6_applies_under_v1_2_0():
    """Pin: R6 DOES apply to validated_mechanism claims in v1.2.0."""
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
    from smoke_test import RULES
    r6 = next(r for r in RULES if r.rule_id == "R6_chemistry_class_collapse")
    assert r6.applies_to(claim), (
        "R6 should apply to validated_mechanism claims under v1.2.0; "
        "applies_to returned False."
    )


def test_r6_paralog_ratio_substantive_caveat_fires():
    """A validated_mechanism claim with paralog pool >= 2x primary must
    produce exactly one substantive caveat from R6."""
    claim = TargetClaim(
        target_symbol="TestKinase",
        indication="test indication",
        mechanism="test mechanism",
        claim_type=ClaimType.VALIDATED_MECHANISM,
        uniprot_id="P00000",
    )
    fixture = {
        "orthology": {"sources_agreeing": 4, "sources_total": 4},
        "genetics": {"gwas_hits": 5, "mendelian_evidence": True},
        "chemistry": {
            "chembl_distinct_compounds": 100,
            "chembl_paralog_compound_counts": {"paralog_X": 250},
        },
        "expression": {"target_tissue_expressed": True},
        "reproducibility": {"retracted": False, "rebuttals_count": 0, "independent_replications": 10},
        "selectivity": {"selectivity_data": True, "off_targets_in_indication_relevant_tissue": False,
                        "selectivity_index_log": 2.0},
    }
    (verdict, _, _, substantive, _), results = run_audit(claim, fixture)
    assert verdict == Verdict.FALSIFIED_WITH_CAVEATS, f"got {verdict.value}"
    r6_caveats = [c for c in substantive if c.rule_id == "R6_chemistry_class_collapse"]
    assert len(r6_caveats) == 1, (
        f"expected exactly 1 R6 substantive caveat, got {len(r6_caveats)}. "
        f"All substantive: {[(c.rule_id, c.text[:60]) for c in substantive]}"
    )
    assert "paralog chemistry pool" in r6_caveats[0].text.lower()


def test_r6_paralog_ratio_below_threshold_does_not_fire():
    """Paralog pool only 1.5x primary -> below 2.0 threshold -> no caveat."""
    claim = TargetClaim(
        target_symbol="TestKinase",
        indication="test indication",
        mechanism="test mechanism",
        claim_type=ClaimType.VALIDATED_MECHANISM,
        uniprot_id="P00000",
    )
    fixture = {
        "orthology": {"sources_agreeing": 4, "sources_total": 4},
        "genetics": {"gwas_hits": 5, "mendelian_evidence": True},
        "chemistry": {
            "chembl_distinct_compounds": 100,
            "chembl_paralog_compound_counts": {"paralog_X": 150},
        },
        "expression": {"target_tissue_expressed": True},
        "reproducibility": {"retracted": False, "rebuttals_count": 0, "independent_replications": 10},
        "selectivity": {"selectivity_data": True, "off_targets_in_indication_relevant_tissue": False,
                        "selectivity_index_log": 2.0},
    }
    (verdict, _, _, substantive, _), _ = run_audit(claim, fixture)
    assert verdict == Verdict.SURVIVED, f"expected SURVIVED, got {verdict.value}"
    assert not any(c.rule_id == "R6_chemistry_class_collapse" for c in substantive)
