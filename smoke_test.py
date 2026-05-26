"""
Smoke test - stdlib only. Runs every sentinel against the v1.0.0 rules
and asserts the expected verdict. If this passes 9/9, v1.0.0 is shippable.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

import yaml


class ClaimType(str, Enum):
    VALIDATED_MECHANISM = "validated_mechanism"
    NOVEL_TARGET = "novel_target"
    CHEMISTRY_SERIES = "chemistry_series"
    EXTRAORDINARY_CLAIM = "extraordinary_claim"


class FalsificationTier(str, Enum):
    PUBLIC_DATA_LOOKUP = "public_data_lookup"
    CHEAP_IN_SILICO = "cheap_in_silico"
    TARGETED_ASSAY = "targeted_assay"
    COHORT_STUDY = "cohort_study"
    CLINICAL_TRIAL = "clinical_trial"


class Verdict(str, Enum):
    SURVIVED = "SURVIVED"
    FALSIFIED_WITH_CAVEATS = "FALSIFIED_WITH_CAVEATS"
    FALSIFIED = "FALSIFIED"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


class CaveatKind(str, Enum):
    OPERATIONAL = "operational"
    SUBSTANTIVE = "substantive"


class RuleStatus(str, Enum):
    PASSED = "passed"
    FALSIFIED = "falsified"
    NOT_APPLICABLE = "not_applicable"
    ABSTAINED = "abstained"


@dataclass(frozen=True)
class Caveat:
    kind: CaveatKind
    rule_id: str
    text: str

    @property
    def is_substantive(self) -> bool:
        return self.kind == CaveatKind.SUBSTANTIVE


@dataclass
class TargetClaim:
    target_symbol: str
    indication: str
    mechanism: str
    claim_type: ClaimType
    uniprot_id: str | None = None


@dataclass
class RuleResult:
    rule_id: str
    status: RuleStatus
    confidence: float
    falsification_tier: FalsificationTier | None = None
    falsification_experiment: str | None = None
    evidence: dict = field(default_factory=dict)
    caveats: list = field(default_factory=list)


class FixtureAdapter:
    def __init__(self, fixture):
        self._f = fixture

    def get(self, section, claim):
        return self._f.get(section, {}) or {}


# ---------------------------------------------------------------------------
# The seven rules of v1.0.0
# ---------------------------------------------------------------------------

class R1_Orthology:
    rule_id = "R1_orthology"
    version = "1.0.0"
    description = (
        "For pathogen/novel-target claims, at least 2 of 4 ortholog DBs must "
        "agree that the protein exists in the claimed organism's lineage and "
        "is distinct from any human ortholog within selectivity-relevant divergence."
    )

    def applies_to(self, claim):
        return claim.claim_type in (ClaimType.NOVEL_TARGET, ClaimType.EXTRAORDINARY_CLAIM)

    def evaluate(self, claim, io):
        d = io.get("orthology", claim)
        st = d.get("sources_total")
        sa = d.get("sources_agreeing")
        if st is None or sa is None:
            return RuleResult(
                self.rule_id, RuleStatus.ABSTAINED, 0.0,
                caveats=[Caveat(CaveatKind.OPERATIONAL, self.rule_id,
                                "no orthology data; rule abstains")],
            )
        if sa == 0 and st >= 3:
            return RuleResult(
                self.rule_id, RuleStatus.FALSIFIED, 0.97,
                FalsificationTier.PUBLIC_DATA_LOOKUP,
                f"Query OrthoDB, OMA, EggNOG, InParanoid for {claim.target_symbol}. All {st} returned zero orthologs.",
                {"sources_agreeing": sa, "sources_total": st},
            )
        majority = (st // 2) + 1
        if sa < majority:
            return RuleResult(
                self.rule_id, RuleStatus.PASSED, 0.5,
                caveats=[Caveat(
                    CaveatKind.SUBSTANTIVE, self.rule_id,
                    f"orthology sources disagree ({sa}/{st}, strict majority requires {majority}); "
                    "selectivity claim weakened until resolved structurally"
                )],
                evidence={"sources_agreeing": sa, "sources_total": st},
            )
        return RuleResult(
            self.rule_id, RuleStatus.PASSED, 0.9,
            evidence={"sources_agreeing": sa, "sources_total": st},
        )


class R2_ChemistrySupport:
    rule_id = "R2_chemistry_support"
    version = "1.0.0"
    description = (
        "For chemistry-series and validated-mechanism claims, ChEMBL must contain "
        "a meaningful number of distinct compounds. Zero compounds for a "
        "chemistry-series claim is a contradiction; for a novel-target claim, "
        "absence of chemistry is expected and the rule does not apply."
    )

    def applies_to(self, claim):
        return claim.claim_type in (ClaimType.CHEMISTRY_SERIES, ClaimType.VALIDATED_MECHANISM)

    def evaluate(self, claim, io):
        d = io.get("chemistry", claim)
        n = d.get("chembl_distinct_compounds")
        if n is None:
            return RuleResult(
                self.rule_id, RuleStatus.ABSTAINED, 0.0,
                caveats=[Caveat(CaveatKind.OPERATIONAL, self.rule_id,
                                "no chemistry data; rule abstains")],
            )
        if n < 1:
            return RuleResult(
                self.rule_id, RuleStatus.FALSIFIED, 0.85,
                FalsificationTier.PUBLIC_DATA_LOOKUP,
                f"Query ChEMBL by target {claim.target_symbol}. Found {n} distinct compounds; chemistry-dependent claim cannot stand.",
                {"chembl_distinct_compounds": n},
            )
        return RuleResult(
            self.rule_id, RuleStatus.PASSED, 0.9,
            evidence={"chembl_distinct_compounds": n},
        )


class R3_GeneticsSupport:
    rule_id = "R3_genetics_support"
    version = "1.0.0"
    description = (
        "For human-disease validated-mechanism claims, the target should have "
        "at least one of: (a) GWAS hits in the indication, (b) Mendelian "
        "evidence with a phenotype consistent with the mechanism, (c) somatic "
        "driver evidence for cancer indications."
    )

    def applies_to(self, claim):
        return claim.claim_type == ClaimType.VALIDATED_MECHANISM

    def evaluate(self, claim, io):
        d = io.get("genetics", claim)
        gwas = d.get("gwas_hits", 0) or 0
        mendelian = d.get("mendelian_evidence", False)
        somatic = d.get("somatic_driver_evidence", False)
        lof = d.get("loss_of_function_phenotype")
        if not any([gwas > 0, mendelian, somatic, lof]):
            return RuleResult(
                self.rule_id, RuleStatus.PASSED, 0.4,
                caveats=[Caveat(CaveatKind.SUBSTANTIVE, self.rule_id,
                                "no genetics support for a validated-mechanism claim")],
                evidence=d,
            )
        if d.get("clinical_outcome_contested", False):
            return RuleResult(
                self.rule_id, RuleStatus.PASSED, 0.7,
                caveats=[Caveat(
                    CaveatKind.SUBSTANTIVE, self.rule_id,
                    "clinical outcome contested; framework does not predict trial outcomes - "
                    "mechanism validity does not entail clinical efficacy"
                )],
                evidence=d,
            )
        return RuleResult(self.rule_id, RuleStatus.PASSED, 0.85, evidence=d)


class R4_Expression:
    rule_id = "R4_expression"
    version = "1.0.0"
    description = (
        "Target should be detectably expressed in the indication-relevant tissue. "
        "Absence of expression is a substantive caveat; presence is supporting."
    )

    def applies_to(self, claim):
        return claim.claim_type in (ClaimType.VALIDATED_MECHANISM, ClaimType.NOVEL_TARGET)

    def evaluate(self, claim, io):
        d = io.get("expression", claim)
        if not d:
            return RuleResult(
                self.rule_id, RuleStatus.NOT_APPLICABLE, 0.0,
                caveats=[Caveat(CaveatKind.OPERATIONAL, self.rule_id,
                                "no expression data available; rule not applicable")],
            )
        expressed = d.get("target_tissue_expressed", None)
        if expressed is False:
            return RuleResult(
                self.rule_id, RuleStatus.FALSIFIED, 0.7,
                FalsificationTier.PUBLIC_DATA_LOOKUP,
                "Query GTEx / Human Protein Atlas / relevant scRNA-seq atlas. Target not detected in indication-relevant tissue.",
                d,
            )
        return RuleResult(self.rule_id, RuleStatus.PASSED, 0.8, evidence=d)


class R5_Replication:
    rule_id = "R5_replication"
    version = "1.0.0"
    description = (
        "The primary claim must not be retracted, and must not have multiple "
        "formal rebuttals without successful independent replication. "
        "Failure here is among the strongest possible falsifications from public data."
    )

    def applies_to(self, claim):
        return True

    def evaluate(self, claim, io):
        d = io.get("reproducibility", claim)
        if not d:
            return RuleResult(
                self.rule_id, RuleStatus.ABSTAINED, 0.0,
                caveats=[Caveat(CaveatKind.OPERATIONAL, self.rule_id,
                                "no reproducibility data; rule abstains")],
            )
        retracted = d.get("retracted", False)
        rebuttals = d.get("rebuttals_count", 0) or 0
        replications = d.get("independent_replications", 0) or 0
        if retracted:
            year = d.get("retraction_year", "?")
            return RuleResult(
                self.rule_id, RuleStatus.FALSIFIED, 0.99,
                FalsificationTier.PUBLIC_DATA_LOOKUP,
                f"PubMed / Retraction Watch lookup. Primary claim is retracted ({year}).",
                d,
            )
        if rebuttals >= 2 and replications == 0:
            return RuleResult(
                self.rule_id, RuleStatus.FALSIFIED, 0.95,
                FalsificationTier.PUBLIC_DATA_LOOKUP,
                f"{rebuttals} formal peer-reviewed rebuttals with zero independent replications. Public-record disproof.",
                d,
            )
        return RuleResult(self.rule_id, RuleStatus.PASSED, 0.9, evidence=d)


class R6_ChemistryClassCollapse:
    rule_id = "R6_chemistry_class_collapse"
    version = "1.2.0"
    description = (
        "If ChEMBL hits for the claimed target collapse onto a single human "
        "Pfam class (e.g. >=80% are inhibitors of a known human paralog), the "
        "chemistry evidence is phantom - it's evidence for the paralog, not "
        "the claimed target. This is the SAT->HDAC4 archetype from the Madurella audit. "
        "v1.2.0 upgrade: scope expanded to validated_mechanism claims, because "
        "the pan-class selectivity question (canonical example: pan-JAK overlap "
        "for TYK2) is exactly the structural risk this rule was written to catch. "
        "In addition to the existing exact-fraction falsification path, v1.2.0 "
        "adds a heuristic substantive caveat when the paralog ChEMBL compound "
        "pool is >= 2x the primary target's: a chemistry community that has "
        "invested more in the paralogs than in the primary is unlikely to have "
        "produced a class-selective story without explicit structural justification."
    )

    def applies_to(self, claim):
        return claim.claim_type in (
            ClaimType.NOVEL_TARGET,
            ClaimType.CHEMISTRY_SERIES,
            ClaimType.VALIDATED_MECHANISM,
        )

    def evaluate(self, claim, io):
        d = io.get("chemistry", claim)
        if not d or (d.get("chembl_distinct_compounds", 0) or 0) == 0:
            return RuleResult(
                self.rule_id, RuleStatus.NOT_APPLICABLE, 0.0,
                caveats=[Caveat(CaveatKind.OPERATIONAL, self.rule_id,
                                "no chemistry to assess for class collapse; rule not applicable")],
            )

        # Path 1: exact class-collapse fraction (when available, e.g. via
        # explicit compound-level join). Unchanged from v1.0.0.
        fraction = d.get("chembl_pfam_class_collapse_fraction", None)
        target = d.get("chembl_pfam_class_collapse_target_symbol", "<unknown>")
        if fraction is not None and fraction >= 0.80:
            pct = f"{fraction:.0%}"
            return RuleResult(
                self.rule_id, RuleStatus.FALSIFIED, 0.95,
                FalsificationTier.PUBLIC_DATA_LOOKUP,
                f"Group ChEMBL hits by Pfam class. {pct} are {target} inhibitors - the chemistry support for {claim.target_symbol} is class-collapsed phantom evidence.",
                {"collapse_fraction": fraction, "collapse_target": target},
            )

        # Path 2 (v1.2.0): paralog-ratio heuristic. If the paralog chemistry
        # pool dwarfs the primary's, the class is much more likely to be a
        # paralog story than a selective primary story. Threshold 2.0 is a
        # heuristic, not a join.
        primary_n = d.get("chembl_distinct_compounds", 0) or 0
        paralog_counts = d.get("chembl_paralog_compound_counts", None)
        if paralog_counts and primary_n > 0:
            # paralog_counts is a dict {paralog_symbol: count}
            max_paralog = max(paralog_counts.values()) if paralog_counts else 0
            if max_paralog >= 2 * primary_n:
                # Find which paralog has the largest pool
                top_paralog = max(paralog_counts.items(), key=lambda x: x[1])
                ratio = max_paralog / primary_n
                return RuleResult(
                    self.rule_id, RuleStatus.PASSED, 0.5,
                    caveats=[Caveat(
                        CaveatKind.SUBSTANTIVE, self.rule_id,
                        f"paralog chemistry pool ({top_paralog[0]}: {top_paralog[1]} compounds) "
                        f"is {ratio:.1f}x the primary target's ({primary_n} compounds); class-collapse "
                        f"risk is structural and the selectivity narrative needs explicit support "
                        f"(e.g. distinct binding mode, allosteric site). Heuristic flag, not falsification."
                    )],
                    evidence={
                        "primary_compounds": primary_n,
                        "paralog_compound_counts": paralog_counts,
                        "max_paralog_ratio": ratio,
                    },
                )

        return RuleResult(self.rule_id, RuleStatus.PASSED, 0.85,
                          evidence={"collapse_fraction": fraction})


class R7_SelectivityCounterscreen:
    rule_id = "R7_selectivity_counterscreen"
    version = "1.1.0"
    description = (
        "For pathogen targets or human targets where indication-relevant "
        "off-target liability matters, selectivity data should exist OR the "
        "claim should explicitly flag the selectivity gap. Hitting an "
        "indication-relevant human paralog at <1-log selectivity is a "
        "computational falsification of the series. v1.1.0 upgrade: for "
        "novel-target claims, the absence of selectivity data is a "
        "substantive (not operational) gap - the selectivity question is "
        "the entire reason such claims are provisional, so no data on it "
        "must surface in the verdict."
    )

    def applies_to(self, claim):
        return claim.claim_type in (
            ClaimType.CHEMISTRY_SERIES, ClaimType.NOVEL_TARGET, ClaimType.VALIDATED_MECHANISM,
        )

    def evaluate(self, claim, io):
        d = io.get("selectivity", claim)
        has = d.get("selectivity_data", False) if d else False
        if not has:
            if claim.claim_type == ClaimType.NOVEL_TARGET:
                return RuleResult(
                    self.rule_id, RuleStatus.PASSED, 0.5,
                    caveats=[Caveat(
                        CaveatKind.SUBSTANTIVE, self.rule_id,
                        "no selectivity data for a novel-target claim; the selectivity "
                        "vs human paralog question is unresolved and the claim remains "
                        "provisional until counter-screen data exists"
                    )],
                    evidence={"selectivity_data_present": False},
                )
            return RuleResult(
                self.rule_id, RuleStatus.ABSTAINED, 0.0,
                caveats=[Caveat(CaveatKind.OPERATIONAL, self.rule_id,
                                "no selectivity data; rule abstains")],
            )
        off = d.get("off_targets_in_indication_relevant_tissue", False)
        sel = d.get("selectivity_index_log")
        if off and sel is not None and sel < 1.0:
            return RuleResult(
                self.rule_id, RuleStatus.FALSIFIED, 0.85,
                FalsificationTier.CHEAP_IN_SILICO,
                f"Counter-screen series against indication-relevant paralog. Observed log selectivity index = {sel} (threshold 1.0). Series is non-selective.",
                d,
            )
        return RuleResult(self.rule_id, RuleStatus.PASSED, 0.8, evidence=d)


RULES = [
    R1_Orthology(), R2_ChemistrySupport(), R3_GeneticsSupport(),
    R4_Expression(), R5_Replication(), R6_ChemistryClassCollapse(),
    R7_SelectivityCounterscreen(),
]


TIER_ORDER = {
    FalsificationTier.PUBLIC_DATA_LOOKUP: 0,
    FalsificationTier.CHEAP_IN_SILICO: 1,
    FalsificationTier.TARGETED_ASSAY: 2,
    FalsificationTier.COHORT_STUDY: 3,
    FalsificationTier.CLINICAL_TRIAL: 4,
}


def aggregate(results):
    applicable = [r for r in results if r.status != RuleStatus.NOT_APPLICABLE]
    evaluated = [r for r in applicable if r.status != RuleStatus.ABSTAINED]
    falsified = [r for r in results if r.status == RuleStatus.FALSIFIED]
    substantive, operational = [], []
    for r in results:
        for c in r.caveats:
            (substantive if c.is_substantive else operational).append(c)

    if len(evaluated) < 2 and not falsified:
        return Verdict.INSUFFICIENT_DATA, 0.0, None, substantive, operational

    cheapest = None
    if falsified:
        def key(r):
            return (TIER_ORDER.get(r.falsification_tier, 99), r.rule_id)
        best = min(falsified, key=key)
        cheapest = (best.falsification_experiment or "", best.falsification_tier, best.rule_id)

    score = len(falsified) / len(applicable) if applicable else 0.0

    if len(falsified) == 0:
        v = Verdict.FALSIFIED_WITH_CAVEATS if substantive else Verdict.SURVIVED
        return v, score, None, substantive, operational
    if len(falsified) >= 2:
        return Verdict.FALSIFIED, score, cheapest, substantive, operational
    only = falsified[0]
    if only.confidence >= 0.95 and only.falsification_tier == FalsificationTier.PUBLIC_DATA_LOOKUP:
        return Verdict.FALSIFIED, score, cheapest, substantive, operational
    return Verdict.FALSIFIED_WITH_CAVEATS, score, cheapest, substantive, operational


def run_audit(claim, fixture):
    io = FixtureAdapter(fixture)
    results = []
    for rule in RULES:
        if not rule.applies_to(claim):
            results.append(RuleResult(rule.rule_id, RuleStatus.NOT_APPLICABLE, 0.0))
            continue
        results.append(rule.evaluate(claim, io))
    return aggregate(results), results


def main():
    here = Path(__file__).parent if "__file__" in globals() else Path.cwd()
    sentinel_path = here / "sentinels" / "v1_sentinels.yaml"
    with open(sentinel_path) as f:
        suite = yaml.safe_load(f)

    passes = 0
    fails = 0
    lines = []
    for s in suite["sentinels"]:
        c = s["claim"]
        claim = TargetClaim(
            target_symbol=c["target_symbol"],
            indication=c["indication"],
            mechanism=c["mechanism"],
            claim_type=ClaimType(c["claim_type"]),
            uniprot_id=c.get("uniprot_id"),
        )
        (verdict, score, cheapest, _, _), _ = run_audit(claim, s["fixture"])
        expected = Verdict(s["expected_verdict"])
        ok = verdict == expected
        if ok:
            passes += 1
        else:
            fails += 1
        marker = "PASS" if ok else "FAIL"
        sid = s["id"]
        extra = ""
        if cheapest:
            extra = f" (cheapest: {cheapest[2]} @ {cheapest[1].value})"
        suffix = "" if ok else f"  EXPECTED {expected.value}"
        lines.append(f"  {marker} {sid} -> {verdict.value}{extra}{suffix}")

    total = passes + fails
    print(f"\nSentinel suite: {passes}/{total} passed\n")
    print("\n".join(lines))
    if fails == 0:
        print("\n*** v1.1.0 IS SHIPPABLE - sentinel calibration passes ***")
        return 0
    print(f"\n*** {fails} sentinel(s) failed ***")
    return 1


if __name__ == "__main__":
    sys.exit(main())
