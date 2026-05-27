# TYK2 / psoriasis audit — v1.3.0 re-audit

**Claim**: TYK2 is a validated therapeutic target for moderate-to-severe
plaque psoriasis. Selective allosteric TYK2 inhibition (via the JH2
pseudokinase domain) suppresses the IL-23/Th17 axis. Deucravacitinib
(BMS-986165, Sotyktu) received FDA approval September 2022 (POETYK
PSO-1 and PSO-2 Phase 3 trials).

**Verdict**: `SURVIVED`
**Score**: 0.00 (0 falsified / **6** applicable — was 5 under v1.1.0)
**Substantive caveats**: 0
**Operational notes**: 0

**Ruleset**: v1.2.0 (SHA `35ef2b2ab5363298097962a0b6ae52c70d551a1edddc341054f75cb6e4fb7221`)
**Claim SHA**: `b663fa07cb1038ac…` (unchanged from v1.2.0 audit)
**Audit date (UTC)**: 2026-05-26

The complete machine-readable record is at
[`../reports/tyk2_psoriasis_audit.json`](../reports/tyk2_psoriasis_audit.json).
The previous v1.1.0 audit narrative remains archived at
[`AUDIT_TYK2_v1_2.md`](AUDIT_TYK2_v1_2.md).

## Why this audit was re-run

The v1.2.x TYK2 audit ([`AUDIT_TYK2_v1_2.md`](AUDIT_TYK2_v1_2.md))
returned `SURVIVED` with R6 `NOT_APPLICABLE`. The
[`AUDIT_LIMITATIONS_v1_2.md`](AUDIT_LIMITATIONS_v1_2.md) document
named this as the audit's most important structural finding: R6's
scope didn't include `validated_mechanism` claims, so the canonical
class-collapse risk for kinase targets (pan-JAK overlap) was not
audited. v1.3.0 closes that scope gap. This re-audit was the
calibration test of the new R6 logic against the real, FDA-approved
mechanism.

The v1.3.0 prediction, written before running the audit: TYK2 would
likely shift to `FALSIFIED_WITH_CAVEATS` because the JAK family
chemistry pools were assumed to dwarf TYK2's.

**The prediction was wrong.** The verdict didn't shift. This document
records why.

## What R6 actually saw

R6 v1.2.0 fired the heuristic substantive-caveat path. The path
requires `chembl_paralog_compound_counts` in the chemistry section,
which the v1.3.0 ChEMBLAdapter populates by fetching paralog compound
counts when `claims/paralog_map.yaml` has an entry for the primary's
UniProt accession.

For TYK2 (P29597), the paralog map entry is:

```yaml
P29597:
  primary_symbol: TYK2
  paralogs:
    - {symbol: JAK1, uniprot_id: P23458}
    - {symbol: JAK2, uniprot_id: O60674}
    - {symbol: JAK3, uniprot_id: P52333}
```

The live ChEMBL adapter fetched compound counts for all four kinases.
Result:

| Kinase | UniProt | ChEMBL target ID | Distinct compounds |
|---|---|---|---|
| TYK2 (primary) | P29597 | CHEMBL3553 | **538** |
| JAK1 (paralog) | P23458 | CHEMBL2835 | 549 |
| JAK2 (paralog) | O60674 | CHEMBL2971 | 772 |
| JAK3 (paralog) | P52333 | CHEMBL2148 | 832 |

The maximum paralog/primary ratio is **832 / 538 = 1.55×**.

R6's heuristic threshold for the paralog-overshadow substantive caveat
is **2.0×**. 1.55× is below threshold. R6 returns `PASSED` with no
caveat and 0.85 confidence. The pre-existing exact-fraction
falsification path is also unused (the fixture supplies `null` for the
exact class-collapse fraction, which is correct under v1.2.x — exact
overlap computation is a v1.4+ feature).

## Interpretation: what the data is telling us

A prior I held going into this audit: "the JAK family is the canonical
class-collapse risk for kinase chemistry; the paralog pools surely
dwarf TYK2's by 5× or more." The data disagreed.

Three possible explanations, in decreasing order of likely correctness:

1. **TYK2 is a genuinely well-studied kinase relative to its
   paralogs.** Of the four JAKs, JAK2 and JAK3 each have inherited
   substantial inhibitor programs (myelofibrosis: ruxolitinib for
   JAK1/2; alopecia and JAK3-selective programs). TYK2-selective
   chemistry is a recent but extensive effort driven by
   deucravacitinib's commercial success and the IL-23 axis interest
   in autoimmune disease. The 538-compound TYK2 pool is *not small*
   relative to JAK1's 549 or JAK3's 832. The field has invested
   comparably across the family.

2. **The adapter's `limit=1000` per fetch caps each kinase at up to
   1000 records.** This caps the *upper bound* of distinct compounds
   per kinase the adapter can count. TYK2 came back with 538 unique
   compounds from 1000 activity records (heavy multi-record reuse
   per compound); JAK3 came back with 832 (less reuse). If JAK1 or
   JAK2 had, say, 4000 unique compounds in ChEMBL but the adapter
   only saw 1000 activity records, the count would be capped. But
   the ratios between observed counts are still informative —
   they reflect the top-activity-density-1000-records signal, which
   is what's tractable without pagination. A future v1.4 with
   pagination would refine the absolute counts but probably not the
   ratio shape.

3. **The 2.0× threshold is poorly calibrated.** This would mean real
   class-collapse cases exist at, say, 1.5×-2.0× and the heuristic
   misses them. To know this would require auditing several more
   kinase-class validated mechanisms and looking at the distribution
   of ratios. v1.4 work, not a v1.3 patch.

The honest reading is some mixture of (1) and (2). The framework's
heuristic gave the right answer for what it actually measures:
*pool-size overshadow*. The pan-JAK selectivity concern that
deucravacitinib was engineered to solve operates at a different
axis — *binding-mode overlap*, not *pool-size overshadow*. A small
molecule that hits TYK2 active site is almost certainly going to hit
JAK1/2/3 active sites at similar potency; the JH2-domain mechanism is
the structural escape from that overlap.

R6 v1.2.0 measures pool-size overshadow because that's what's
tractable from public chemistry counts. The framework does not yet
audit compound-level overlap (a v1.4+ feature). What it can say from
its current input is: TYK2's chemistry pool is *not* dwarfed by its
paralogs. That is real, evidence-based information. It is not a full
selectivity audit — and the framework should not pretend it is.

## What this audit does and does not prove

**Proves:**

- R6 v1.2.0's scope expansion is operative on a `validated_mechanism`
  claim. The rule applies, the rule evaluates, the rule returns a
  PASSED status with input data attached to the report.
- The live ChEMBL paralog-count fetch works end-to-end on Kaggle for
  a real four-target paralog group.
- Under the framework's actual input — ChEMBL distinct-compound counts
  per UniProt accession, capped at the adapter's 1000-record fetch
  limit — TYK2 does not exhibit a 2×-overshadow class-collapse signal.
- The v1.3.0 SURVIVED verdict for TYK2 is therefore *qualitatively
  stronger* than the v1.2.x SURVIVED: the v1.2.x verdict was silent
  on the pan-class question; v1.3.0 verdicts confirm the framework
  asked it and found no overshadow at the pool-size resolution.

**Does not prove:**

- That TYK2 chemistry would not exhibit compound-level overlap with
  JAK1/2/3 if the join were actually computed. It almost certainly
  would. Active-site TYK2 inhibitors hit pan-JAK by structural
  necessity; deucravacitinib's value is that it does not bind the
  active site.
- That R6's 2.0× threshold is the right calibration. The threshold
  was set before this audit and not adjusted after. Whether it
  catches *more subtle* class-collapse cases at, e.g., 1.5× is
  unknown until more kinase-class targets are audited.
- That the framework would catch deucravacitinib's predecessors that
  failed for pan-JAK toxicity. Tofacitinib (JAK1/3) and ruxolitinib
  (JAK1/2) are pan-JAK by design — they would all show paralog
  ratios near 1.0 (each is approximately equally well-developed
  across its targeted set). R6's pool-size heuristic would not
  flag them. The pool-size and overlap axes are different things,
  and only the latter would distinguish them under this rule.

## Comparison: v1.1.0 vs v1.3.0 TYK2 verdict

| | v1.1.0 | v1.3.0 |
|---|---|---|
| Verdict | `SURVIVED` | `SURVIVED` |
| Substantive caveats | 0 | 0 |
| Rules applicable | 5 | **6** (R6 now in scope) |
| R6 status | `NOT_APPLICABLE` (scope) | `PASSED` (heuristic checked, ratio 1.55×) |
| R6 input data | (rule not invoked) | `{primary_compounds: 538, paralog_counts: {JAK1: 549, JAK2: 772, JAK3: 832}}` |
| Reader weight | "framework didn't check pan-JAK" | "framework checked pan-JAK by pool-size; no overshadow" |

The verdict is the same. The *evidence behind the verdict* is
strictly stronger. That is the v1.3.0 calibration result for the TYK2
audit: not a verdict shift, but an evidence depth increase.

## Honest v1.4 questions

This audit surfaces three v1.4-worthy questions:

1. **Compound-level overlap.** Currently R6 measures pool sizes. A
   compound-level overlap computation (for each TYK2-active compound,
   what fraction also hits a JAK paralog with comparable potency?)
   would be a substantively different rule. The deucravacitinib
   case is precisely where this would matter: pool sizes are
   comparable, but binding-mode-level overlap (active-site competition
   vs. JH2 allosteric) is what differentiates selective from pan-JAK
   compounds.

2. **Pagination beyond 1000 records.** The current 1000-record cap
   per ChEMBL fetch is a known limitation. Whether pagination would
   change ratios meaningfully is unknown until measured.

3. **Threshold calibration distribution.** With only one external
   audit (TYK2) at the pool-ratio axis, we don't know whether 2.0× is
   too coarse, just right, or too aggressive. Auditing GPR75/obesity,
   NLRP3/inflammasome, or other validated mechanisms with named
   paralog risks would expand the calibration sample.

These are real v1.4 priorities, not v1.3 patches. The v1.3.0 result
stands as the published audit for now.

## Reproduce locally

```bash
git clone https://github.com/crisprking/falsifiable-targets.git
cd falsifiable-targets
git checkout v1.3.0
python run_audit.py claims/tyk2_psoriasis.yaml \
    --json-out reports/my_tyk2_v1_3_0.json
```

The result will match this document's verdict and per-rule status as
long as the ruleset SHA is `35ef2b2a…` (the v1.2.0 lock). The exact
paralog counts may drift as ChEMBL adds records; the cached version
in `.ae_cache/` from the audit date freezes the numbers reported here.
