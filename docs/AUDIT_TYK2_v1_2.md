# TYK2 / psoriasis audit

**Claim**: TYK2 is a validated therapeutic target for moderate-to-severe
plaque psoriasis. Selective allosteric TYK2 inhibition (via the JH2
pseudokinase domain) suppresses the IL-23/Th17 axis that drives
psoriatic inflammation. Deucravacitinib (BMS-986165, Sotyktu) received
FDA approval in September 2022 on the basis of the POETYK PSO-1 and
PSO-2 Phase 3 trials.

**Verdict**: `SURVIVED`
**Score**: 0.00 (0 falsified / 5 applicable)
**Substantive caveats**: 0
**Operational notes**: 0

**Ruleset**: v1.1.0 (SHA `2f9aab7d0ebc209f62c16eb35be31bc5b65fa2eb09adc02bea5bff5176269b32`)
**Claim SHA**: `b663fa07cb1038ac...` (full hash in JSON report)
**Audit date (UTC)**: 2026-05-25

The complete machine-readable record is at
[`../reports/tyk2_psoriasis_audit.json`](../reports/tyk2_psoriasis_audit.json).
The claim YAML with all encoded inputs and primary-literature citations
is at [`../claims/tyk2_psoriasis.yaml`](../claims/tyk2_psoriasis.yaml).

## Why TYK2 was chosen

TYK2 is the first non-self audit performed by this framework. It was
chosen with three criteria:

1. **Validated mechanism with FDA approval.** Deucravacitinib's 2022
   approval makes the predicted verdict (`SURVIVED`) externally
   verifiable. If the framework returned anything other than
   `SURVIVED`, the framework would be the thing under scrutiny, not
   TYK2.
2. **Live-adapter exercise.** Both UniProt (P29597) and ChEMBL
   (CHEMBL3553) contain rich live data for TYK2. The audit therefore
   stress-tests the live read-side adapters introduced in v1.0.1.
3. **Known calibration risk.** The pan-JAK selectivity question - the
   reason deucravacitinib needed JH2-domain binding to be clinically
   viable - is the textbook class-collapse risk for a kinase target.
   It is exactly the kind of failure mode R6 was written to catch.
   The audit was set up to find out whether R6 fires on this case.

The third criterion turned out to be the audit's most interesting
finding. See §3 below.

## Per-rule findings

### R1_orthology - `NOT_APPLICABLE`

TYK2 is a `validated_mechanism` claim. R1 only fires on `novel_target`
and `extraordinary_claim` claims, where the question "does this protein
even exist outside the claimant's preferred lineage" is live. For a
human protein with a 2022 FDA approval, R1 does not apply by design.

### R2_chemistry_support - `PASSED` (confidence 0.9)

Live ChEMBL fetched 538 distinct molecules with bioactivity against
CHEMBL3553 (TYK2). The rule's threshold is `n >= 1`; 538 substantially
clears it. Caveat: my adapter fetches the first 1000 activity records
(pagination is a v1.3 improvement), so 538 is a lower bound on the
distinct-compound count rather than the full enumeration.

### R3_genetics_support - `PASSED` (confidence 0.85)

R3's required inputs (`gwas_hits`, `mendelian_evidence`, etc.) come
from the fixture, not live data. The values encoded:

- 24 GWAS hits in IL-23/Th17 axis disorders (TYK2 is among the
  strongest GWAS signals for psoriasis)
- Mendelian evidence: complete TYK2 deficiency produces hyper-IgE
  syndrome with mycobacterial susceptibility (Minegishi et al.,
  *Immunity* 2006); partial loss-of-function via the P1104A variant
  is protective against psoriasis
- `clinical_outcome_contested = false` - Phase 3 readouts and FDA
  approval are the standard of evidence

R3 passes without caveat. **Reader weight**: This is fixture data,
not live evidence. v1.3.0 will add an Open Targets adapter to make
this a live-evidence rule.

### R4_expression - `PASSED` (confidence 0.8)

Live UniProt provided a tissue-specificity annotation:
"Observed in all cell lines analyzed. Expressed in a variety of
lymphoid and non-lymphoid cell lines." The fixture independently
supplies `target_tissue_expressed = true`.

**Reader weight**: R4's live check is structurally weak - it tests
"does UniProt have any expression annotation," not "is the target
expressed in lesional psoriatic skin specifically." For TYK2 the
broader expression is sufficient for the rule to pass, but the rule
is not semantically auditing the indication-specific tissue match.
A more rigorous check (GTEx by tissue, or scRNA-seq of psoriatic skin)
is v1.3+.

### R5_replication - `PASSED` (confidence 0.9)

Fixture inputs: `retracted = false`, `rebuttals_count = 0`,
`independent_replications = 50` (POETYK PSO-1, POETYK PSO-2, multiple
mechanistic replications across independent labs).

**Reader weight**: No live retraction-check adapter. The fixture
asserts no retraction; a retraction occurring after the audit date
would not be caught.

### R6_chemistry_class_collapse - `NOT_APPLICABLE`

**This is the audit's most important structural finding.** R6 was
written to catch the failure mode where ChEMBL compounds claimed as
evidence for a target turn out to be inhibitors of a paralog in the
same Pfam class - the SAT/HDAC4 archetype.

Under ruleset v1.1.0, R6's `applies_to()` is restricted to
`novel_target` and `chemistry_series` claims. TYK2 is
`validated_mechanism`, so R6 returns `NOT_APPLICABLE` regardless of
the actual class-collapse fraction in the data.

This is a real gap. The pan-JAK selectivity question is the canonical
class-collapse risk for kinase-family validated mechanisms.
Deucravacitinib's entire scientific story is that it solves the
class-collapse problem (by binding JH2 instead of the active site, a
binding mode JAK1/2/3 do not share). The framework's job *should*
include checking this. Under v1.1.0 it does not.

**Status**: queued as the headline change for v1.3.0. See
[`AUDIT_LIMITATIONS_v1_2.md`](AUDIT_LIMITATIONS_v1_2.md) §2 for the
v1.3.0 implementation plan.

### R7_selectivity_counterscreen - `PASSED` (confidence 0.8)

Fixture inputs:
- `selectivity_data = true`
- `off_targets_in_indication_relevant_tissue = false`
- `selectivity_index_log = 2.0` (sourced from Wrobleski et al.,
  *J Med Chem* 2019, reporting >100-fold cellular selectivity vs
  JAK1/2/3)

The R7 threshold is `selectivity_index_log >= 1.0`. 2.0 is two log
units above threshold.

**Reader weight**: This is fixture data, encoded from primary
literature. A live kinome-selectivity adapter is v1.4+.

## What this audit does and does not prove

**Proves:**

- Live UniProt + ChEMBL adapters work end-to-end against a real,
  well-characterised target on Kaggle's network environment
- The v1.1.0 ruleset is not biased toward `FALSIFIED_WITH_CAVEATS` -
  a fully-validated mechanism cleanly passes with zero substantive
  caveats
- The framework distinguishes between the in-house novel-target case
  (Ipi1, 2 substantive caveats) and a validated external mechanism
  (TYK2, 0 substantive caveats) using the same ruleset and engine
- The SHA-locked report format is reproducible: re-running the audit
  with `AE_OFFLINE=1` against the archived cache will produce
  byte-identical verdict and per-rule status

**Does not prove:**

- That TYK2 has been audited on its most-relevant axis. R6 (the
  class-collapse rule, which is the natural place to assess pan-JAK
  overlap) was `NOT_APPLICABLE` under v1.1.0 by design. The
  `SURVIVED` verdict is correct under v1.1.0 logic but does not
  reflect a check of the most class-specific risk for a kinase target.
- That the genetics, replication, and selectivity verdicts reflect
  independent live evidence. They reflect fixture values paraphrased
  from primary literature, not data the framework fetched and
  parsed.
- That the framework would catch a covert TYK2 retraction or a new
  contradicting GWAS finding - those rules use fixture data only
  under v1.2.x.

## Comparison with the Ipi1 inaugural audit

The framework's two completed audits, side by side:

| Axis | Ipi1 (v1.1.0) | TYK2 (v1.1.0) |
|---|---|---|
| Claim type | `novel_target` | `validated_mechanism` |
| Verdict | `FALSIFIED_WITH_CAVEATS` | `SURVIVED` |
| Substantive caveats | 2 (R1 orthology, R7 selectivity) | 0 |
| Live-data rules | R5 (no retraction) | R2 (ChEMBL count), R4 (UniProt tissue) |
| Fixture-data rules | R1, R7 (the load-bearing rules for this claim) | R3, R5, R7 |
| Not-applicable rules | R2, R3, R4, R6 | R1, R6 |

The asymmetry is the calibration result. Same engine, same SHA-locked
rule logic, different verdicts in the direction the structure of the
two claims predicts.

## Next steps

1. **v1.3.0**: expand R6 to `validated_mechanism`, implement live
   class-collapse computation across JAK family, re-audit TYK2.
   Whatever the new verdict, it will be reported honestly.
2. **v1.3.0**: live Open Targets adapter for R3, converting the
   genetics verdict from fixture-paraphrased to live-fetched.
3. **v1.3.0**: live Europe PMC retraction-check for R5.

The TYK2 audit is the v1.2.x calibration milestone, not the project's
endpoint.
