# Backtest 2026 — audit results (B1, B2, B6)

_Audit run: 2026-05-28T01:58:55.459173+00:00_  
_Pre-registration merge commit: 264a246_  
_Audited against main HEAD: 264a246d8bbfa9a9552bff0e523a55d91fd7c719_  
_Ruleset SHA: 35ef2b2ab5363298097962a0b6ae52c70d551a1edddc341054f75cb6e4fb7221_  

## Summary

- Predictions HELD: **2 / 3**
- DIVERGED (actual verdict differs from prediction): 1

Predictions were locked in commit 
[264a246](https://github.com/crisprking/falsifiable-targets/commit/264a246) 
BEFORE this audit ran. The pre-registered YAMLs were not modified during 
the audit. The git log between this commit and the pre-registration commit 
proves the order.

## Per-claim results

### STAP cells at 2014-01-29

- **Predicted:** `FALSIFIED_WITH_CAVEATS or INSUFFICIENT_DATA` — no rule fires; R7 substantive caveat expected for extraordinary_claim with no selectivity data
- **Actual:** `SURVIVED` (exit code 0, score 0.0)
- **Match:** **DIVERGED**
- **Claim SHA:** `707a0e26495d5d65...`
- **Report:** [`studies/backtest_2026/audit_outputs/b1_stap_cells_2014_at_submission.json`](studies/backtest_2026/audit_outputs/b1_stap_cells_2014_at_submission.json)
- **All applicable rules cleared.**

### CETP/torcetrapib at 2005-08-24

- **Predicted:** `FALSIFIED_WITH_CAVEATS` — R3 substantive caveat via clinical_outcome_contested
- **Actual:** `FALSIFIED_WITH_CAVEATS` (exit code 1, score 0.0)
- **Match:** **HELD**
- **Claim SHA:** `1e904be887969649...`
- **Report:** [`studies/backtest_2026/audit_outputs/b2_cetp_torcetrapib_2005_at_submission.json`](studies/backtest_2026/audit_outputs/b2_cetp_torcetrapib_2005_at_submission.json)
- **Substantive caveats:**
  - `R3_genetics_support`: clinical outcome contested; framework does not predict trial outcomes - mechanism validity does not entail clinical efficacy
- **Rules abstained (no data):** R2_chemistry_support

### BACE1 at 2012-11-01

- **Predicted:** `SURVIVED` — no rule fires — calibration anchor
- **Actual:** `SURVIVED` (exit code 0, score 0.0)
- **Match:** **HELD**
- **Claim SHA:** `3b40faefe54eff83...`
- **Report:** [`studies/backtest_2026/audit_outputs/b6_bace1_alzheimers_2012_at_submission.json`](studies/backtest_2026/audit_outputs/b6_bace1_alzheimers_2012_at_submission.json)
- **Rules abstained (no data):** R2_chemistry_support

## What these results mean for the article

**B6 (BACE1@2012) held.** The framework did not predict the future. 
Given evidence as of November 2012 — APP A673T published 3 months 
earlier, zero failed Phase 3 trials — the framework correctly verdicts 
`SURVIVED`. Paired with the as-of-today BACE1 YAML which correctly 
demotes the same target to `FALSIFIED_WITH_CAVEATS`, this is the 
trajectory-tracking story: the framework's verdict moves as evidence 
accumulates; it does not predict trial outcomes.

**B2 (CETP@2005) held.** This is the killer demonstration: using only 
evidence published before ILLUMINATE Phase 3 began (Boekholdt 2005 and 
de Grooth 2004, both null on CV outcome for the TaqIB variant despite 
HDL effect), the framework would have raised a substantive caveat that 
named the exact failure mode that materialized 15 months later.

**B1 (STAP@submission) diverged: predicted FWC/INSUFFICIENT_DATA, got `SURVIVED`.** 
Investigation note below.

## Divergence investigation (per the pre-registered methodology)

The pre-registration committed: "a claim where the audit's verdict 
clearly disagrees with the historical outcome" gets reported as-is, 
then in a footnote we note any evidence missed and what a re-audit 
would say. The YAML for the diverged case is NOT modified in this commit.

### STAP cells at 2014-01-29 — DIVERGED

Predicted `FALSIFIED_WITH_CAVEATS or INSUFFICIENT_DATA`; got `SURVIVED`.

Possible explanations (to investigate, not commit-now):

- The YAML field semantics may differ from my reading of the rule code.
- The at-submission evidence may have been over- or under-specified.
- The framework may simply behave differently than predicted on this 
  claim type, which is itself a finding worth publishing.

Next step: read the rule that should have fired, compare against the 
YAML field values, and either (a) confirm the framework's behavior is 
correct given the inputs and update the article, or (b) propose a 
YAML adjustment in a separate, clearly-marked commit on a separate 
branch. NEVER silently retune the pre-registered YAML.

## Reproducibility

Anyone can reproduce these results:

```
git clone https://github.com/crisprking/falsifiable-targets.git
cd falsifiable-targets
pip install -e '.[all]'
for f in studies/backtest_2026/claims/*.yaml; do
  python run_audit.py "$f" --no-live --json-out "/tmp/$(basename $f .yaml).json"
done
```

Ruleset SHA `35ef2b2ab5363298097962a0b6ae52c7...` is content-
addressed; the same audit on any machine produces the same JSON modulo 
timestamp and tool.python_version fields.
