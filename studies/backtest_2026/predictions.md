# Locked predictions — backtest 2026

This document records verdict predictions made BEFORE any audit runs.
This commit contains the predictions but no audit results. The audit is
a separate commit in a separate PR; the git log proves the order.

## Conventions

- "Locked at" = the SHA of the commit that introduces this file. Predictions
  made before this commit's timestamp are pre-registered; predictions made
  after are not.
- "Confidence" is informal and self-reported. It is not a calibrated
  probability.

## Predictions

### B1: STAP cells at submission (2014-01-29)

- **Predicted verdict:** FALSIFIED_WITH_CAVEATS or INSUFFICIENT_DATA.
- **Predicted cheapest falsifier:** none — no rule fires. R5 cannot fire
  (no retraction yet); R2 abstains (no chemistry); R7 likely raises a
  substantive caveat for an extraordinary_claim with no selectivity data.
- **Confidence:** medium.
- **Why this matters:** if the framework hedges from day one before any
  retraction, the article story becomes "the framework flagged STAP as
  under-evidenced at submission, not in retrospect." A stronger claim than
  "the framework caught a known fraud."

### B2: CETP/torcetrapib at submission (2005-08-24)

- **Predicted verdict:** FALSIFIED_WITH_CAVEATS via R3.
- **Predicted cheapest falsifier:** R3 caveat — clinical_outcome_contested,
  given conflicting human genetics in 2004-2005: Boekholdt 2005 (IPD meta-analysis, n=13,677) found the TaqIB variant ASSOCIATED with reduced CAD risk via HDL; de Grooth 2004 (CARE cohort) found the variant did NOT predict cardiovascular events. The conflict is what justifies clinical_outcome_contested.
- **Confidence:** medium-high.
- **Why this matters:** if this holds, the framework would have flagged a
  Phase 3 program pharma later spent ~$1B on, using only evidence available
  the day Phase 3 was registered. This is the killer demonstration if it
  fires. If it doesn't fire, that itself is informative — the YAML may be
  under-specified for this case.

### B6: BACE1 at submission (2012-11-01)

- **Predicted verdict:** SURVIVED.
- **Predicted cheapest falsifier:** none expected — no rule should fire.
- **Confidence:** high.
- **Why this matters:** this is the calibration anchor. The framework
  must NOT predict the future. Given APP A673T (just published) and zero
  failed Phase 3 trials at Nov 2012, the claim was well-evidenced. If the
  framework correctly says SURVIVED at this date, and the existing
  as-of-2026 BACE1 YAML correctly demotes the same target to
  FALSIFIED_WITH_CAVEATS, that is the trajectory-tracking story:
  "the framework's verdict moves as evidence accumulates; it does not
  predict trial outcomes."

## What a "good result" looks like for this study

The framework is the kind of tool the article claims it is if:

- B1 hedges (not SURVIVED, not FALSIFIED) — honest about ambiguity.
- B2 raises a substantive caveat in 2005 — flags contested Mendelian
  evidence before the trial failure.
- B6 SURVIVES in 2012 and DEMOTES to FALSIFIED_WITH_CAVEATS as-of-today.

If any of the three diverges, that is itself a finding. We report it
honestly and ask why.
