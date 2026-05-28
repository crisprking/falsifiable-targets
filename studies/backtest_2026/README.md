# Backtest 2026 — at-submission retrospective audit

**Status (this commit):** pre-registration. Predictions locked, NO audit run yet.

## What this study is

A blind backtest of the falsifiable-targets engine against historical
target-validation claims, using only evidence published before each
claim's "submission date" (defined per-claim below).

The hypothesis: the framework can distinguish well-evidenced targets from
under-evidenced ones at the time the claim was made, NOT in retrospect
after outcomes are known.

This is NOT a validation study — n is too small and selection isn't
random. It is a transparent backtest where the methodology is committed
to the repo before any audit runs.

## Methodology

1. **Claim selection.** 10 historical target-validation claims with
   unambiguous outcomes known by 2026. Mix of approvals, failures,
   retractions, and contested cases. Full list below.

2. **As-of-date discipline.** Each YAML reflects evidence published at or
   before a stated cutoff date. PMIDs are cited inline. Evidence that
   emerged after the cutoff date is excluded, even if it later changed
   the field's consensus.

3. **Pre-registered predictions.** A locked verdict prediction for each
   claim, written before any audit runs. Committed to this repo in
   `predictions.md` in the same commit as the YAML, in a separate commit
   from the audit run.

4. **Honest failure handling.** Every claim is reported regardless of
   outcome. A claim where the audit disagrees with the historical result
   is reported as-is, then in a footnote we note any evidence that may
   have been missed and what a re-audit would say. No silent retries.

5. **Separation of pre-registration and execution.** The pre-registration
   commit contains predictions but no audit results. The audit commit
   runs the engine and records results. The git log proves the order.

## Claims in this study (this commit covers B1, B2, B6)

| # | Target | Indication | As-of date | Reason for inclusion |
|---|---|---|---|---|
| B1 | OCT4 / STAP | Pluripotent stem cells | 2014-01-29 | Retraction case; tests R5 path when no replication exists yet |
| B2 | CETP / torcetrapib | CV disease (statin add-on) | 2005-08-24 | Contested Mendelian evidence at submission; killer test if framework flags |
| B3 | PCSK9 | Hypercholesterolemia | TBD | Validated success, R3 strong |
| B4 | HMGCR / statins | Hypercholesterolemia | TBD | Baseline validated success |
| B5 | Anti-amyloid mAbs (aducanumab) | Alzheimer's | TBD | Contested mechanism, controversial approval |
| B6 | BACE1 / verubecestat | Alzheimer's | 2012-11-01 | Mechanism-valid clinical-failed; tests trajectory tracking |
| B7 | GFAJ-1 arsenic life | (extraordinary claim) | 2010-12-02 | Extraordinary_claim path |
| B8 | CB1 / rimonabant | Obesity | TBD | On-target tox case |
| B9 | Anti-amyloid (lecanemab) | Alzheimer's | TBD | Modest-success recent case |
| B10 | PAR-1 / vorapaxar | Cardiovascular | 2007 (TRACER start) | Blind case — outcome was unknown to author until literature search |

This commit ships B1, B2, B6 with full at-submission YAMLs and locked
predictions. The remaining 7 will be added in a follow-up
pre-registration commit before any audit runs.

## Limitations

- I cannot reconstruct exactly what was "known" at the submission date.
  The defensible approximation is "evidence indexed in PubMed at or before
  the as-of date." Evidence is cited by PMID in each YAML comment.
- I am one encoder. There is no independent double-encoding for
  inter-rater reliability.
- Sample size (n=10) is too small to claim statistical significance over
  a base rate. The study is a demonstration of methodology applied to
  named cases, not a validation paper.
- B10 is the only fully outcome-blind case. The other 9 are at-submission
  blind (I didn't use post-submission evidence) but I knew the outcomes
  going in. This is documented per-claim.
