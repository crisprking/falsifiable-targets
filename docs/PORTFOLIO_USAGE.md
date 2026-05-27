# Portfolio audit usage guide

_Generated 2026-05-27T16:35:29.655379+00:00_

## What this is

A workflow for batch-auditing drug-discovery target claims against the
falsifiable-targets v1.4.x framework. Turns a directory of claim YAMLs
into a triaged decision sheet with verdict, risk tier, fired rules,
substantive caveats, and next-experiment recommendations per claim.

## Output contract

Each run produces three files under `reports/`:

- `portfolio_decision.csv`: for spreadsheets and program managers.
- `portfolio_decision.md`: for PR descriptions and wiki pages, sorted
  by risk tier with highest-risk claims first.
- `portfolio_manifest.json`: for downstream automation, including
  per-claim rows, tier and verdict summaries, and the run timestamp.

## Risk tier semantics

| Tier | Source verdict | What to do |
|---|---|---|
| `TIER_1_DROP` | FALSIFIED | Rule fired with concrete public-data evidence. Cut from portfolio unless contested with new evidence. |
| `TIER_2_INVESTIGATE` | FALSIFIED_WITH_CAVEATS or INSUFFICIENT_DATA | Substantive caveats present; mechanism may be valid but concerns must be addressed before committing budget. |
| `TIER_3_PROCEED` | SURVIVED | Passes all applicable rules. Standard due diligence still applies. |
| `TIER_0_FIX_INPUT` | ERROR | YAML failed schema validation. Fix the input file and re-run. |

## Reading the action_items column

Priority order for what the column surfaces:

1. **Fired rules**: `R5_replication FALSIFIED @ public_data_lookup: <experiment>`.
   Which rule failed, at what evidence tier, what experiment resolves it.
2. **Substantive caveats**: `R3_genetics_support: clinical outcome contested; ...`.
   Rule did not fail outright but raised a concern that demoted the verdict.
3. **Data gaps**: `data gap: R2_chemistry_support abstained (no input)`.
   No data to evaluate. Fill the gap and re-audit.
4. **All clear**: `all rules cleared`. Passes every applicable rule.

## Exit-code contract

```
0 - SURVIVED (or portfolio: all TIER_3_PROCEED)
1 - FALSIFIED_WITH_CAVEATS (or portfolio: TIER_2_INVESTIGATE present)
2 - FALSIFIED (or portfolio: TIER_1_DROP present)
3 - INSUFFICIENT_DATA
5 - ERROR (missing/empty/malformed claim file or unhandled exception)
```

Treat exit 5 differently from exit 1 or 2 in CI: 1 and 2 are scientific
findings; 5 indicates a broken input file that needs to be fixed before
the audit is meaningful.

## CLI usage

```
python run_audit.py claims/my_target.yaml --no-live --json-out reports/my_target.json
python run_audit.py claims/my_target.yaml --debug   # full traceback on error
python scripts/portfolio_audit.py claims/ reports/  # batch
```

## Required claim YAML structure

```
claim:
  target_symbol: "PCSK9"                # required
  uniprot_id: "Q8NBP7"                  # required
  indication: "Hypercholesterolemia"    # required
  mechanism: |                          # required, multiline OK
    Mechanism description.
  claim_type: validated_mechanism       # required; one of:
                                        #   validated_mechanism,
                                        #   novel_target,
                                        #   extraordinary_claim,
                                        #   chemistry_series,
                                        #   repurposing
  citation_pmid: "16554528"             # optional but recommended
fixture:                                # optional offline test data
  orthology: {sources_agreeing: 4, sources_total: 4}
  genetics:  {gwas_hits: 31, mendelian_evidence: true}
```

## Verified contracts (acceptance gate)

Six contracts, all passing:

1. Verdict + exit-code mapping across five reference scenarios.
2. Error paths exit 5 with actionable stderr; no tracebacks bleed out.
3. Ruleset SHA pinned across all reports.
4. Determinism: two back-to-back runs produce byte-identical reports.
5. `--debug` flag accepted by argparse, preserves exit codes, adds
   tracebacks on exceptions that escape `main()`.
6. Substantive caveat text propagated into portfolio output with
   non-empty text per caveat.

## Workflow recommendations

**Program-level triage.** Run weekly. The markdown decision sheet is
the meeting agenda. TIER_1_DROP needs a drop/contest decision.
TIER_2_INVESTIGATE needs an owner assigned to the caveats.
TIER_3_PROCEED continues under standard due diligence.

**Individual target deep-dive.** Re-run with `--debug` if YAML fails.
Read `action_items`: a fired rule's `falsification_experiment` is your
next experiment; a substantive caveat is the question you need to
answer to upgrade the verdict.

**Literature claims and external candidates.** Audit before reviewers
do. A claim that survives the audit comes with its own pre-built
defense.

**Published audits.** Pin the ruleset SHA in the citation, not the
version number. Two `v1.4.0` releases could diverge; two identical
`ruleset_sha256` values cannot.

## What this tool does NOT do

The audit does not predict clinical success. A SURVIVED verdict means
the claim is internally coherent and not falsified by public-data
lookups at the indicated evidence tier. It does not mean the target
will succeed in trials.

BACE1 illustrates this. The mechanism is validated by human genetics
(APP A673T), so the framework does not falsify it. Every Phase 3
BACE1 inhibitor trial has failed. The framework reports the
`clinical_outcome_contested` caveat to flag this gap. The researcher
makes the clinical bet; the framework keeps the mechanism claim honest.

