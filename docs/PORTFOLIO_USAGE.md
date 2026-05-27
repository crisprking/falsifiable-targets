# Portfolio audit usage guide

## What this is

A workflow for batch-auditing drug-discovery target claims against the
falsifiable-targets framework. Turns a directory of claim YAMLs into a
triaged decision sheet with verdict, risk tier, fired rules, substantive
caveats, and next-experiment recommendations per claim.

## Output contract

Each run produces three files under `reports/`:

- `portfolio_decision.csv` for spreadsheets and program managers.
- `portfolio_decision.md` for PR descriptions and wiki pages, sorted by
  risk tier with highest-risk claims first.
- `portfolio_manifest.json` for downstream automation, including
  per-claim rows, tier and verdict summaries, and the run timestamp.

## Risk tier semantics

| Tier | Source verdict | What to do |
|---|---|---|
| `TIER_1_DROP` | FALSIFIED | A rule fired with concrete public-data evidence. Cut from portfolio unless the falsification is contested with new evidence. |
| `TIER_2_INVESTIGATE` | FALSIFIED_WITH_CAVEATS or INSUFFICIENT_DATA | Substantive caveats present; mechanism may be valid but specific concerns must be addressed before committing budget. |
| `TIER_3_PROCEED` | SURVIVED | Passes all applicable rules. Standard due diligence still applies. |
| `TIER_0_FIX_INPUT` | ERROR | YAML failed schema validation. Fix the input file and re-run. |

## Reading the action_items column

Priority order for what the column surfaces:

1. **Fired rules**: `R5_replication FALSIFIED @ public_data_lookup: <recommended experiment>`. Which rule failed, at what evidence tier, and what experiment would resolve the question.
2. **Substantive caveats**: `R3_genetics_support: clinical outcome contested; ...`. The rule did not fail outright, but raised a concern significant enough to demote the verdict. Read the text; it explains the concern.
3. **Data gaps**: `data gap: R2_chemistry_support abstained (no input)`. The rule had no data to evaluate. This is a gap, not a finding. Fill the gap (ChEMBL pull, selectivity counter-screen, etc.) and re-audit.
4. **All clear**: `all rules cleared`. The claim passes every applicable rule.

## Exit-code contract

```
0 - SURVIVED (or portfolio: all TIER_3_PROCEED)
1 - FALSIFIED_WITH_CAVEATS (or portfolio: TIER_2_INVESTIGATE present)
2 - FALSIFIED (or portfolio: TIER_1_DROP present)
3 - INSUFFICIENT_DATA
5 - ERROR (missing/empty/malformed claim file or unhandled exception)
```

For CI integration, treat exit 5 differently from exit 1 or 2: 1 and 2
are scientific findings (the audit ran successfully and the verdict is
what it is), while 5 indicates a broken input file that needs to be
fixed before the audit is meaningful.

## CLI usage

```
# Single claim
python run_audit.py claims/my_target.yaml --no-live --json-out reports/my_target.json

# Show full traceback on error paths (default: one-line stderr)
python run_audit.py claims/my_target.yaml --debug

# Portfolio batch
python scripts/portfolio_audit.py claims/ reports/

# Acceptance gate (six contracts)
python scripts/acceptance_gate.py
```

## Required claim YAML structure

```yaml
claim:
  target_symbol: "PCSK9"               # required
  uniprot_id: "Q8NBP7"                 # required
  indication: "Hypercholesterolemia"   # required
  mechanism: |                         # required, multiline OK
    Mechanism description goes here.
  claim_type: validated_mechanism      # required; one of:
                                       #   validated_mechanism,
                                       #   novel_target,
                                       #   extraordinary_claim,
                                       #   chemistry_series,
                                       #   repurposing
  citation_pmid: "16554528"            # optional but recommended
fixture:                               # optional offline test data
  orthology: {sources_agreeing: 4, sources_total: 4}
  genetics:  {gwas_hits: 31, mendelian_evidence: true}
  # ... other rule-specific fields
```

## Verified contracts (acceptance gate)

The tool passes six contracts, verified by `scripts/acceptance_gate.py`:

1. **Verdict + exit-code mapping** across five reference scenarios
   (PCSK9, BACE1, STAP, paralog collapse, weak fungal kinase).
2. **Error paths exit 5** with actionable stderr messages and no Python
   tracebacks bleeding into output.
3. **Ruleset SHA pinned**: all reports reference the same
   content-addressed ruleset.
4. **Determinism**: two back-to-back runs on the same claim produce
   byte-identical reports (modulo audit timestamp).
5. **`--debug` flag** accepted by argparse, preserves exit codes, adds
   tracebacks on exceptions that escape `main()`.
6. **Substantive caveat text** propagated from rule-internal caveats
   into portfolio output, with non-empty text per caveat.

Three additional complex scenarios verified manually (s7-s9):

- **KRAS-G12C / NSCLC**: FALSIFIED_WITH_CAVEATS (mechanism real;
  clinical outcome contested - the partial-success target case).
- **CCR5 / HIV cure**: SURVIVED (extraordinary_claim with extraordinary
  evidence - Mendelian protection plus two confirmed cures).
- **ZIKV NS5**: FALSIFIED_WITH_CAVEATS (R1 orthology disagreement flagged
  - the ambiguous-evidence case).

## Workflow recommendations

**Program-level triage.** Run weekly. The markdown decision sheet is
the meeting agenda. TIER_1_DROP claims need a decision (drop, or
contest the falsification with new evidence). TIER_2_INVESTIGATE
claims need an owner assigned to address the listed caveats.
TIER_3_PROCEED claims continue under standard due diligence.

**Individual target deep-dive.** If your YAML fails validation, re-run
with `--debug` to get a traceback when the error is an escaped
exception. Read the `action_items` field: when a rule fires, the
`falsification_experiment` text is your next experiment. When a
substantive caveat appears, the caveat text is the question you need
to answer to upgrade the verdict.

**Literature claims and external candidates.** Audit before reviewers
do. The framework is designed to surface the things a skeptical
reviewer would surface, but at writing time rather than rebuttal time.

**Published audits.** Pin the ruleset SHA in the citation, not the
version number. Two `v1.4.0` releases could diverge; two identical
`ruleset_sha256` values cannot.

## What this tool does NOT do

The audit does not predict clinical success. A `SURVIVED` verdict means
the claim is internally coherent and not falsified by public-data
lookups at the indicated evidence tier. It does not mean the target
will succeed in trials.

BACE1 illustrates this. The mechanism is validated by human genetics
(APP A673T), so the framework does not falsify it. Every Phase 3 BACE1
inhibitor trial has failed. The framework reports the
`clinical_outcome_contested` caveat to flag this gap. The researcher
makes the clinical bet; the framework keeps the mechanism claim honest.
