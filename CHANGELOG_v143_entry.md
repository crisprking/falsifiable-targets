## [1.4.3] — 2026-05-27 — Exit-code contract hardening + portfolio audit + hard-test suite

**Patch release. No rule changes. Ruleset SHA unchanged.**

### `run_audit.py` hardening

- Wrap `main()` invocation in a `_run()` entry point with explicit
  exception-to-exit-code translation. Unhandled exceptions now exit
  5 (per documented contract) instead of leaking Python default
  exit code 1.
- Add `--debug` flag for full tracebacks on error paths. Default
  behavior unchanged: one-line `ERROR:` message to stderr.
- `KeyboardInterrupt` now produces exit 130 with a clean message.

### New: production portfolio audit (`scripts/portfolio_audit.py`)

Batch-audits a directory of claim YAMLs and emits three output
formats: CSV (for spreadsheets), Markdown (for PR descriptions /
wiki), and JSON manifest (for downstream automation).

The `action_items` column surfaces signal in priority order:

1. **Fired rules**: rule ID + falsification tier + recommended experiment
2. **Substantive caveats**: rule ID + caveat text
3. **Data gaps**: abstained rules (the "need more input" signal)
4. **All clear**: passes every applicable rule

Risk-tier classification (`TIER_1_DROP` / `TIER_2_INVESTIGATE` /
`TIER_3_PROCEED` / `TIER_0_FIX_INPUT`) gives leadership a single
sortable column for portfolio triage.

Portfolio exit code reflects worst tier across claims:
2 (drop present), 5 (broken input), 1 (investigate present), 0 (clean).

### New: six-contract acceptance gate (`scripts/acceptance_gate.py`)

Verifies before release:

1. Verdict + exit-code mapping across 5 reference scenarios
2. Error paths exit 5 with no tracebacks
3. Ruleset SHA pinned across all reports
4. Determinism: back-to-back runs byte-identical modulo timestamp
5. `--debug` accepted, preserves exit codes, adds tracebacks
6. Substantive caveat text propagated into report JSON

### New: hard-test scenario suite (`tests/hardtest_claims/`)

12 verified scenario YAMLs covering every verdict tier:

- **s1 PCSK9** — canonical SURVIVED (validated_mechanism)
- **s2 BACE1** — FALSIFIED_WITH_CAVEATS (mechanism real, clinical contested)
- **s3 OCT4/STAP** — FALSIFIED (retracted paper)
- **s4 SAT/HDAC4** — FALSIFIED (R6 paralog-class collapse)
- **s5 fungal kinase** — FALSIFIED_WITH_CAVEATS (orthology + selectivity gaps)
- **s6a-d** — validator stress: bad enum, missing field, empty file, parse error
- **s7 KRAS-G12C** — FALSIFIED_WITH_CAVEATS (partial-success target)
- **s8 CCR5** — SURVIVED (extraordinary claim with extraordinary evidence)
- **s9 ZIKV NS5** — FALSIFIED_WITH_CAVEATS (ambiguous evidence)

### Docs

- `docs/PORTFOLIO_USAGE.md` — workflow guide with risk-tier semantics,
  exit-code contract, and action_items reference.

### Verified

All six contracts pass on Python 3.10+, Linux. CCR5/HIV correctly
SURVIVES despite `extraordinary_claim` classification because the
evidence matches the bar. BACE1 correctly demoted to FALSIFIED_WITH_CAVEATS
without any rule firing, driven entirely by substantive caveat from R3.

---

