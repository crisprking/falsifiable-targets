# Changelog

All releases tag a git ref and stamp every audit produced under that release
with a ruleset SHA. The SHA-stability test in `tests/test_sentinels.py`
prevents silent drift; any rule-logic change must bump the ruleset version
and the corresponding lock value.

## [1.3.1] - 2026-05-26 — v1.3.0 audit narrative + honest finding

**Documentation release. No code changes. No rule changes.**

Ruleset SHA: `35ef2b2ab5363298097962a0b6ae52c70d551a1edddc341054f75cb6e4fb7221` (unchanged from v1.3.0 / ruleset v1.2.0)
Tests: 38 (unchanged from v1.3.0)

Added:

- `docs/AUDIT_TYK2_v1_3.md` — narrative for the TYK2 re-audit under expanded R6
- `docs/AUDIT_LIMITATIONS_v1_3.md` — updated limitations doc; R6 scope gap from v1.2.x is closed; new gap (compound-level overlap vs pool-size overshadow) is named

**Calibration finding from the v1.3.0 TYK2 re-audit:**

A prediction was published in the v1.3.0 release notes: TYK2 would
likely shift from `SURVIVED` to `FALSIFIED_WITH_CAVEATS` because the
JAK family chemistry pools were assumed to dwarf TYK2's by ≥ 2×.
**The prediction was wrong.** Live ChEMBL data showed:

| Kinase | Distinct compounds |
|---|---|
| TYK2 (primary) | 538 |
| JAK1 (paralog) | 549 |
| JAK2 (paralog) | 772 |
| JAK3 (paralog) | 832 |

Maximum paralog/primary ratio: **1.55×**. Below the heuristic
threshold of 2.0×. R6 returns PASSED without caveat.

The honest interpretation, written in [`docs/AUDIT_TYK2_v1_3.md`](docs/AUDIT_TYK2_v1_3.md):
TYK2 cleared the heuristic at the pool-size axis. R6 v1.2.0 measures
pool-size overshadow, not compound-level overlap. The pan-JAK selectivity
concern that deucravacitinib's JH2-domain mechanism resolves operates at
the latter axis, which is a v1.4+ feature. The v1.3.0 SURVIVED is
qualitatively stronger than the v1.1.0 SURVIVED — R6 was actually
evaluated against live paralog data instead of returning NOT_APPLICABLE.

No threshold adjustment was made after seeing the result. Lowering the
threshold to "force" R6 to fire would be motivated reasoning of exactly
the kind the framework is supposed to prevent.

## [1.3.0] - 2026-05-26 — R6 scope expansion + paralog-ratio heuristic

**Ruleset change. New SHA lock.**

Ruleset SHA: `35ef2b2ab5363298097962a0b6ae52c70d551a1edddc341054f75cb6e4fb7221` (new, v1.2.0 ruleset)
Previous SHA: `2f9aab7d0ebc209f62c16eb35be31bc5b65fa2eb09adc02bea5bff5176269b32` (v1.1.0)
Tests: 38

Changed:

- **R6_chemistry_class_collapse** version bumped from 1.0.0 to 1.2.0. Two changes: (1) scope expanded to `validated_mechanism` claims (previously only `novel_target` and `chemistry_series`), and (2) new heuristic substantive-caveat path: when `chembl_paralog_compound_counts` is populated and any paralog has >= 2x the primary's compound count, R6 emits a substantive caveat naming the overshadow. The existing exact-fraction falsification path (>= 0.80) is preserved. Rationale: the pan-class selectivity question (canonical example: pan-JAK overlap for TYK2) is exactly the structural risk R6 was written to catch.

Added:

- `claims/paralog_map.yaml` — manually curated paralog map keyed by primary UniProt accession. Initial entries: P29597 (TYK2) → JAK1/JAK2/JAK3. PCSK9 and TNF intentionally omitted with rationale notes.
- `ChEMBLAdapter.__init__` now accepts a `paralog_map` parameter. When supplied and the claim's UniProt accession has an entry, the adapter fetches compound counts for each paralog and attaches them to the chemistry section as `chembl_paralog_compound_counts`.
- `adapters.load_paralog_map(path)` — convenience loader for paralog_map.yaml.
- `run_audit.py` automatically loads `claims/paralog_map.yaml` (if present) and passes it to the live composite adapter.
- Sentinel `KINASE_CLASS_COLLAPSE_VALIDATED` exercising the new heuristic path directly. Predicts `FALSIFIED_WITH_CAVEATS` for a synthetic validated-mechanism claim with paralog pool 5x primary.
- Four new tests in `tests/test_tyk2_audit.py` and four new paralog-map tests in `tests/test_adapters.py`.

Compatibility:

- Existing v1.1.0 sentinels all still pass with unchanged verdicts. The R6 changes are strictly additive: the heuristic path only fires when `chembl_paralog_compound_counts` is present in the chemistry section, which v1.1.0 fixtures never supplied.
- The TYK2 fixture-only test under v1.2.0 still returns SURVIVED (no paralog counts in fixture). The live-mode TYK2 audit on Kaggle is expected to now fire the R6 substantive caveat — that's the v1.3.0 calibration finding.

Impact on existing audits (verdict shape comparison):

- **Ipi1** under v1.2.0 ruleset: unchanged. `FALSIFIED_WITH_CAVEATS`, 2 substantive caveats (R1 + R7). No paralog map entry for Ipi1.
- **TYK2** under v1.2.0 ruleset, live mode: expected to shift from `SURVIVED` (v1.1.0) to `FALSIFIED_WITH_CAVEATS` because the JAK family paralog chemistry pool dwarfs TYK2's. This is a real finding: the framework now catches the pan-JAK selectivity concern that v1.2.x missed. The deucravacitinib defense (JH2-domain binding) is a *structural narrative* not encoded in any rule input; the v1.2.0 caveat reads honestly that the heuristic flags the risk and structural support is required.

## [1.2.1] - 2026-05-25 — Public artifact

**Documentation release. No code changes. No rule changes.**

Ruleset SHA: `2f9aab7d0ebc209f62c16eb35be31bc5b65fa2eb09adc02bea5bff5176269b32` (unchanged from v1.1.0)
Tests: 31 (unchanged from v1.2.0)

Added:

- `README.md` — public entry point with headline result and reproduce-locally instructions
- `docs/AUDIT_LIMITATIONS_v1_2.md` — honest per-rule accounting of live-data vs fixture-data sources, R6 scope gap, v1.3.0 roadmap
- `docs/AUDIT_IPI1_v1_1.md` — human-readable narrative of the inaugural self-audit
- `docs/AUDIT_TYK2_v1_2.md` — human-readable narrative of the first external audit
- `docs/MADURELLA_M7_HOOK.md` — short doc that drops into the Madurella repo as M7, cross-references this repo by SHA
- `CHANGELOG.md` — this file
- `.gitignore` — Python / pytest / cache ignores

## [1.2.0] - 2026-05-25 — First external audit

Ruleset SHA: `2f9aab7d0e...` (unchanged from v1.1.0)
Tests: 31 (10 sentinels + 3 v1.1 regression pins + 15 adapter tests + 2 TYK2 pins + 1 SHA stability)

Added:

- `claims/tyk2_psoriasis.yaml` — the first external audit subject. Methodology notes on every hand-set field.
- `run_audit.py` — CLI runner. Builds the live+fixture composite adapter, runs the audit, emits console + JSON report with claim SHA and ruleset SHA stamps.
- `tests/test_tyk2_audit.py` — regression pins for TYK2 SURVIVED + R6 NOT_APPLICABLE under v1.1.0.

Fixed:

- `ChEMBLAdapter` now distinguishes "live fetch failed / offline" (returns `{}` so the composite falls back to fixture) from "API responded with zero compounds" (returns `{"chembl_distinct_compounds": 0}`). Caught by an offline-mode dry run; would have caused spurious R2 falsification if not fixed.

Audit run: TYK2 / psoriasis → `SURVIVED`, 0 substantive caveats, 0 falsifications. Live UniProt returned `pdb_count=62`. Live ChEMBL returned `chembl_distinct_compounds=538`.

## [1.1.0] - 2026-05-25 — R7 substantive-caveat upgrade

Ruleset SHA: `2f9aab7d0ebc209f62c16eb35be31bc5b65fa2eb09adc02bea5bff5176269b32` (new)
Previous SHA: `6470f7523a373d018cbc367e13d6d26c3e585613f8b29fe27ca3c57763abe787`
Tests: 28

Changed:

- **R7_selectivity_counterscreen** version bumped from 1.0.0 to 1.1.0. For `novel_target` claims with no selectivity data, the rule now emits a **substantive** caveat (was operational). Rationale: for a novel target, the selectivity-vs-paralog question is the entire reason the claim is provisional, so its absence is the headline gap, not a tooling gap. Other claim types unchanged.

Added:

- Sentinel `NOVEL_TARGET_SELECTIVITY_GAP` exercising the new R7 behavior directly.
- `test_ipi1_inaugural_audit_two_substantive_caveats` regression test pinning Ipi1's caveat shape at `[R1_orthology, R7_selectivity_counterscreen]`.
- `test_novel_target_gap_has_r7_substantive_caveat` regression test for the new sentinel.

Impact on existing audit: Ipi1 verdict unchanged (`FALSIFIED_WITH_CAVEATS`, score 0.0), substantive caveats went 1 → 2. The framework now surfaces both axes of provisionality the Madurella M6 record flagged in prose.

## [1.0.1] - 2026-05-25 — Live adapters

Ruleset SHA: `6470f7523a...` (unchanged from v1.0.0)
Tests: 25

Added:

- `adapters/io.py` — `FixtureAdapter`, `CompositeAdapter`, `UniProtAdapter`, `ChEMBLAdapter`, `default_composite`. UniProt and ChEMBL adapters support both mock mode (hermetic, used in tests) and live mode (HTTP with on-disk cache under `.ae_cache/`). `AE_OFFLINE=1` env var forces cache-only.
- `tests/test_adapters.py` — 15 hermetic tests against mock UniProt and ChEMBL payloads.

Verified live against PCSK9 (Q8NBP7 / CHEMBL2929): UniProt returned 62 PDB refs, real tissue text; ChEMBL returned 495 distinct compounds. End-to-end live + fixture composite produces correct merged data for PCSK9.

This release intentionally changes no rule logic. The SHA stability test continues to lock the v1.0.0 ruleset.

## [1.0.0] - 2026-05-25 — Initial ship

Ruleset SHA: `6470f7523a373d018cbc367e13d6d26c3e585613f8b29fe27ca3c57763abe787`
Tests: 10 (9 sentinels + 1 SHA stability)

Initial release:

- Seven rules: R1 orthology, R2 chemistry support, R3 genetics support, R4 expression, R5 replication, R6 chemistry class collapse, R7 selectivity counter-screen
- Nine sentinel calibration cases covering positive, negative, and edge categories (PCSK9, BCR-ABL, TNF; SAT/HDAC4 phantom, STAP-cells, arsenic-life; cruzain peptidomimetics, CETP/LDL, TEX10/Ipi1 orthology threshold)
- Inaugural audit: Ipi1 (KXX81897.1) → `FALSIFIED_WITH_CAVEATS`, 1 substantive caveat under v1.0.0 (later 2 under v1.1.0)
- Stdlib-only smoke test (no third-party deps beyond PyYAML)
- Verdict aggregator with cheapest-falsification ranking across the five falsification tiers

The framework's first audit was its creator's own published headline target. This was a deliberate choice about the order of legitimacy.
