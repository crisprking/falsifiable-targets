# Changelog

All releases tag a git ref and stamp every audit produced under that release
with a ruleset SHA. The SHA-stability test in `tests/test_sentinels.py`
prevents silent drift; any rule-logic change must bump the ruleset version
and the corresponding lock value.

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
