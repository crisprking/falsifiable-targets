# Changelog



## [1.4.2] — 2026-05-26 — Documentation alignment + ft-audit input validation

**Patch release. No rule changes. Ruleset SHA unchanged from v1.4.0/v1.4.1.**

### Documentation fixes (high impact)

`docs/ADAPTER_PROTOCOL.md` and `docs/CLAIM_SCHEMA.md` field-name reference
tables now match what the engine actually reads. Pre-v1.4.2 users who
followed the docs literally would have hit silent calibration drift on:

- Oncology claims using `somatic_driver` (doc name) instead of
  `somatic_driver_evidence` (engine reads this) silently lost R3 genetics
  support → unexpected FALSIFIED_WITH_CAVEATS on EGFR, BRAF, etc.
- Clinically-failed-but-mechanistically-plausible claims (CETP, BACE1,
  SYK family) had no way to encode the failure mode in docs. The engine
  reads `clinical_outcome_contested: true` under `genetics`; now documented.
- `pubpeer_serious_concerns`, `max_phase`, `open_targets_score`, others
  documented but never read; now flagged as v1.5 placeholders.

### UX fix: `ft-audit` now surfaces validation errors

Previously raised Python tracebacks on malformed input. Now exits with
code 5 and a friendly message for empty file / bad YAML / missing
required fields / invalid claim_type.

### Verified

ft-smoke 11/11 sentinels still pass. Ruleset SHA unchanged
(`35ef2b2ab5363298...`). TYK2 audit produces byte-identical report to
v1.4.1.


## [1.4.1] - 2026-05-26 — Packaging hotfix

**Hotfix release. No rule changes.** Fixes a packaging bug in v1.4.0 where
`sentinels/v1_sentinels.yaml` and `claims/*.yaml` were not shipped with the
wheel, causing `FileNotFoundError` on clean installs (Kaggle, Colab, fresh
CI environments).

Ruleset SHA: `35ef2b2ab5363298097962a0b6ae52c70d551a1edddc341054f75cb6e4fb7221` (unchanged from v1.4.0 / ruleset v1.2.0)
Tests: 58 (unchanged)

### Fixed

- **Packaging: ship `sentinels/` and `claims/` YAML data files.** v1.4.0
  declared `sentinels/` and `claims/` as data directories but they weren't
  packages (no `__init__.py`), so setuptools silently dropped them from the
  wheel. Result: `ft-smoke` and `ft-audit` failed at runtime with
  `FileNotFoundError` on any non-editable install. Fixed by:
  - Adding empty `sentinels/__init__.py` and `claims/__init__.py`
  - Declaring both in `[tool.setuptools.packages.find].include`
  - Adding `package-data` entries for each
  - Adding a `MANIFEST.in` for the sdist path
- This unblocks `pip install falsifiable-targets` (future PyPI release),
  `pip install git+https://github.com/crisprking/falsifiable-targets.git`,
  and `pip install ./falsifiable-targets-1.4.1.tar.gz`.

### Notes

- No rule logic was touched. The ruleset SHA is byte-identical to v1.4.0.
- Editable installs (`pip install -e .`) worked in v1.4.0 because they
  read directly from the working tree; only built distributions were
  broken. This is why the bug wasn't caught earlier — local development
  loops always used `-e .`.
- Tests in `tests/test_portability.py` are extended in this release to
  catch the same bug in the future: a new test installs into a temp
  prefix and confirms the data files arrive at the expected site-packages
  locations.



All releases tag a git ref and stamp every audit produced under that release
with a ruleset SHA. The SHA-stability test in `tests/test_sentinels.py`
prevents silent drift; any rule-logic change must bump the ruleset version
and the corresponding lock value.

## [1.4.0] - 2026-05-26 — Packaging, portability, the moat release

**Infrastructure release. No rule changes.** This release turns the framework
from "audit scripts on Kaggle" into "a tool other researchers can use."

Ruleset SHA: `35ef2b2ab5363298097962a0b6ae52c70d551a1edddc341054f75cb6e4fb7221` (unchanged from v1.3.0 / ruleset v1.2.0)
Tests: 58 (was 38; +4 portability, +16 schema-validation)

### Added

- **`pyproject.toml`** — modern PEP 621 packaging. Install with
  `pip install -e ".[all]"`. Optional extras: `[validate]` (pydantic),
  `[dev]` (pytest, ruff).
- **CLI entry points**:
  - `ft-audit claim.yaml [--no-live] [--json-out FILE]`
  - `ft-smoke` — runs the sentinel calibration suite
  - `ft-validate claim.yaml` — schema check before running an audit
- **`validate_claim.py`** — pydantic-preferred schema validator with
  manual fallback. Catches typos, missing required fields, unknown
  claim types, malformed UniProt IDs, unknown fixture sections.
- **`_version.py`** — single source of truth for tool version and
  ruleset SHA. CI verifies the git tag matches.
- **`adapters/protocol.py`** — `Adapter` Protocol (runtime_checkable)
  documenting the contract for writing new data adapters. Built-in
  adapters verified conforming.
- **Worked example claims**:
  - `claims/example_hmgcr_statin.yaml` — canonical SURVIVED
  - `claims/example_synthetic_retracted.yaml` — canonical FALSIFIED via R5
  - `claims/example_novel_caveats.yaml` — canonical FALSIFIED_WITH_CAVEATS
- **`LICENSE`** — Apache License 2.0
- **`CITATION.cff`** — GitHub-rendered "Cite this repository" metadata
- **`CONTRIBUTING.md`** — full contribution guide: audits, sentinels,
  rule-change RFC process, adapter writing
- **`docs/WHY_THIS_TOOL.md`** — the moat doc: what makes the framework
  different from a checklist
- **`docs/ADAPTER_PROTOCOL.md`** — contributor guide for new adapters,
  including canonical fixture field names per section
- **`docs/CLAIM_SCHEMA.md`** — claim YAML field reference
- **`tests/test_portability.py`** — 4 tests that catch hardcoded
  machine-specific paths (`/kaggle/`, `/Users/foo/`) before they ship
- **`tests/test_validation.py`** — 16 tests pinning the schema validator
  against shipped claim files and sentinels
- **`.github/workflows/ci.yml`** — matrix CI: Python 3.10/3.11/3.12 on
  Ubuntu + Python 3.12 on macOS. Runs smoke + pytest + reproducibility
  checks on every push. Lint via ruff. Tag/version consistency check
  on release tags.

### Fixed

- **Critical portability bug**: `tests/test_sentinels.py`,
  `tests/test_adapters.py`, `tests/test_tyk2_audit.py` had hardcoded
  `/kaggle/working/falsifiable-targets` as the project root, causing
  the test suite to fail on any non-Kaggle environment. All three now
  use `Path(__file__).resolve().parent.parent`. The portability test
  guarantees this regression cannot recur.
- **Cosmetic**: `smoke_test.py` shipped banner now reads the ruleset
  version from `_version.py` (was hardcoded "v1.1.0" through v1.3.x).
- **Cosmetic**: unused `last_err` variable in `adapters/io.py`
  retry loop, flagged by ruff.

### Changed

- **Audit JSON report schema bumped to v1.1** (was v1.0): now includes
  `tool: {name, version, python_version}` and `adapter_inventory`.
  Old v1.0 reports still parse; new audits stamp v1.1.
- **`adapters/io.py`** User-Agent string updated to `falsifiable-targets/1.4.0`.

### Verified

Full local verification (mirrors what CI runs):
- `python smoke_test.py` — 11/11 sentinels pass under ruleset v1.2.0
- `python -m pytest tests/` — 58/58 tests pass (offline mode, AE_OFFLINE=1)
- `ft-validate claims/*.yaml` — all shipped claim files validate
- `ft-audit claims/example_hmgcr_statin.yaml --no-live` — SURVIVED
- `ft-audit claims/example_synthetic_retracted.yaml --no-live` — FALSIFIED
- `ft-audit claims/example_novel_caveats.yaml --no-live` — FALSIFIED_WITH_CAVEATS
- `ft-audit claims/ipi1_madurella.yaml --no-live` — FALSIFIED_WITH_CAVEATS (2 substantive caveats)
- `ft-audit claims/tyk2_psoriasis.yaml --no-live` — SURVIVED
- `ruff check .` — clean

### What this release does NOT do

- No rule logic changes. Ruleset SHA stays locked at v1.2.0's value.
- No new audits. The Ipi1 and TYK2 audits stand under v1.2.0 ruleset.
- No live adapter behavior changes (the User-Agent string change is
  cosmetic; cache keys are URL-based, not header-based).

The next release (v1.5.0) is queued for substantive content changes:
either R3 live Open Targets adapter, or v1.4 compound-level overlap
work for R6 (the gap surfaced by the v1.3.0 TYK2 re-audit).

---

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
