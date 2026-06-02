# Contributing to falsifiable-targets

Thanks for your interest. This document is the operating manual for
contributing audits, sentinels, rule changes, or adapters.

The framework's value proposition is **deterministic, SHA-locked
verdicts**. That cuts both ways: useful contributions move the
framework forward; sloppy contributions corrupt the ruleset hash and
break every prior audit's reproducibility. The processes below exist
to make the former easy and the latter impossible.

---

## What you can contribute

| Type | Difficulty | Process |
|---|---|---|
| **A new audit** (claim YAML + JSON report) | Easy | PR with claim, JSON, narrative doc |
| **A new sentinel** (calibration test case) | Medium | PR with sentinel YAML entry + justification |
| **A new adapter** (Open Targets, AlphaFold, PRIDE...) | Medium | PR with `adapters/<name>.py`, mock tests |
| **A rule fix / scope change** | Hard | RFC issue first, then PR; ruleset version bumps |
| **A new rule** | Hard | RFC issue first, then PR; requires sentinel proof |
| **Documentation, typos, examples** | Easy | Direct PR |
| **Tooling, CI, packaging** | Easy-Medium | Direct PR for small; issue first for large |

---

## Process: contributing a new audit

This is the most common contribution path. You believe a published
target claim is either well-supported or fragile; you want the
framework to weigh in.

1. **Fork & clone**, then create a branch: `audit/<symbol>-<indication>`
2. **Write the claim YAML** at `claims/<symbol>_<indication>.yaml`.
   See `claims/tyk2_psoriasis.yaml` for the canonical template.
3. **Validate it**:
   ```bash
   ft-validate claims/your_claim.yaml
   ```
4. **Run the audit** (live mode if you have internet, `--no-live`
   otherwise):
   ```bash
   ft-audit claims/your_claim.yaml --json-out reports/your_audit.json
   ```
5. **Write the narrative** at `docs/AUDIT_<SYMBOL>_v<ruleset>.md`.
   Required sections: claim summary, fixture provenance (cite every
   hand-set value), verdict, what the framework caught, what it
   *didn't* and why, "what would change the verdict."
6. **Open the PR**. The CI runs the existing test suite; we don't
   require the live audit to pass in CI (no network).

Honesty discipline: if your audit produces a verdict you disagree with,
**ship the verdict, then queue a rule change in a separate PR**.
Audits that "look right because we tweaked R6 in the same PR" are how
ruleset SHAs drift and how the framework loses credibility.

---

## Process: contributing a new sentinel

Sentinels are calibration anchors. Each one says: "This claim, in this
exact fixture state, *must* produce this verdict, or the calibration
is wrong."

Adding a sentinel is an assertion about how the framework *should*
behave. Process:

1. Identify a real-world calibration gap. Examples:
   - A published falsification mode the framework doesn't catch
   - A class of legitimate claim the framework currently abstains on
   - A specific failure mode (e.g. phantom-cofactor like SAT/HDAC4)
2. Encode it in `sentinels/v1_sentinels.yaml` with a unique `id`,
   complete `claim`, full `fixture`, and `expected_verdict`.
3. Run `python smoke_test.py`. Two outcomes:
   - **Passes**: the framework already handles this case. Your
     sentinel still has value as a regression pin — submit the PR.
   - **Fails**: the framework doesn't yet handle this case. This is
     the more interesting contribution. Submit the sentinel, and in
     the PR description note this is a *calibration gap*. We'll
     either accept the sentinel as a known-failing case (with a
     future-rule-fix milestone) or accept your accompanying rule
     change (see next section).

---

## Process: rule changes

This is the most consequential type of contribution because rule
changes bump the **ruleset SHA**, which invalidates the byte-level
reproducibility of every prior audit run under the old ruleset.

Required process:

1. **Open an RFC issue first**. Title: `RFC: <one-line description>`.
   Body must contain:
   - The problem (which sentinel or audit surfaced the gap)
   - The proposed rule change (exact text of `applies_to` and
     `evaluate` semantics)
   - Backward-compat analysis: which existing sentinels' verdicts
     change? Which existing audits' verdicts change?
   - Version-bump proposal (`patch`/`minor`/`major` on the ruleset)
2. **Discussion in the issue.** Don't write code until the RFC is
   accepted. RFCs that change a rule's *verdict* on more than one
   existing audit get extra scrutiny.
3. **Implementation PR** referencing the RFC issue. Required artifacts:
   - The rule change in `smoke_test.py`
   - A version bump on the rule (e.g., `version: "1.0.0"` → `"1.1.0"`)
   - A version bump on the ruleset descriptor in `_version.py`
   - The new ruleset SHA recomputed and pinned in `tests/test_sentinels.py`
   - At least one new or modified sentinel exercising the new
     behavior
   - Re-running any audits whose verdicts the change affects, with
     updated narrative docs
4. **CHANGELOG entry.** Required for every rule change.

Anti-pattern we explicitly reject: *quietly* changing a rule because
an audit produced an inconvenient verdict. If TYK2 v1.3.0 had
returned `FALSIFIED_WITH_CAVEATS` and we tuned the R6 threshold from
2.0 to 3.0 to make it pass, the framework would be useless. The
threshold is locked at 2.0 by sentinel KINASE_CLASS_COLLAPSE_VALIDATED;
moving it requires moving the sentinel, which requires justifying
why the calibration shifted, which requires an RFC. That's the moat.

---

## Process: contributing an adapter

Adapters bring more data into the audit. New adapters DO NOT change
the ruleset SHA (no rule logic changes), but they can change verdicts
by populating previously-empty fixture sections.

1. Implement in `adapters/<name>.py`. Conform to the `Adapter`
   protocol in `adapters/protocol.py` — just a `get(section, claim)`
   method returning a dict.
2. Write hermetic mock tests in `tests/test_adapters.py` (or a new
   `tests/test_<name>_adapter.py`). Live API tests go in
   `tests/live/` and are skipped by CI; use `pytest.mark.live`.
3. Cache to `.ae_cache/` so repeated runs are deterministic.
4. Respect `AE_OFFLINE=1` env var.
5. Update `docs/ADAPTER_PROTOCOL.md` if your adapter exercises a new
   convention.

---

## Local development

```bash
git clone https://github.com/crisprking/falsifiable-targets.git
cd falsifiable-targets
python -m venv .venv && source .venv/bin/activate
pip install -e ".[all]"

# The basic checks every PR must pass:
python smoke_test.py              # sentinel calibration
python -m pytest tests/           # test suite
ft-validate claims/*.yaml         # all claims pass schema

# Style:
ruff check .
ruff format --check .
```

---

## Code style

- Python 3.10+ syntax. Type hints encouraged but not required.
- `ruff` for linting and formatting (config in `pyproject.toml`).
- Comments explain *why*, not *what*. The code already shows what.
- Tests must be hermetic (no network) unless explicitly marked `live`.
- Imports: stdlib, then third-party, then local. Within each group,
  alphabetical.

---

## Releases

The maintainer cuts releases. The process:

1. Update `_version.py` with the new tool version and (if changed)
   the new ruleset version + SHA.
2. Update `CHANGELOG.md` with the changes.
3. Update `CITATION.cff` `version` and `date-released`.
4. Commit, tag `v<version>`, push tag.
5. GitHub Actions auto-publishes to PyPI (when configured).

---

## Code of conduct

Be kind. Engage with the technical substance. Disagreement about
whether R6's threshold is well-calibrated is welcome; ad hominem isn't.

---

## Questions

Open an issue. The maintainer is `crisprking` on GitHub.
