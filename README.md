# falsifiable-targets

[![CI](https://github.com/crisprking/falsifiable-targets/actions/workflows/ci.yml/badge.svg)](https://github.com/crisprking/falsifiable-targets/actions/workflows/ci.yml)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

A **deterministic, SHA-locked target-validation audit framework** for
drug-discovery claims. Encodes seven rules that consume a structured
claim YAML and emit a verdict: `SURVIVED`, `FALSIFIED_WITH_CAVEATS`,
`FALSIFIED`, or `INSUFFICIENT_DATA` — with the cheapest available
falsification experiment named for each flagged rule.

Every audit is content-addressed: the ruleset is hashed (SHA-256), the
claim is hashed, the report stamps both. Reproducing an audit
byte-identically is one `git clone` and one CLI call.

> **Why this exists:** too many published target claims are defended
> by evidence that, on inspection, points at a different target — the
> SAT/HDAC4 archetype, the STAP-cells failure mode, the CETP
> clinical-genetics gap. The seven rules detect those failure modes
> from public data alone, at the cheapest falsification tier available.
> See [`docs/WHY_THIS_TOOL.md`](docs/WHY_THIS_TOOL.md) for the full
> argument.

---

## Quick start

```bash
# Install from a clone:
git clone https://github.com/crisprking/falsifiable-targets.git
cd falsifiable-targets
pip install -e ".[all]"

# Validate your claim YAML against the schema:
ft-validate claims/example_hmgcr_statin.yaml

# Run an audit (offline / fixture-only):
ft-audit claims/example_hmgcr_statin.yaml --no-live

# Run an audit with live UniProt + ChEMBL lookups + JSON report:
ft-audit claims/tyk2_psoriasis.yaml --json-out my_audit.json
```

The CLI exits with a code reflecting the verdict (0=SURVIVED,
1=FALSIFIED_WITH_CAVEATS, 2=FALSIFIED, 3=INSUFFICIENT_DATA), so it
chains into CI / data pipelines.

---

## Headline results

Three external audits in the public record, two ruleset versions, one
prediction made and falsified:

| Claim | Type | Ruleset | Verdict | Substantive caveats |
|---|---|---|---|---|
| Ipi1 / Madurella (in-house novel antifungal) | `novel_target` | v1.1.0 | `FALSIFIED_WITH_CAVEATS` | 2 |
| TYK2 / psoriasis (deucravacitinib, FDA 2022) | `validated_mechanism` | v1.1.0 | `SURVIVED` | 0 |
| TYK2 / psoriasis (re-audit, R6 expanded) | `validated_mechanism` | v1.2.0 | `SURVIVED` | 0 |

The Ipi1 vs TYK2 asymmetry is the **calibration finding**: same engine,
same SHA-locked rule logic, opposite verdicts in the direction each
claim's structure predicts. Most validation frameworks built by a
research group never produce an inconvenient verdict about that group's
own work; the framework's first action was to demote its creator's
headline target.

---

## What's in v1.4.0

v1.4.0 is the **packaging release** that turns the framework into a
tool other researchers can use:

- **Installable**: `pip install -e .` works on Mac, Linux, Windows
- **Three CLI entry points**: `ft-audit`, `ft-smoke`, `ft-validate`
- **Schema validation**: `ft-validate claim.yaml` catches typos before
  a 30-second live audit. Pydantic-preferred with manual fallback.
- **Adapter protocol** (`adapters/protocol.py`): write your own data
  source in 20 lines of Python, drop it into `CompositeAdapter`
- **Portability tests**: any hardcoded machine-specific path
  (`/kaggle/`, `/Users/foo/`) fails the CI test
- **Worked examples**: three additional claim files showing each
  non-`SURVIVED` verdict pattern
- **Apache 2.0 license**, `CITATION.cff`, `CONTRIBUTING.md`
- **GitHub Actions CI** across Python 3.10/3.11/3.12 + macOS

The ruleset is **unchanged from v1.3.0** (SHA still
`35ef2b2a...`) — this release is infrastructure only, no rule logic
touched.

---

## Reading an audit report

```json
{
  "schema_version": "1.1",
  "tool": {"name": "falsifiable-targets", "version": "1.4.0", "python_version": "3.12.3"},
  "ruleset_version": "1.2.0",
  "ruleset_sha256": "35ef2b2ab5363298097962a0b6ae52c70d551a1edddc341054f75cb6e4fb7221",
  "claim": {
    "target_symbol": "TYK2",
    "uniprot_id": "P29597",
    "claim_sha256": "..."
  },
  "verdict": "SURVIVED",
  "score": 0.0,
  "cheapest_falsification": null,
  "per_rule": [{"rule_id": "R1_orthology", "status": "not_applicable", "...": "..."}],
  "substantive_caveat_count": 0,
  "operational_note_count": 0,
  "adapter_inventory": [{"name": "FixtureAdapter", "live": false}]
}
```

Re-running the audit on any machine with the same ruleset SHA and the
same claim SHA must produce a byte-identical verdict. If it doesn't,
the adapter layer changed but the rules did not — see
[`docs/AUDIT_LIMITATIONS_v1_3.md`](docs/AUDIT_LIMITATIONS_v1_3.md) §3.

---

## The seven rules

| Rule | Audits | Cheapest tier |
|---|---|---|
| R1_orthology | Pathogen/novel targets: do ≥ majority of ortholog DBs agree on the protein's presence/absence? | public data lookup |
| R2_chemistry_support | Chemistry-series/validated: are there meaningful ChEMBL compounds? | public data lookup |
| R3_genetics_support | Validated-mechanism: GWAS, Mendelian, or somatic-driver evidence? | public data lookup |
| R4_expression | Target detectably expressed in indication-relevant tissue? | public data lookup |
| R5_replication | Not retracted; no overwhelming rebuttals without replication | public data lookup |
| R6_chemistry_class_collapse | Do ChEMBL hits collapse onto a paralog class (phantom evidence)? | public data lookup |
| R7_selectivity_counterscreen | Does selectivity data exist; for novel targets, absence is substantive | cheap in silico |

The aggregator always names the *cheapest* available experiment to
falsify a flagged claim, because experiments cheap enough to run from
a laptop are where the framework's value compounds.

---

## Limitations

Read [`docs/AUDIT_LIMITATIONS_v1_3.md`](docs/AUDIT_LIMITATIONS_v1_3.md)
before drawing strong conclusions from any audit. In particular:

- R6 audits **pool-size overshadow** (ChEMBL compound counts across
  paralogs). It does NOT audit **compound-level overlap** (whether
  individual binders are pan-family). The deucravacitinib JH2-domain
  selectivity argument lives at the latter axis, which is queued for
  v1.5.
- R3 currently relies on hand-encoded fixture data. A live Open Targets
  adapter is the v1.5 milestone (see [`CONTRIBUTING.md`](CONTRIBUTING.md))
- R5's retraction signal requires the claim file to declare retraction
  status. A live Retraction Watch adapter is queued.

The framework's gaps are documented openly. None of these limitations
are surprises to the authors.

---

## Contributing

PRs welcome. See [`CONTRIBUTING.md`](CONTRIBUTING.md) for:

- How to add a new audit (a claim YAML + a narrative doc)
- How to add a new sentinel (calibration test case)
- How to write a new adapter (live data source)
- How to propose a rule change (RFC process)
- The **honesty discipline**: rule changes that "happen to" produce a
  more convenient verdict are rejected by design.

---

## Project structure

```
falsifiable-targets/
├── smoke_test.py            # The 7-rule engine. Single file, self-contained.
├── run_audit.py             # CLI entry point: ft-audit
├── validate_claim.py        # Schema validator: ft-validate
├── _version.py              # Single source of truth for versions
├── pyproject.toml           # Modern packaging (PEP 621)
├── adapters/                # Data-source adapters
│   ├── io.py                # FixtureAdapter, UniProtAdapter, ChEMBLAdapter, Composite
│   └── protocol.py          # The Adapter Protocol (for contributors)
├── claims/                  # Claim YAMLs (machine-readable input)
│   ├── ipi1_madurella.yaml
│   ├── tyk2_psoriasis.yaml
│   ├── example_*.yaml       # Worked examples for each verdict pattern
│   └── paralog_map.yaml     # R6 v1.2.0 paralog set per primary target
├── sentinels/               # The calibration suite
│   └── v1_sentinels.yaml
├── tests/                   # 58 hermetic tests (run with pytest)
├── reports/                 # Audit JSON artifacts (SHA-stamped)
├── docs/
│   ├── WHY_THIS_TOOL.md         # The moat doc
│   ├── ADAPTER_PROTOCOL.md      # How to write an adapter
│   ├── CLAIM_SCHEMA.md          # Claim YAML field reference
│   └── AUDIT_*.md               # Narrative audit reports
└── .github/workflows/ci.yml     # CI: lint + tests + reproducibility checks
```

---

## Release history

See [`CHANGELOG.md`](CHANGELOG.md) for the full record. Brief summary:

- **v1.0.0** — seven rules, 9 sentinels, inaugural Ipi1 audit
- **v1.0.1** — live UniProt + ChEMBL adapters; ruleset unchanged
- **v1.1.0** — R7 substantive-caveat upgrade; new sentinel; new SHA
- **v1.2.0** — first external audit (TYK2); CLI runner; adapter bug fix
- **v1.2.1** — public artifact release; documentation
- **v1.3.0** — R6 scope expansion + paralog-ratio heuristic; new SHA
- **v1.3.1** — TYK2 v1.3.0 re-audit narrative; honest negative result
- **v1.4.0** — *(this release)* packaging, portability, CLI, contributor docs

---

## Citation

If you use this in research, cite the specific **ruleset SHA** alongside
the tool version. See [`CITATION.cff`](CITATION.cff) for BibTeX-equivalent
metadata.

```
falsifiable-targets v1.4.0 (ruleset SHA-256 35ef2b2a...)
https://github.com/crisprking/falsifiable-targets
```

---

## License

Apache License 2.0. See [`LICENSE`](LICENSE).
