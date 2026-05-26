# falsifiable-targets

A target-validation audit framework. Encodes seven rules that take a target
claim (validated mechanism, novel target, chemistry series, or extraordinary
claim) and returns a deterministic verdict: `SURVIVED`,
`FALSIFIED_WITH_CAVEATS`, `FALSIFIED`, or `INSUFFICIENT_DATA`. Every audit is
SHA-locked: the ruleset is content-addressed, the claim is content-addressed,
the report names both hashes.

The framework exists because too many published target claims are
defended by evidence that, on inspection, points at a different target -
the SAT/HDAC4 archetype from the Madurella audit, the STAP-cells failure
mode, the CETP clinical-genetics gap. The seven rules are designed to
detect those failure modes from public data alone, at the cheapest
falsification tier available.

## Headline result (v1.3.0)

Three audits have been run end-to-end:

| Claim | Type | Ruleset | Verdict | Substantive caveats |
|---|---|---|---|---|
| **Ipi1 / Madurella** (in-house novel antifungal target) | `novel_target` | v1.1.0 | `FALSIFIED_WITH_CAVEATS` | 2 (orthology contested, no selectivity data) |
| **TYK2 / psoriasis** (deucravacitinib, FDA-approved 2022) | `validated_mechanism` | v1.1.0 | `SURVIVED` | 0 |
| **TYK2 / psoriasis** (re-audit under expanded R6) | `validated_mechanism` | v1.2.0 | `SURVIVED` | 0 |

The Ipi1 vs TYK2 asymmetry is the calibration finding: same engine,
same SHA-locked rule logic, different verdicts in the direction the
structure of each claim predicts. The v1.3.0 TYK2 re-audit is
qualitatively stronger than the v1.1.0 one — R6 now applies to
`validated_mechanism` claims and was actually evaluated against live
ChEMBL paralog data (TYK2: 538 compounds; JAK1: 549; JAK2: 772; JAK3:
832; max ratio 1.55× < 2.0× threshold). The pan-JAK class-collapse
risk that v1.2.x didn't even check is now structurally audited and
returns no overshadow at the pool-size axis.

The v1.3.0 TYK2 SURVIVED is *not* a clean bill of selectivity — see
[`docs/AUDIT_TYK2_v1_3.md`](docs/AUDIT_TYK2_v1_3.md) for the honest
distinction between *pool-size overshadow* (what R6 audits) and
*compound-level overlap* (v1.4 milestone).

## Limitations

Read [`docs/AUDIT_LIMITATIONS_v1_3.md`](docs/AUDIT_LIMITATIONS_v1_3.md)
before drawing conclusions from any audit. For TYK2 specifically under
v1.3.0, 3 of 6 applicable rules used live data; 3 used hand-encoded
fixture values from primary literature. The v1.3.0 R6 closes the
class-collapse scope gap from v1.2.x but operates on *pool-size
overshadow*, not *compound-level overlap* — the deucravacitinib
JH2-domain selectivity argument lives at the latter axis, which is
queued for v1.4.

## Quick reproduce

```bash
# In any Python 3.10+ environment with PyYAML and pytest installed
python smoke_test.py         # 10/10 sentinel calibration suite
python -m pytest tests/      # full test suite (31 tests)
python run_audit.py claims/ipi1_madurella.yaml --no-live
python run_audit.py claims/tyk2_psoriasis.yaml --json-out reports/my_run.json
```

The `--no-live` flag runs against fixture data only (hermetic, no
network). Without `--no-live`, the audit fetches UniProt and ChEMBL
live and caches responses under `.ae_cache/`.

## Reproducing the audits exactly

The v1.2.0 ruleset SHA is locked at:

```
35ef2b2ab5363298097962a0b6ae52c70d551a1edddc341054f75cb6e4fb7221
```

Every audit report (`reports/*.json`) stamps this SHA alongside the
claim's own SHA. Two audits with matching ruleset SHA and matching
claim SHA must produce byte-identical verdicts; if they do not, the
adapter layer changed but the rules did not - read
[`docs/AUDIT_LIMITATIONS_v1_2.md`](docs/AUDIT_LIMITATIONS_v1_2.md) §3.

## Release history

See [`CHANGELOG.md`](CHANGELOG.md). Brief summary:

- **v1.0.0** - seven rules, 9 sentinels, inaugural Ipi1 audit
- **v1.0.1** - live UniProt + ChEMBL adapters added, ruleset unchanged
- **v1.1.0** - R7 substantive-caveat upgrade for novel-target claims; new sentinel; new ruleset SHA
- **v1.2.0** - first external audit (TYK2 / psoriasis); CLI runner; adapter bug fix
- **v1.2.1** - documentation release; public artifact; no code changes
- **v1.3.0** - R6 scope expansion to validated_mechanism + paralog-ratio heuristic; new ruleset SHA; new sentinel; TYK2 expected to shift to FALSIFIED_WITH_CAVEATS under live mode

## Audits

- [`docs/AUDIT_IPI1_v1_1.md`](docs/AUDIT_IPI1_v1_1.md) - inaugural self-audit, ruleset v1.1.0
- [`docs/AUDIT_TYK2_v1_2.md`](docs/AUDIT_TYK2_v1_2.md) - first external audit, ruleset v1.1.0
- [`docs/AUDIT_TYK2_v1_3.md`](docs/AUDIT_TYK2_v1_3.md) - TYK2 re-audit under v1.2.0 (R6 scope expanded; verdict unchanged but evidence strictly deeper)
- [`reports/`](reports/) - JSON artifacts (machine-readable, SHA-stamped)

## The seven rules at a glance

| Rule | Audits | Cheapest tier |
|---|---|---|
| R1_orthology | Pathogen/novel targets: do >= majority of ortholog DBs agree the protein is present? | public data lookup |
| R2_chemistry_support | Chemistry-series/validated claims: are there meaningful ChEMBL compounds? | public data lookup |
| R3_genetics_support | Validated-mechanism claims: GWAS, Mendelian, or somatic-driver evidence? | public data lookup |
| R4_expression | Target detectably expressed in indication-relevant tissue? | public data lookup |
| R5_replication | Not retracted; no overwhelming rebuttals without replication | public data lookup |
| R6_chemistry_class_collapse | Do ChEMBL hits collapse onto a paralog class (phantom evidence)? | public data lookup |
| R7_selectivity_counterscreen | Does selectivity data exist; for novel targets, absence is substantive | cheap in silico |

Falsification tiers ordered: `public_data_lookup` < `cheap_in_silico` <
`targeted_assay` < `cohort_study` < `clinical_trial`. The aggregator
always names the *cheapest* available experiment that would resolve a
falsification, because falsifications cheap enough to run from a laptop
are the framework's primary value proposition.
