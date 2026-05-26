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

## Headline result (v1.2.0)

Two audits have been run end-to-end:

| Claim | Type | Verdict | Substantive caveats |
|---|---|---|---|
| **Ipi1 / Madurella** (in-house novel antifungal target) | `novel_target` | `FALSIFIED_WITH_CAVEATS` | 2 (orthology contested, no selectivity data) |
| **TYK2 / psoriasis** (deucravacitinib, FDA-approved 2022) | `validated_mechanism` | `SURVIVED` | 0 |

These two audits use the same ruleset, the same SHA-locked rule logic, and
the same engine. The framework's first action was to attach substantive
caveats to its creator's own published headline target. The first external
audit then cleanly passed a validated, approved mechanism. That asymmetry
is the calibration finding.

## Limitations

Read [`docs/AUDIT_LIMITATIONS_v1_2.md`](docs/AUDIT_LIMITATIONS_v1_2.md)
before drawing conclusions from any audit. For TYK2 specifically, only 2
of 7 rules used live data; 3 used hand-encoded fixture values from
primary literature; 2 were not applicable to the claim type. The
`SURVIVED` verdict is correct under v1.1.0 rules but reflects a partial
audit. v1.3.0 will expand live-data coverage and the scope of `R6`.

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

The v1.1.0 ruleset SHA is locked at:

```
2f9aab7d0ebc209f62c16eb35be31bc5b65fa2eb09adc02bea5bff5176269b32
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
- **v1.2.1** - documentation release (this one); public artifact; no code changes

## Audits

- [`docs/AUDIT_IPI1_v1_1.md`](docs/AUDIT_IPI1_v1_1.md) - inaugural self-audit, ruleset v1.1.0
- [`docs/AUDIT_TYK2_v1_2.md`](docs/AUDIT_TYK2_v1_2.md) - first external audit, ruleset v1.1.0
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
