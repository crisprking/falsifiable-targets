# M7 — Falsifiability audit of the Ipi1 headline finding

This module records the result of auditing this repo's headline novel-
target finding (Ipi1, KXX81897.1, Rix1 complex) through an external,
SHA-locked target-validation framework.

## Result

| Field | Value |
|---|---|
| Verdict | `FALSIFIED_WITH_CAVEATS` |
| Score | 0.00 (0 rules falsified / 2 rules applicable) |
| Substantive caveats | 2 |
| Operational notes | 2 |

## Audit traceability

- **Framework**: [falsifiable-targets](https://github.com/crisprking/falsifiable-targets)
- **Ruleset version**: v1.2.0 (active in v1.3.0+ releases)
- **Ruleset SHA-256**:
  `35ef2b2ab5363298097962a0b6ae52c70d551a1edddc341054f75cb6e4fb7221`
- **Claim file**: `claims/ipi1_madurella.yaml` in the framework repo
- **Full narrative**: [`docs/AUDIT_IPI1_v1_1.md`](https://github.com/crisprking/falsifiable-targets/blob/main/docs/AUDIT_IPI1_v1_1.md)

Note: The Ipi1 audit verdict is unchanged across ruleset versions
(v1.1.0 → v1.2.0). The v1.2.0 ruleset's changes affect R6 scope and
paralog-overshadow detection, which do not apply to the Ipi1 claim
(novel target, no chemistry). All Ipi1 verdict and caveat results
under v1.2.0 are byte-identical to v1.1.0.

## What the framework caught

The audit reproduced, from a six-field structured fixture alone, the
two open questions this repo's M5/M6 narrative had flagged in prose:

1. **R1 orthology** — substantive caveat. Source-disagreement on the
   Ipi1↔TEX10 ortholog call (2 of 4 databases agree; strict majority
   requires 3) means the selectivity argument cannot be assumed from
   orthology thresholds alone and must be resolved structurally.
2. **R7 selectivity counter-screen** — substantive caveat. No
   selectivity data exists for the Ipi1 series (no series exists yet),
   so the selectivity-vs-TEX10 question is unresolved at the
   chemistry-axis level until counter-screen data is produced.

The framework recovered exactly the conditional structure the original
narrative described. No prior knowledge of the narrative was supplied
to the framework; only the structured fixture in
`claims/ipi1_madurella.yaml`.

## Interpretation

The headline result of this repo stands as **conditional, not
invalidated**. The two conditions are now named in machine-readable,
SHA-locked form. Resolving them - via structural comparison of the
fungal Ipi1 binding pocket with the human TEX10 binding pocket, and
via selectivity counter-screening of any future Ipi1 chemistry against
TEX10 - is the work that converts the headline from
`FALSIFIED_WITH_CAVEATS` to `SURVIVED` under the same ruleset.

## Reproducing the audit

```bash
git clone https://github.com/crisprking/falsifiable-targets.git
cd falsifiable-targets
python run_audit.py claims/ipi1_madurella.yaml --no-live \
    --json-out /tmp/m7_audit.json
```

The audit is fully deterministic given the locked ruleset SHA. Output
will match the values above as long as the ruleset SHA is
`2f9aab7d0e...`.

## Cross-reference

The framework that produced this audit is itself audited: the
framework's first action was to audit *this* repo's claim. The first
external audit it performed (TYK2 / psoriasis) cleanly survived the
same ruleset. Under v1.3.0 the framework was upgraded to also audit
the pan-JAK class-collapse risk for TYK2; the re-audit found no
pool-size overshadow (TYK2: 538 ChEMBL compounds; JAK1/2/3: 549/772/832;
max ratio 1.55× < 2.0× threshold). See [`docs/AUDIT_TYK2_v1_3.md`](https://github.com/crisprking/falsifiable-targets/blob/main/docs/AUDIT_TYK2_v1_3.md)
for the honest distinction between pool-size overshadow (what R6 v1.2.0
audits) and compound-level overlap (v1.4 milestone).

The asymmetry remains the framework's headline result: substantive
caveats on the in-house novel target, zero substantive caveats on the
FDA-approved external mechanism under deeper-audit conditions.
