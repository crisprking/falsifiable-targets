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
  *(repo URL will be live after v1.2.1 push)*
- **Ruleset version**: v1.1.0
- **Ruleset SHA-256**:
  `2f9aab7d0ebc209f62c16eb35be31bc5b65fa2eb09adc02bea5bff5176269b32`
- **Claim file**: `claims/ipi1_madurella.yaml` in the framework repo
- **Full narrative**: [`docs/AUDIT_IPI1_v1_1.md`](https://github.com/crisprking/falsifiable-targets/blob/main/docs/AUDIT_IPI1_v1_1.md)

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
same ruleset. The asymmetry — substantive caveats on the in-house
novel target, zero substantive caveats on the FDA-approved external
mechanism — is documented as the framework's v1.2.x calibration result.
