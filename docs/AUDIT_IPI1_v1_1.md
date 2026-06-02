# Ipi1 / Madurella mycetoma audit

**Claim**: Ipi1 (KXX81897.1, a Rix1-complex component required for
pre-60S ribosome biogenesis in *Madurella mycetomatis*) is a viable
novel antifungal target for Madurella mycetoma. The fungal protein
presents a high-confidence AlphaFold model (pLDDT 93.2) with a
ligand-accessible pocket (fpocket score 0.749, volume 497 Å³). No
prior chemistry exists against this gene in ChEMBL as of the audit
date (May 2026).

**Verdict**: `FALSIFIED_WITH_CAVEATS`
**Score**: 0.00 (0 falsified / 2 applicable)
**Substantive caveats**: 2 (R1 orthology, R7 selectivity)
**Operational notes**: 2

**Ruleset**: v1.1.0 (SHA `2f9aab7d0ebc209f62c16eb35be31bc5b65fa2eb09adc02bea5bff5176269b32`)
**Source of underlying claim**: [crisprking/madurella-target-discovery](https://github.com/crisprking/madurella-target-discovery), M5/M6 record

The complete machine-readable record (regenerable via `python
run_audit.py claims/ipi1_madurella.yaml --no-live --json-out
reports/ipi1_audit.json`) and the claim YAML
[`../claims/ipi1_madurella.yaml`](../claims/ipi1_madurella.yaml) are
the audit's machine artifacts; this document is the human-readable
narrative.

## Why this audit exists

The framework's first action was to audit its creator's own published
headline target. The reasoning, in one paragraph: a target-validation
framework whose inaugural public demonstration is killing other
people's work - before it has demonstrated willingness to demote its
creator's own work - has the order of legitimacy backwards. The
honest sequence is self-audit first, external audit second. The
Madurella M5/M6 record had already flagged the Ipi1 selectivity-vs-
TEX10 question as the open issue in narrative prose; this audit
encodes that question structurally and lets the framework decide
whether it can be reproduced from the encoded fixture alone.

## Verdict shape

`FALSIFIED_WITH_CAVEATS` at score 0.00 means:

- Zero rules outright falsified the claim
- At least one substantive caveat was attached - in this case two
- The claim is held to be **provisional but not invalidated**

This is exactly the verdict shape the original Madurella M6 record
described in prose. The framework reproduced the original narrative's
conditional structure from the encoded inputs, without ever reading
the narrative. That convergence is the inaugural audit's calibration
finding.

## Per-rule findings

### R1_orthology - `PASSED` with **substantive caveat** (confidence 0.5)

Fixture inputs: `sources_agreeing = 2`, `sources_total = 4`,
`ortholog_count_human_pathogen = 1`.

R1's threshold logic: strict majority of orthology sources must
agree. For 4 sources, strict majority is 3. With 2/4 agreeing, the
rule does not falsify (a falsification requires `sa == 0` and
`st >= 3`), but it does attach a substantive caveat.

The caveat text:
> orthology sources disagree (2/4, strict majority requires 3);
> selectivity claim weakened until resolved structurally

This is the framework's recovery of the M5/M6 narrative's core
observation: OrthoFinder missed the fungal-Ipi1 ↔ human-TEX10
ortholog because TEX10 has a large C-terminal extension absent from
the fungal protein, and the shared region is ~20% identical at
sequence level. Different ortholog databases handle this divergent
homology differently; some call it a clear ortholog (because the
core domain is recognisable), others miss it (because sequence-
identity thresholds reject it).

The substantive caveat is the framework's structural way of saying
"the selectivity question between fungal Ipi1 and human TEX10 needs
to be resolved by a structural argument, not by an orthology-
threshold call."

### R2_chemistry_support - `NOT_APPLICABLE`

`novel_target` claims do not invoke R2. By design: a novel target
with zero ChEMBL chemistry is *expected*, not a falsification.

### R3_genetics_support - `NOT_APPLICABLE`

Human GWAS evidence is not the relevant axis for a fungal pathogen
target.

### R4_expression - `NOT_APPLICABLE` (operational caveat)

The Ipi1 claim's fixture supplies no expression data. R4 abstains
operationally rather than substantively because expression-in-
indication-tissue is a *supporting* axis for a novel pathogen target,
not the load-bearing one. (For a human target, R4 would be more
load-bearing.)

### R5_replication - `PASSED` (confidence 0.9)

Not retracted, no formal rebuttals. The novel-target claim is from
the project's own M5/M6 work; there has been no replication attempt
because the target is novel.

### R6_chemistry_class_collapse - `NOT_APPLICABLE`

The `applies_to()` triggers for `novel_target` and `chemistry_series`,
so this rule *could* fire for Ipi1. However, R6 requires
`chembl_distinct_compounds > 0`, and the fixture supplies 0. The
rule returns NOT_APPLICABLE with an operational note - "no chemistry
to assess for class collapse."

This is the correct behavior. R6's class-collapse signal is a
falsification axis (the chemistry is phantom because it collapses
onto a paralog); zero chemistry can't collapse onto anything.

### R7_selectivity_counterscreen - `PASSED` with **substantive caveat** (confidence 0.5)

This is the v1.1.0 rule change in action. Under v1.0.0, R7 would have
returned `ABSTAINED` with only an operational caveat ("no selectivity
data; rule abstains"). Under v1.1.0, for a `novel_target` claim,
R7 emits a substantive caveat instead:

> no selectivity data for a novel-target claim; the selectivity vs
> human paralog question is unresolved and the claim remains
> provisional until counter-screen data exists

The rationale for the v1.1.0 change: for a novel target whose
selectivity vs a human paralog is the named open question, the
absence of selectivity data is not a tooling gap - it's the entire
reason the claim is provisional. The framework should surface that
in the verdict, not abstain politely.

The Ipi1 audit was the case that surfaced this gap during v1.0.0
review; the change shipped in v1.1.0 with this audit as the
regression test.

## Comparison: v1.0.0 vs v1.1.0 Ipi1 verdict

| | v1.0.0 | v1.1.0 |
|---|---|---|
| Verdict | `FALSIFIED_WITH_CAVEATS` | `FALSIFIED_WITH_CAVEATS` |
| Score | 0.00 | 0.00 |
| Substantive caveats | 1 (R1 only) | 2 (R1 + R7) |
| Operational notes | 3 | 2 |
| Reader impact | "orthology contested" | "orthology contested AND selectivity-data gap is substantive" |

The verdict didn't change. The *evidence shape* sharpened. Under
v1.1.0 the framework explicitly says both load-bearing axes of
provisionality, not just one.

## What this audit means

- The framework's calibration on a novel-target claim with explicit
  open questions produces the expected verdict shape.
- The Madurella headline result stands as conditional, not
  invalidated. The conditions are named structurally rather than
  buried in prose.
- The v1.1.0 R7 change was driven by the gap this audit surfaced,
  which is the correct order: audits drive rule changes, never the
  reverse.
- The first action the framework took was to attach substantive
  caveats to its creator's own work. Every subsequent audit (TYK2
  and beyond) inherits the methodological standing that decision
  bought.

## Reproduce locally

```bash
python run_audit.py claims/ipi1_madurella.yaml --no-live \
    --json-out reports/ipi1_audit.json
```

The audit completes in well under a second. The verdict, score, and
caveat shape will be identical to this document as long as the
ruleset SHA equals `2f9aab7d...` (the v1.1.0 lock).
