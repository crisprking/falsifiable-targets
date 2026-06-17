# falsifiable-targets · R3 — direction-of-effect audit

**Human genetic support predicts drug *targets*. It does not, on its own, predict their *direction*.**

This repository holds the direction-of-effect (DOE) component of the `falsifiable-targets` auditor: an abstention-first caller that reads Open Targets directional evidence and tries to answer the question genetic support *doesn't* — **agonize or inhibit?** — and a reproducible audit suite showing exactly where that question can be answered and where it can't.

It is a **boundary / methods result**, not a predictor. It corroborates, from a falsification-first direction, the published state of the art (Chen et al., *npj Drug Discovery*, 2025) and the Open Targets DOE assessment. It does not supersede them, and it should not be cited as a competitive direction predictor — on unselected targets it commits on 2 of 132 pairs and one of those is wrong. The contribution is the map, the named failure mode, and the discipline of refusing where the signal isn't there.

---

## The result in one table

The same 78 unselected approved-drug targets (11 common diseases, drug mechanism used *only* as an answer key) scored under three trait-anchorings, source set and mapping held fixed. All figures from live Open Targets runs.

| Anchoring | Coverage | Accuracy (committed) | Always-inhibitor baseline | Discrimination |
|---|---:|---:|---:|---|
| **Exact indication** (the tool's logic) | **1.5 %** (2/132) | 50 % (1/2), Wilson [10, 90] | 50 % | none — MCC 0.00 |
| Indication + descendants | 4 % | 40 % | — | none — MCC 0.00 |
| **Pooled across all traits** (the shortcut) | **34 %** (45/132) | **33 %**, Wilson [21, 48] | **73 %** | one-directional, *below chance* |
| *Curated biomarker panel* | *≈ 7 %* | *8/8, Wilson [0.68, 1.00]* | — | *both directions* |

Read top to bottom, that is the whole story. At the disease indication, **130 of 132 pairs carry no directional genetic row at all** — the tool isn't adjudicating ambiguity, it's finding nothing to adjudicate. The one place coverage climbs (pooling each target's genetics across every associated trait) buys it by dredging loss-of-function rare-disease evidence that points *opposite* the inhibitor drug, landing **40 points below the baseline of always guessing "inhibitor."** Direction is recoverable with both coverage and correctness only where the genetic and therapeutic arms share a trait — biomarker-proximal targets, the favourable corner of a small domain.

---

## Why pooling inverts the sign (attributed three ways)

The pooled caller is, to first order, an **"always say activator"** machine: per-class recall **{inhibitor 0.09, activator 1.00}**. The cause is that *direction is trait-relative, not target-intrinsic.* The directional evidence pooled across a target's traits is dominated by *loss-of-function-causes-rare-disease* rows annotated at a Mendelian condition — not the common-disease indication — whose sign is correct *there* and inverts when read against an inhibitor drug treating overactivity at the common disease.

The cleanest teaching example is **SOST**: its only directional evidence sits at *sclerosteosis*, so its sign is right for sclerosteosis and wrong for osteoporosis — and the indication-anchored tool correctly goes **silent** on it. Same pattern for ACE (renal tubular dysgenesis), TNFSF11/RANKL (osteopetrosis), CETP (CETP deficiency), P2RY12 (platelet bleeding disorder), and JAK2 (whose drug exploits the GoF V617F, not the LoF-disease arm).

Three independent checks prove this is the data's structure, not a bug in the mapping:

- **(A) The mapping is Open Targets' own convention** — `desired_from_label` matches OT's documented logic 4-for-4, anchored on their PCSK9 example (LoF protective → inhibit).
- **(B) Gene-burden assumes loss-of-function by construction** — **80/80 burden rows in the set are `directionOnTarget = LoF`**, so any risk-increasing burden gene is *forced* to "activator." The skew is structural to OT's methodology, not chosen here.
- **(C) The same mapping recovers truth where the frame aligns** — at the LDL biomarker (EFO_0004611), four inhibitor-drug lipid anchors return **unanimous inhibitor, zero conflict**: PCSK9 (32 burden votes), ANGPTL3 (15), NPC1L1 (11), HMGCR (1).

A disease-cluster bootstrap (B = 2000, resample the 11 diseases with replacement, seed 20240617) shows the residual is fragile: pooled accuracy **33.3 % [18.9, 44.0]** sits entirely below the baseline **73.3 % [60.9, 86.3]** (intervals don't overlap — below-chance is not a panel artifact), and the **MCC interval [+0.000, +0.249] includes zero.** The faint residual is carried by three mechanistically frame-aligned genes (PCSK9, PDE3A, EDNRA); drop their diseases and it collapses.

---

## What the tool actually does

It reads only Open Targets sources with intrinsically directional study design — **gene-burden** (`directionOnTarget` ∈ {LoF, GoF}, `directionOnTrait` ∈ {risk, protect}) and **curated clinical-variant** sources (eva/ClinVar, ClinGen, Genomics England, Orphanet, Gene2Phenotype, UniProt). GWAS credible sets carry no direction label on the evidence row and are excluded by construction.

Each directional row is mapped to a therapeutic call by **genetic mimicry** — a drug should mimic the protective genetic direction:

| directionOnTarget | directionOnTrait | call |
|---|---|---|
| LoF | protect | **inhibitor** |
| LoF | risk | **activator** |
| GoF | risk | **inhibitor** |
| GoF | protect | **activator** |

The caller commits only when its strongest *speaking* tier (burden, then clinical) is internally unanimous; it abstains `CONFLICTED` on any in-tier split and `NO_LABEL` when no directional row exists. A circularity guard is asserted at runtime — `assert GT_SOURCE.isdisjoint(NONCIRC)` — so the answer key (the approved drug's mechanism, `clinical_precedence`) is provably **disjoint** from the predictor source set. Approval can never vouch for itself.

---

## Reproduce

The audit suite lives in the notebook (`notebookef00f61107.ipynb`); every result maps to a version-pinned cell with a runtime assertion where one applies.

| Cell | Produces | Saved artifact |
|---|---|---|
| 76 | trait-scope × decision-rule grid | `scope_rule_grid.json` |
| 77 / 78 | mapping audit A–D (convention, LoF fraction, LDL positive control, residual) | `mapping_audit.json` |
| 79 | disease-cluster bootstrap (B = 2000) | `panel_bootstrap.json` |
| 75 | indication-anchored score | `indication_anchored.json` |

```bash
# the figure (reads scope_rule_grid.json — nothing hardcoded, fails loud if absent)
python fig_boundary_surface.py
```

Cells are self-contained: each re-bootstraps the engine and re-harvests before scoring, so any one runs from a cold kernel (Kaggle wipes `/kaggle/working` between sessions). Public data only — Open Targets, ChEMBL (via `clinical_precedence`), ClinVar. **Pin the Open Targets release** before a camera-ready run; the numbers above are tied to the release queried.

---

## What is and isn't established

| Claim | Status |
|---|---|
| The tool abstains rather than miscommit when aligned directional evidence is absent | **Supported** (130/132 `NO_LABEL`; 17/18 known traps abstain; bootstrap-robust) |
| Naive trait-pooling is below the majority baseline and systematically one-directional | **Supported** (acc [18.9, 44.0] entirely below baseline [60.9, 86.3]) |
| The inversion is a LoF-disease / trait-frame artifact, not a bug | **Supported** (mapping ≡ OT convention; 80/80 burden LoF; LDL positive control) |
| Direction is recoverable where genetic + therapeutic arms align at the trait | **Supported (narrow)** (curated panel 8/8; ≈ 5–7 % of targets) |
| The tool accurately calls direction on **unselected** targets | **Unmeasurable** (coverage 1.5 %, n = 2; only 1797/2000 resamples commit ≥ 1) |
| Anchoring to the indication **fixes** the inversion | **Refuted as stated** (IL17RA commits wrong at psoriasis vs brodalumab) |

---

## Positioning

This is a **selective-prediction / cautionary** contribution, orthogonal to the predictor literature rather than competitive with it. The target-level enrichment of genetic support is established (Nelson 2015; King 2019; Minikel 2024, 2.6×). For the *direction* problem, Chen et al. (2025) report macro-AUROC **0.59** for gene–disease-specific DOE across 47,822 pairs — and independently report the same ~76 % inhibitor class imbalance. Those numbers don't threaten this story; they **validate** it. The unoccupied niche is the auditor's, not the predictor's: a tool that **refuses where the state of the art is at chance** is the complement their 0.59 implies.

The two genuinely new, methods-level contributions: **(i)** a named and quantified failure mode — naive trait-pooling of platform directional evidence is below a majority baseline and one-directional; **(ii)** an interpretable mechanism for it — the trait-relative loss-of-function-disease confound, attributed three ways and named down to specific genes — packaged as an abstention-first design with a hard circularity guard, where the coverage map *is* the deliverable.

Full write-up: [`direction_of_effect_boundary.md`](./direction_of_effect_boundary.md).

---

## Cite

> Genetic support predicts drug targets, not their direction: a calibrated-abstention audit of Open Targets direction-of-effect, and the trait-pooling inversion. *Preprint, 2026.* github.com/crisprking/falsifiable-targets

### Key references

- Nelson MR, et al. *Nat Genet* 47, 856–860 (2015).
- King EA, Davis JW, Degner JF. *PLoS Genet* 15, e1008489 (2019).
- Minikel EV, Painter JL, Dong CC, Nelson MR. *Nature* 629, 624–629 (2024). doi:10.1038/s41586-024-07316-0
- Duffy Á, et al. *Nat Genet* 56, 51–59 (2024).
- Chen R, Duffy Á, Park JK, et al. (Do R, senior author). *npj Drug Discovery* 2, 24 (2025). doi:10.1038/s44386-025-00027-0
- Buniello A, et al. *Nucleic Acids Res* 53, D1467–D1475 (2025).

---

*Public-data-only by design. The auditor demotes its author's own headline target before it vouches for an FDA approval — the credibility is in what it refuses, not what it claims.*
