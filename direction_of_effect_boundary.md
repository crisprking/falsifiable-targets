# Genetic support predicts drug targets, not their direction

### A calibrated-abstention audit of Open Targets direction-of-effect, and the trait-pooling inversion

**A boundary/methods result. Corroborates Chen et al. (*npj Drug Discovery*, 2025) and the Open Targets direction-of-effect assessment; does not supersede them.**

---

## Abstract

Human genetic support for a drug *target* is one of the most reliable predictors of clinical success — but support is not direction. Whether to **agonize or inhibit** a supported target is a separate question, and it is the one that determines therapeutic success. We audit how far genetic associations alone answer it, using an abstention-first caller that reads only intrinsically directional, non-circular evidence (rare-variant burden plus curated Mendelian direction-of-effect), anchors to a single indication, and refuses when aligned evidence is absent.

On an unselected set of 132 approved-drug (target, indication) pairs across 11 common diseases (78 unique targets, drug mechanism used only as an answer key), the caller **commits on 1.5 %** of pairs and abstains on the rest, because at the disease indication a directional genetic signal is essentially absent: of 132 pairs, 130 carry no directional burden/clinical row at the indication trait. A disease-cluster bootstrap (B = 2000) puts this at coverage 1.5 % [0.0, 3.7], with only 1797/2000 resamples committing even a single pair — so accuracy on unselected targets is not low, it is **unmeasurable**.

The naive shortcut — pooling Open Targets directional evidence across *all* of a target's associated traits to manufacture coverage — raises coverage to 34 % but yields **33 % accuracy [18.9, 44.0]**, which sits entirely below the always-inhibitor majority baseline of **73 % [60.9, 86.3]** and is systematically one-directional (true-activator recall 1.00, true-inhibitor recall 0.09). We attribute the inversion three independent ways: the direction mapping is Open Targets' own documented convention; gene-burden evidence carries a built-in loss-of-function assumption (100 % of burden rows in our set), which forces any risk-increasing burden gene to "activator"; and the same mapping correctly recovers "inhibitor" for four lipid anchors (PCSK9, ANGPTL3, NPC1L1, HMGCR) at the LDL biomarker where the genetic and therapeutic arms share a trait. The residual signal at pooled scope is carried by three mechanistically frame-aligned targets (PCSK9, PDE3A, EDNRA) and its MCC interval includes zero.

The usable result is a **boundary, not a predictor**: human genetics speaks to drug-target direction only where the genetic and therapeutic arms align at the queried trait — biomarker-proximal targets (≈ 5–7 % of the set, where a curated panel scores 8/8) — and a tool that abstains elsewhere is doing the correct thing, while one that pools to manufacture coverage scores below chance. This independently reproduces, from a falsification-first direction, the state of the art for the disease-anchored problem (macro-AUROC 0.59; Chen et al. 2025).

---

## 1. The question, and why it is not the usual one

The headline result the field trusts is target-level: drug mechanisms with human genetic support succeed at roughly 2.6× the rate of those without, an enrichment that rises with confidence in the causal gene and is largely independent of effect size or allele frequency (Minikel et al. 2024); when causal genes are unambiguous — Mendelian traits and coding-linked GWAS — the approval enrichment exceeds two-fold and holds prospectively for Mendelian associations (King et al. 2019), extending the foundational doubling first reported by Nelson et al. (2015).

That number is real on average and it is about **whether a target is causal**. It does not tell you the **direction of modulation** — whether the drug should be an agonist or an antagonist — for any single target. Direction is the harder problem, and the recent state of the art quantifies exactly how hard: a dedicated framework predicting direction-of-effect (DOE) reports macro-AUROC 0.95 for DOE-specific druggability across 19,450 genes and 0.85 for isolated DOE among 2,553 druggable genes, but only **0.59 for gene–disease-specific DOE across 47,822 pairs**, with performance improving where more genetic evidence is available (Chen et al. 2025). The disease-anchored problem — the one a drug program actually faces — is near a coin even for a full machine-learning model with embeddings.

This work asks a complementary question. Not "can we push the 0.59 higher?" — that is an arms race against an embeddings team and the wrong contribution for a public auditor — but: **on an unselected target set, where does genetic direction exist at all, how does the obvious coverage shortcut fail, and is the correct behavior, when aligned evidence is absent, to abstain?**

---

## 2. The caller, in one paragraph

The tool reads only Open Targets evidence sources with intrinsically directional study design — gene-burden (which carries `directionOnTarget` ∈ {LoF, GoF} and `directionOnTrait` ∈ {risk, protect}) and curated clinical-variant sources (`eva`/ClinVar, ClinGen, Genomics England, Orphanet, Gene2Phenotype, UniProt). GWAS credible sets, which on the evidence row carry no direction label, are excluded by construction. Each directional row is mapped to a therapeutic direction by **genetic mimicry** — a drug should mimic the protective genetic direction — and the caller commits only when its strongest *speaking* tier (burden, then clinical) is internally unanimous, abstaining (`CONFLICTED`) on any in-tier bidirectional split and `NO_LABEL` when no directional row exists. A hard circularity guard is asserted at runtime: the answer key (the approved drug's mechanism, `clinical_precedence`) is **disjoint** from the predictor source set, so approval can never vouch for itself.

The genetic-mimicry mapping, encoding GoF = +1 / LoF = −1 and risk = +1 / protect = −1 with inhibitor iff the product is positive:

| Direction on target | Direction on trait | Therapeutic call | Canonical example |
|---|---|---|---|
| LoF | protect | **inhibitor** | PCSK9 (LoF lowers LDL, protective → inhibit) |
| LoF | risk | **activator** | tumour-suppressor / restore-function |
| GoF | risk | **inhibitor** | block the gained activity |
| GoF | protect | **activator** | |

This is Open Targets' own convention, not ours — established in §4 below.

---

## 3. Data and method

**Unselected set.** Phase 1 harvests, from a fixed 11-disease panel (coronary artery disease, hypertension, rheumatoid arthritis, Crohn disease, inflammatory bowel disease, psoriasis, COPD, breast carcinoma, osteoporosis, ulcerative colitis, atopic dermatitis), the top-20 Open Targets–associated targets per disease that also carry an approved-drug mechanism — yielding **132 (target, indication) pairs over 78 unique targets**. The drug mechanism is read once, from `clinical_precedence`, as the **answer key**; it is never an input to the prediction.

**Three anchorings, one set.** The same 78 targets are scored under three trait-scope choices, holding the source set and the mapping fixed:

- **Exact indication** — direction read at the drug's own indication trait (the tool's real logic).
- **Indication + descendants** — broadened to sub-phenotype EFO descendants.
- **Pooled** — each target's associated diseases auto-discovered and genetic direction pooled across all of them (the coverage shortcut).

crossed with three decision rules — tier-precedence (the tool's), strict-unanimous, and majority-vote — giving a 3 × 3 surface.

**Metrics.** Coverage (fraction committed), accuracy on committed pairs with Wilson 95 % intervals, balanced accuracy, and Matthews correlation coefficient (MCC) — the last two because the population is ≈ 73 % inhibitor, so a constant predictor scores high raw accuracy on the majority class and MCC ≈ 0 is the honest floor for "no discrimination." All numbers below are from live Open Targets runs; the disease-cluster bootstrap (§7) resamples the 11 diseases with replacement (B = 2000, fixed seed) to test robustness to panel choice.

**Provenance.** Every claim maps to a specific, version-pinned notebook cell with a runtime assertion where one applies: the mapping (engine `desired_from_label`, unit-tested); the circularity guard (`assert GT_SOURCE.isdisjoint(NONCIRC)`); the pooling vs. anchoring contrast (the pooled cell auto-discovers `associatedDiseases`; the indication cell anchors to a single EFO). Public data only: Open Targets, ChEMBL (via `clinical_precedence`), ClinVar.

---

## 4. Result 1 — the landscape: direction is sparse and biomarker-proximal

Scored on the same unselected 78-target set, the three anchorings give a clean monotone in coverage at a steep cost in correctness:

**Table 1. The three-anchoring landscape (exact figures from live runs).**

| Anchoring | Coverage | Accuracy (committed) | Majority baseline | Discrimination |
|---|---:|---:|---:|---|
| **Exact indication** (tool's logic) | **1.5 %** (2/132) | 50 % (1/2), Wilson [10, 90] | 50 % | none (MCC 0.00) |
| Indication + descendants | 4 % | 40 % | — | none (MCC 0.00) |
| **Pooled** (coverage shortcut) | **34 %** (45/132) | **33 %**, Wilson [21, 48] | **73 %** | one-directional (below) |
| *Curated biomarker panel* | *≈ 7 %* | *8/8, Wilson [0.68, 1.00]* | — | *both directions* |

The exact-indication corner — the tool as designed — commits on **2 of 132 pairs**. At the disease indication, **130 of 132** pairs carry no directional burden/clinical row at all (`NO_LABEL`); zero are internal conflicts. The tool is not adjudicating ambiguity, it is finding nothing to adjudicate. This is the principled silence the design claims, now shown on an unselected set rather than asserted: where the genetic and therapeutic arms do not share a trait, the directional signal is simply not annotated.

The only regime where direction is recoverable with both coverage and correctness is **biomarker-proximal**: a curated, drug-sourced panel of lipid/metabolic targets scores 8/8 (PCSK9, APOB, ANGPTL3, HMGCR, NPC1L1, LDLR, APOC3, GCK), discriminating **both** directions (LDLR and GCK are activator targets; the rest inhibitor). That panel is ≈ 5–7 % of target space once you anchor to the disease — the favourable corner of a small domain, exactly as a curated-panel critique would predict.

---

## 5. Result 2 — the trait-pooling inversion: coverage bought below chance

Pooling is the obvious way to manufacture the coverage the indication-anchored tool lacks, and it is actively harmful. At pooled scope the caller commits on 45 pairs at **33 % accuracy [18.9, 44.0]**, **40 points below** the always-inhibitor majority baseline of 73 %. The error is not noise but a direction: the per-class recall at the pooled × majority corner is **{inhibitor: 0.09, activator: 1.00}**. Pooling is, to first order, an "always say activator" caller — it is right on essentially every true activator and wrong on ≈ 91 % of true inhibitors. The 54.5 % balanced accuracy and +0.16 MCC are the *signature* of near-constant-class calling plus a faint residual, not evidence of discrimination.

The decision rule barely matters. Across every row of the 3 × 3 surface, tier-precedence, strict-unanimous, and majority-vote are nearly identical (2/2/2 % at exact, 4/4/4 % at descendants, 33/31/34 % at pooled). The rule bites only when a target's evidence is internally split — and targets almost never are: each is a one-directional flood, so there is nothing to adjudicate. **Trait scope is essentially the only knob.** This is simpler and more defensible than a "the vote hides the minority" account: broadening the trait frame, not the aggregation rule, is what inverts the sign.

This is the load-bearing cautionary result. The existing literature builds DOE predictors that work to the extent the data allows; a representative multi-omics benchmark states plainly that it does not account for the directionality of therapeutic versus genetic effects. **None characterize naive trait-pooling of platform directional evidence as below a majority baseline and systematically one-directional.** Anyone querying Open Targets directly for direction will reach for this shortcut; this quantifies how, and by how much, it breaks.

---

## 6. Result 3 — the mechanism, attributed three ways

The inversion has a single, teachable cause: **direction is trait-relative, not target-intrinsic.** The directional rows pooled across a target's associated traits are dominated by *loss-of-function-causes-rare-disease* evidence annotated at a **different** trait — usually a Mendelian condition — than the common-disease indication. Its sign is correct *for that trait* and inverts when read against an inhibitor drug treating overactivity at the common disease. Worked examples from the diagnostic dump, each a true-inhibitor target that pooling miscalls "activator": ACE (directional evidence at renal tubular dysgenesis), TNFSF11/RANKL (recessive osteopetrosis), SOST (sclerosteosis), CETP (CETP deficiency), P2RY12 (platelet bleeding disorder), JAK2 (the drug exploits the GoF V617F, not the LoF-disease arm). SOST is the sharpest: its only directional evidence sits at sclerosteosis, so its sign is correct for sclerosteosis and inverted for osteoporosis — and the indication-anchored tool correctly goes silent on it.

Three independent checks establish that this is the data's structure, not our code:

**(A) The mapping is Open Targets' convention.** The engine's `desired_from_label` matches, four-for-four, the convention stated in Open Targets' own documentation — anchored on their PCSK9 example, where loss of function is protective and therefore implies inhibition — with the contrapositive and the gain-of-function cases completing the table. The Open Targets 2025 platform paper states the therapeutic logic directly: where increased protein activity confers protection, an agonist is the suitable hypothesis.

**(B) Gene-burden assumes loss-of-function — so the activator skew is structural.** Open Targets' gene-burden direction-of-effect method explicitly assumes all variants are loss-of-function, scoring β > 0 (or OR > 1) as risk and β < 0 (or OR < 1) as protective. In our set, **100 % of gene-burden rows (80/80)** carry `directionOnTarget = LoF`. Because burden assumes LoF by construction, **any gene whose rare variants raise disease risk is forced to "activator"** — the skew is baked into the platform's methodology, not chosen by us. (Open Targets further reports burden directions are stable, agreeing across collapsing models in 98 % of continuous and 87 % of binary traits — so this is the trait frame, not sampling noise.)

**(C) The same mapping recovers truth where the frame aligns.** At the LDL biomarker trait (EFO_0004611), four inhibitor-drug lipid anchors return **unanimous inhibitor** with zero conflict — PCSK9 (32 burden votes), ANGPTL3 (15), NPC1L1 (11), HMGCR (1). Where genetics and therapy share a trait, the mapping is correct. The inversion everywhere else is therefore provably the trait frame, not the engine.

---

## 7. Result 4 — the residual is three frame-aligned genes, and it is fragile

The faint +0.16 MCC at pooled scope is not weak discrimination; it is three targets where the disease-causing genetic direction happens to coincide with the therapeutic direction at the queried trait. Of 64 true-inhibitor targets, pooling calls only **3 correctly** (PCSK9, PDE3A, EDNRA) and inverts 21 to activator. The three are **mechanistically frame-aligned**, in two flavours: PCSK9 is loss-of-function-protective at a biomarker; PDE3A and EDNRA are gain-of-function disease genes whose inhibitor drugs reverse the gain (activating PDE3A mutations cause hypertension with brachydactyly; EDNRA is an activation-driven vascular receptor blocked by endothelin antagonists). In each, the genetic and therapeutic arms point the same way at the trait. (The PDE3A interpretation is literature-grounded; EDNRA's exact driving rows should be dumped for the camera-ready version.)

So the residual and the curated 8/8 are the **same phenomenon viewed twice**: the tool has signal exactly and only where the arms align. A disease-cluster bootstrap confirms the residual is not robust — resamples that drop these genes' diseases collapse it:

**Disease-cluster bootstrap (B = 2000, resample 11 diseases with replacement, fixed seed):**

- **Exact × tier:** coverage **1.5 % [0.0, 3.7]**; accuracy **unmeasurable** — only 1797/2000 resamples commit ≥ 1 pair. Indication-anchored accuracy on unselected targets cannot be estimated, and that is the honest statement, with a number behind it.
- **Pooled × majority:** coverage 34 % [26, 46]; accuracy **33 % [18.9, 44.0]**; majority baseline **73 % [60.9, 86.3]**; balanced accuracy 54.5 % [50.0, 60.0]; **MCC +0.16 [+0.000, +0.249]**.
- **Below-baseline holds across all resamples:** the upper bound of pooled accuracy (44.0 %) sits *below* the lower bound of the baseline (60.9 %). The two intervals do not overlap. "Below chance" is not a panel artifact.
- **MCC interval includes zero**, and balanced accuracy bottoms at the 50 % constant-predictor floor — no robust direction discrimination, only frame coincidence.

---

## 8. The at-indication failure mode: IL17RA

Anchoring to the indication is not a guarantee of correctness when the tool *does* commit. IL17RA is the one pooled-miss that did not abstain at its own indication: at psoriasis, the burden/clinical evidence reads (LoF → risk → activator) while brodalumab is an inhibitor. The loss-of-function-disease-arm inversion can therefore occur **at the indication itself**, not only at remote rare diseases. It is a single case in this set, but it bounds the claim: indication-anchoring greatly reduces the inversion (17 of the 18 pooled-traps now correctly abstain) without eliminating it. This belongs in the limitations, not hidden.

---

## 9. Claim ledger

The discipline of this work is to state exactly what the data supports and no more.

**Table 2. What is and is not established.**

| Claim | Status | Evidence |
|---|---|---|
| The tool abstains rather than miscommit when aligned directional evidence is absent | **Supported** | 130/132 `NO_LABEL` at the indication; 17/18 known pooling-traps abstain; robust under bootstrap |
| Naive trait-pooling is below the majority baseline and systematically one-directional | **Supported** | accuracy 33 % [18.9, 44.0] entirely below baseline [60.9, 86.3]; per-class recall {inhibitor 0.09, activator 1.00} |
| The inversion is a loss-of-function-disease / trait-frame artifact, not a bug | **Supported** | mapping ≡ OT convention (A); 80/80 burden rows LoF (B); positive control recovers inhibitor at LDL (C) |
| Direction is recoverable where genetic + therapeutic arms align at the trait | **Supported (narrow)** | curated biomarker panel 8/8, Wilson [0.68, 1.00]; ≈ 5–7 % of targets |
| The tool accurately calls direction on **unselected** targets | **Unmeasurable** | coverage 1.5 %, n = 2 commits; only 1797/2000 resamples commit ≥ 1 — uninterpretable |
| Anchoring to the indication **fixes** the inversion | **Refuted as stated** | IL17RA commits wrong at its own indication (psoriasis) |

The landscape is now complete across all three anchorings without further queries: disease-anchored 1.5 %, biomarker-anchored ≈ 7 %, pooled 34 % but below chance. That is the whole map. Going further would be hunting for a friendlier harness — the one move this project's ethic forbids.

---

## 10. Related work and positioning

This is a **selective-prediction / cautionary** contribution, and it is deliberately **orthogonal** to the predictor literature rather than competitive with it.

The target-level enrichment of genetic support is established (Nelson et al. 2015; King et al. 2019; Minikel et al. 2024). The direction-of-effect problem has a current state of the art that this work corroborates rather than challenges: Chen et al. (2025) achieve macro-AUROC 0.59 for gene–disease-specific DOE across 47,822 pairs (vs. 0.85 disease-agnostic), with performance rising where genetic evidence is denser. Our runs land squarely in that regime from a falsification-first direction — indication-anchored coverage ≈ 1.5 % with near-total abstention, pooled accuracy below its own baseline. The 0.59 ceiling and the independently reported inhibitor-majority class imbalance do not threaten this story; they **validate** it. A direction-matched score already exists (the genetic priority score and its direction-of-effect variant, Duffy et al. 2024; Chen et al. 2024), and the directional fields we consume come from Open Targets' direction-of-effect assessment over eight evidence sources (Buniello et al. 2025).

The unoccupied niche is the auditor's, not the predictor's. The Do-lab products answer "what is the predicted direction for all 47,822 pairs?"; this work answers "which of those can a human trust, and which is 0.59 noise wearing a number?" A tool that **refuses where the state of the art is at chance** is the complement their own 0.59 implies. The two genuinely new, methods-level contributions are: (i) a named and quantified failure mode — naive trait-pooling of platform directional evidence is below a majority baseline and one-directional; and (ii) an interpretable mechanism for it — the trait-relative loss-of-function-disease confound, attributed three independent ways and named down to specific genes — packaged as an abstention-first design with a hard circularity guard, where the coverage map *is* the deliverable.

This is **not** a competitive DOE predictor and should not be positioned as one. It commits on two unselected targets and one is wrong. What it is — a reproducible boundary plus a cautionary harness — is genuinely under-supplied: the field publishes positive predictors; very few quantify the naive-pooling inversion or argue, with receipts, that silence is the correct output where aligned evidence is absent.

---

## 11. Limitations

- **n is pairs, not independent genes.** 132 pairs over 78 targets; committed n is tiny by construction. The exact-indication accuracy is unmeasurable, not merely low, and we say so.
- **Ground truth is mostly winners.** `clinical_precedence` reflects approved/clinical mechanisms; failed programs are under-represented, which is why coverage feels thin. The principled next step is re-anchoring to clinical outcomes *including failures* (public approximation: ChEMBL phase + Open Targets stopped/withdrawn + ClinicalTrials.gov), as the foundational enrichment studies do via proprietary Pharmaprojects.
- **Eleven diseases, common-disease-weighted.** The bootstrap shows the headline numbers are robust to which of these 11 are chosen, but the panel is not exhaustive and is biased toward cardiometabolic/immune indications.
- **No molQTL / Mendelian-randomization tier in this benchmark.** The caller here is the Tier-1 burden/curated layer; a functional-biomarker cis-MR tier (weighted on pQTL, never eQTL-abundance) is the documented next commit-tier and is not exercised in these numbers.
- **Indication-anchoring reduces but does not eliminate the inversion** (IL17RA).
- **Public-data snapshot.** Results are tied to the Open Targets release queried; numbers should be regenerated against a pinned release for the camera-ready version.

---

## 12. Conclusion

Human genetic support tells you a drug *target* is worth pursuing; it rarely tells you which *way* to push it. On unselected approved-drug targets, a caller built on rare-variant burden and curated Mendelian direction commits on 1.5 % at the disease indication and abstains on the rest, because the directional signal there is essentially absent. Where the genetic and therapeutic arms align at the trait — biomarker-proximal targets, ≈ 5–7 % of the set — direction is recoverable (8/8 on a curated panel). Pooling genetics across a target's associated traits raises coverage to 34 % only by dredging loss-of-function rare-disease evidence that opposes the inhibitor drug, scoring 33 % — below the 73 % always-inhibitor baseline, with errors running one way.

The contribution is a boundary, not a predictor: **genetics speaks to drug direction only at biomarker-proximal targets; a tool that abstains elsewhere is correct, and one that pools to manufacture coverage is below chance.** That reproduces the state of the art (0.59 disease-specific DOE) from the opposite direction, and it delivers the refusal-first thesis with data — while stating plainly that the positive utility on unselected targets is, for now, near zero.

---

## References

1. Nelson MR, et al. The support of human genetic evidence for approved drug indications. *Nature Genetics* 47, 856–860 (2015).
2. King EA, Davis JW, Degner JF. Are drug targets with genetic support twice as likely to be approved? Revised estimates of the impact of genetic support for drug mechanisms on the probability of drug approval. *PLoS Genetics* 15, e1008489 (2019).
3. Minikel EV, Painter JL, Dong CC, Nelson MR. Refining the impact of genetic evidence on clinical success. *Nature* 629, 624–629 (2024). doi:10.1038/s41586-024-07316-0.
4. Duffy Á, et al. Development of a human genetics-guided priority score for 19,365 genes and 399 drug indications. *Nature Genetics* 56, 51–59 (2024).
5. Chen R, et al. Expanding drug targets for 112 chronic diseases using a machine learning-assisted genetic priority score. *Nature Communications* 15, 8891 (2024).
6. Chen R, Duffy Á, Park JK, et al. (Do R, senior author). Genetic evidence informs the direction of therapeutic modulation in drug development. *npj Drug Discovery* 2, 24 (2025). doi:10.1038/s44386-025-00027-0.
7. Buniello A, et al. Open Targets Platform: facilitating therapeutic hypotheses building in drug discovery. *Nucleic Acids Research* 53, D1467–D1475 (2025).
8. Open Targets Platform Documentation — Target–disease evidence: gene-burden direction-of-effect assessment (assumption of LoF; β/OR → risk/protective). platform-docs.opentargets.org/evidence.
9. Open Targets blog — How the Open Targets Platform integrates gene burden analyses (direction consistency 98 %/87 %; LoF enrichment; PCSK9 example). blog.opentargets.org.
10. Minikel EV, et al. Evaluating drug targets through human loss-of-function genetic variation. *Nature* 581, 459–464 (2020).

*Reference details to be finalized against the originals at submission.*

---

*Reproducibility: every result maps to a version-pinned notebook cell — the genetic-mimicry mapping (`desired_from_label`, unit-tested and audited against the Open Targets convention); the runtime circularity guard (`assert GT_SOURCE.isdisjoint(NONCIRC)`); the pooling-vs-anchoring contrast (pooled cell auto-discovers `associatedDiseases`, indication cell anchors to one EFO); the scope×rule grid; the mapping audit (A–D); and the disease-cluster bootstrap. Saved artifacts: `scope_rule_grid.json`, `mapping_audit.json`, `panel_bootstrap.json`, `indication_anchored.json`. Public data only (Open Targets, ChEMBL, ClinVar); pin the Open Targets release for the camera-ready run.*
