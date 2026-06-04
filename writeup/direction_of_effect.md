# Most genetic evidence can't tell you which way to push. `direction_audit` gets half of it back — and flags the one place it can't.

*A direction-of-effect audit for drug targets, built on public data, validated against approved drugs, and explicit about the one structural blind spot that shows up twice. Part of [falsifiable-targets](https://github.com/crisprking/falsifiable-targets).*

---

PCSK9 is one of the most validated drug targets in existence. Its inhibitor, evolocumab, is FDA-approved and lowers LDL by blocking the protein. So ask Open Targets for the genetic direction of effect of PCSK9 on coronary artery disease — should you inhibit it or activate it? — and it hands you back nothing. Not "inhibit." Not "unsure." The `directionOnTarget` field is null across all 38 GWAS evidence rows.

That is the gap. The field has been treating association as approval, and the bill comes due in the clinic, where pointing the right way at the wrong direction is the most expensive mistake a program can make.

## What exists, and what none of it does

Genetic priority scores are not new. Nelson and colleagues showed targets with human genetic support succeed roughly twice as often; Minikel and colleagues sharpened that to about a 2.6-fold enrichment; Duffy and colleagues built a continuous genetic priority score; Open Targets publishes its own harmonized direction-of-effect layer. They are useful and I lean on the same data.

What none of them do is **refuse**. They emit a number, or a direction, for every target — and where the underlying direction is actually missing, that output is resting on association alone. Which is precisely where the confident, wrong answers live. The contribution here is not another score. It is a tool that says *I don't know* exactly when it doesn't, and recovers a real direction exactly when it can.

## The labels exist only where direction was free

Open Targets can harmonize direction-of-effect only when the evidence is *intrinsically* directional — rare-variant burden tests (collapse loss-of-function alleles, see which way the trait moves) and curated clinical variants. Common-variant GWAS is not directional, and it is the overwhelming majority of genetic support. Pulled across a panel of well-studied targets, source by source:

| genetic source | rows | directional | % |
|---|---:|---:|---:|
| `gwas_credible_sets` | 1442 | 0 | 0% |
| `gene_burden` | 117 | 117 | 100% |
| `eva` (ClinVar) | 2 | 0 | 0% |

Every directional row came from burden testing. Not one came from GWAS. On the GWAS rows, `directionOnTarget`, `directionOnTrait`, `beta`, and `oddsRatio` are all null — the direction simply isn't carried. The consequence is worse than missing data: a tool consuming the labels cannot distinguish *no direction available* from *no association*. Both look like silence, and it can be blind on ~95% of its evidence while looking confident.

## The method: recover from colocalisation, or refuse

The direction a GWAS locus won't give you is often recoverable one layer down. Where a disease credible set **colocalises** with a *cis* eQTL or pQTL for the target, the two signals share a causal variant and the sign of the effect-size relationship (`betaRatioSignAverage`) becomes readable:

- **sign > 0** — the allele that raises the target's product raises risk → more is bad → **inhibit**
- **sign < 0** — the allele that raises the product lowers risk → more is good → **activate**

I don't trust one colocalisation; molecular-QTL data is noisy and *trans* signals impersonate *cis*. The rule is a **locus consensus**: each of a target's GWAS loci votes ±1 from its strong cis colocalisations (`|sign| ≥ 0.8`, `H4 ≥ 0.8`), and a direction is called only when **≥80% of ≥3 decided loci agree**. Thin or split evidence returns nothing, deliberately.

So for any target–disease claim the audit does one of three things: use a **burden-derived label** if present (the functional axis — what losing the gene does, which is what an inhibitor mimics; drug-derived labels are excluded as circular); else attempt **colocalisation recovery**; else **refuse** — `INSUFFICIENT_DIRECTION` or `RECOVERY_CONFLICTED`, and say which. Every verdict is content-addressed: the inputs hash to a sha256, so any claim below is checkable by anyone running the same query against public data.

I calibrated on two targets whose direction isn't in question. PCSK9 recovered `inhibitor` at 97% consensus over its *cis*-locus (matching evolocumab) — and as the locus-independence gate below shows, that 97% is agreement across 21 distinct GWAS studies, not 29 independent loci. LPL, where genetics and physiology both say *more is better*, recovered `activator` in both directions tested. Both ways, right, from colocalisation alone.

## Validation against approved drugs

A method that rescues PCSK9 isn't validated — I might have picked targets where it works. So the panel was committed *before* the run (hash `5ac82f6268a65c38`): 18 target–disease pairs whose mechanism is settled by an approved drug, both directions, cardiometabolic and immune, plus one deliberately contested control. Mechanism is cited from the drug. Every target is reported; nothing is dropped.

| Target | Disease | Approved drug | Mechanism | Evidence | Verdict | SHA |
|---|---|---|---|---|---|---|
| PCSK9 | CAD | evolocumab | inhibitor | coloc-recovered | **VOUCHES** | `eb6730582cee` |
| HMGCR | CAD | statins | inhibitor | thin coloc | refused (insufficient) | `0b4b359b15f2` |
| NPC1L1 | LDL | ezetimibe | inhibitor | burden-label | **VOUCHES** | `b42a33df4d7d` |
| ANGPTL3 | triglycerides | evinacumab | inhibitor | burden-label | **VOUCHES** | `85c32c4d3c75` |
| APOC3 | triglycerides | volanesorsen | inhibitor | burden-label | **VOUCHES** | `9fce8f07c60d` |
| LPL | triglycerides | (LPL activation) | activator | burden-label | **VOUCHES** | `24409ed52442` |
| CETP | CAD | anacetrapib | inhibitor | coloc-recovered | **VOUCHES** | `275524228dce` |
| PPARG | T2D | thiazolidinediones | activator | burden-label | **VOUCHES** | `d16c597bc405` |
| GLP1R | T2D | semaglutide | activator | coloc-recovered | **VOUCHES** | `0062a7410fa1` |
| DPP4 | T2D | sitagliptin | inhibitor | no coloc | refused (insufficient) | `88d8d0db2752` |
| SLC22A12 | gout | lesinurad | inhibitor | burden-label | **VOUCHES** | `59912f59c4b3` |
| IL6R | RA | tocilizumab | inhibitor | coloc-recovered | **CAVEAT — decoy** | `515834d00c3d` |
| TNF | RA | adalimumab | inhibitor | no coloc (MHC) | refused (insufficient) | `d6932ef7db72` |
| IL23R | IBD | anti-IL23 axis | inhibitor | coloc-recovered | **VOUCHES** | `83eaab5c000b` |
| IL12B | IBD | ustekinumab | inhibitor | coloc-recovered | **VOUCHES** | `c520001a2dfc` |
| IL4R | atopic eczema | dupilumab | inhibitor | no coloc | refused (insufficient) | `748b1acb5e25` |
| TYK2 | psoriasis | deucravacitinib | inhibitor | coloc-recovered | **WRONG — abundance decoupled** | `baff0f4b9717` |
| TSLP | asthma | tezepelumab | inhibitor | *invalid EFO* | excluded | `0aaf96e7d25e` |
| GIPR | T2D | (contested) | contested | coloc | refused (conflicted, 73%) | `b984b4efa127` |

The numbers, on the 12 GWAS-only targets — the claims a label-based tool is blind to:

- **Coverage** — recovery returned a direction for **6 of 12 = 50%** (Wilson 95% CI 25–75%). The rest refused, on absent or thin cis-colocalisation. (Drop TSLP, whose disease ontology is invalid, and it reads 6 of 11 ≈ 55%; I lead with the figure the committed manifest records.)
- **Accuracy** — of those six, **five matched the approved drug = 83%** (CI 44–97%). PCSK9, CETP, IL23R, IL12B → inhibitor; GLP1R → activator. One was wrong (TYK2).
- **Reliability by confidence tier** — the load-bearing number, not the headline accuracy. Each recovered call is graded by whether a *protein* QTL corroborates it: **HIGH** (pQTL agrees) was **3 of 3** correct, **STANDARD** (expression-QTL only) was **2 of 3**, and **CAVEAT** (pQTL conflicts) did not arise on this panel. **Zero calls were wrong among the pQTL-corroborated tier.** The single miss falls in the tier the engine pre-flags as weaker — TYK2, below.
- **Burden labels** matched the drug **6 of 6**. The **contested control (GIPR) refused**, at 73% consensus.

Two integrity gates underwrite these numbers. **Sign-orientation calibration** runs nine controls whose direction is textbook (six inhibitor, three activator); NPC1L1 abstains at 76%, below the consensus floor, and of the eight that do recover, **all eight point the correct way — zero reversals**. A flipped sign convention would have made PCSK9 read `activator`; it didn't, in either direction. **Locus independence** then asks whether "N loci" means N independent regions — it doesn't, and shouldn't: *cis*-colocalisation localises to the target's own locus, so each target is a single genomic region, and the consensus percentage is cross-study agreement *at* that locus. PCSK9's 97% is 21 distinct GWAS studies, TYK2's 100% is 14, IL23R's 93% is 14 — real replication across cohorts, reported as such rather than as an inflated independent-loci count.

## The tool publishes its own misses

The accuracy number isn't the integrity signal — anyone can report a number. The signal is what the audit did with its failure. **TYK2 is in the table marked WRONG, in the same pre-committed run, under its own SHA.** The framework did not quietly drop the target it got wrong or retrofit a flag to launder it; it recorded the miss next to the hits, in the direction the claim's structure predicts, where anyone can diff it against the code. A target audit that hides its own errors is a marketing asset. One that prints them is an instrument. It happened twice: TYK2 stands in the table as a published error, and — below — the contested control GIPR stopped the tool from overturning its own refusal. The framework's failures are load-bearing, not buried.

## Where it breaks — one blind spot, twice

Two targets failed. They are the same lesson from two angles.

**IL6R.** Tocilizumab inhibits it; recovery returned `activator`. The Asp358Ala variant raises *soluble* IL6R while *lowering* signaling and reducing RA and cardiovascular risk — it phenocopies the drug. The pQTL measures the soluble pool, which runs opposite to the signaling the drug acts on. Flagged in advance as a documented decoy → `RECOVERED_CAVEAT_DECOY`. And a diagnostic confirms the eQTL and pQTL *agree*, so no internal consistency check can catch it.

**TYK2.** Deucravacitinib inhibits it, and the protective genetics (the P1104A partial-loss-of-function allele) agree. Recovery returned `activator`, 100% consensus — which the locus-independence gate shows is 14 GWAS studies agreeing at one *cis*-locus, so the miss is robust, not a noisy fluke. It is reported as a genuine error, and it lands in the STANDARD (expression-QTL-only) tier: the engine never had protein corroboration for it, which is exactly why the tier system pre-flags such calls as weaker. The reason is the same family: P1104A acts on enzyme *activity*, a coding-function variant; the cis-eQTL coloc reads *expression abundance*; those axes come apart for TYK2 (which also sits in a gene-dense, high-LD stretch of 19p13 where colocalisation is fragile).

The conceptual core, pulled out of the weeds:

> **Colocalisation recovers the direction of the *abundance* → trait relationship. For most targets, abundance is function, and the recovered direction is the therapeutic one. It diverges exactly where abundance is decoupled from function** — soluble decoys whose measured pool is inverse to signaling, and targets whose causal genetics are coding-activity variants invisible to expression QTLs.

This is also why burden-derived labels are the better evidence where they exist: a burden test reads the *functional* axis, which is what a drug mimics. Recovery is a good proxy for that — about four times in five here — but a proxy, and it fails in a characterizable way, not a random one. Critically, the failure is not detectable from colocalisation data alone, so it is handled by curation and disclosure, not a cleverer filter.

And the flag is no longer merely curated — it is now *discoverable*. A burden test reads function; colocalisation reads abundance; so wherever a target carries both, a disagreement between them is a candidate abundance≠function decoupling. Run that comparison across 44 candidate genes and the two directions agree, at the gene level, for **6 of 8** that overlap (**75%**, CI 41–93%). The disagreements are not noise. **GCK** is the worked example: loss-of-function *GCK* variants cause MODY2 hyperglycemia, so the functional direction is *activator* — and glucokinase activators (dorzagliatin, approved in China) are a real drug class — so burden correctly says `activator`, while coloc reads `inhibitor`, because glucokinase activity is set post-translationally by GKRP sequestration, invisible to an expression QTL. The same failure as TYK2, surfaced automatically rather than by hand. (One further gene-level disagreement, IRS1, is flagged as a likely trait-mismatch artifact, not asserted as a real decoupling.) The curated flag has become a measurement.

## Coverage where the disease is mute: the biomarker bridge

Half of those refusals aren't a genetics problem — they're an *endpoint* problem. HMGCR refuses on coronary artery disease because the disease-level cis-colocalisation is thin, but it's unmissable on LDL cholesterol, where it has loci to spare. So when a disease endpoint comes back direction-mute, the audit re-runs the identical locus-consensus recovery on a causally-signed upstream biomarker and maps the sign through (higher LDL means higher coronary risk, so the map is identity).

A bridge is only trusted after it passes a **dual-endpoint concordance** test: for every target that resolves on *both* the disease and the biomarker, do the two directions agree? The LDL→CAD bridge passed **4 of 4**, and in both directions — PCSK9 and CETP inhibitor on each endpoint, LDLR and ABCG5 activator on each. That's the bridge being *faithful*, not just convenient. Folded into the benchmark it moves coverage from 50% to 58% and accuracy from 83% to 86%; one target moves (HMGCR, rescued via LDL), and every direct call, label, the IL6R flag, and the TYK2 error stay exactly as they were. Rescued-via-biomarker verdicts carry a weaker confidence grade and their own falsifier — an MR of the biomarker on the disease, or a disease-endpoint coloc, upgrades them to direct.

The rule that keeps this honest came from the contested control. When I extended the idea to a glycemic→type-2-diabetes bridge, it tried to "rescue" GIPR — a target whose disease-level genetics are genuinely conflicted — into a confident `inhibitor` off HbA1c. That's the failure mode a biomarker bridge invites: overriding a meaningful refusal with a cleaner-looking proxy. So the fallback fires *only* when the disease endpoint is direction-mute (`INSUFFICIENT_DIRECTION`), never when it is conflicted (`RECOVERY_CONFLICTED`). A conflict is information, not a gap. GIPR stays refused — which, for a target where agonists and antagonists both have a metabolic rationale, is the only honest answer.

## What I'd claim

A direction-of-effect audit can speak to the GWAS majority that label-based tools cannot — recovering a mechanism direction for roughly half of otherwise-silent claims directly, rising to 58% once a validated biomarker bridge fills the endpoint-sparse cases, at 83–86% accuracy against approved drugs, with zero wrong among its protein-corroborated calls — **if** it prefers burden-derived functional labels where they exist, refuses on thin or contested evidence instead of guessing, only trusts a biomarker bridge after it passes head-to-head concordance, and openly discloses its one structural blind spot. Not a black box that emits a confident mechanism for every target; a tool that rescues what it can defend, refuses what it can't, names where its own method is the wrong instrument, and hashes every verdict so you can check it.

## Before your next target

Before you commit to a mechanism direction from genetic evidence alone, find out *which kind* of evidence it is. Is the direction burden-derived (functional, trustworthy), recovered from colocalisation (the abundance proxy — good, with the decoupling caveat), or merely assumed from an association? If it's the third, the audit returns `INSUFFICIENT_DIRECTION` — and so should you. The tool is public, the panel is pre-committed, the hashes are checkable. Use it, or build something better and beat the number.

## Limits and what's next

This is a curated benchmark, not a population scan; turning "50% of these" into a platform-wide coverage statistic needs the bulk colocalisation layer and is the next build. The abundance≠function flag began curated (IL6R, then TYK2) but is now a discovery method: the burden×coloc concordance scan above turns disagreements into candidate decouplings automatically, and surfaced GCK without prompting. And recovery is hostage to molecular-QTL coverage; most refusals are simply where cis-colocalisation doesn't exist yet.

## Code and data

Everything is public and reproducible at **[github.com/crisprking/falsifiable-targets](https://github.com/crisprking/falsifiable-targets)**: the engine (`direction_audit.py`), the pre-specified benchmark (`validation/systematic_benchmark.py`), the committed run with every SHA in this article (`validation/results/systematic_run.json`), and the development trail — including the diagnostic that *falsified* an earlier auto-detection idea. Public data only (Open Targets GraphQL), no private inputs, no model in the loop. If a direction here is wrong, it is wrong reproducibly — paste the SHA, run the query, and check.

---

### References

1. Nelson MR, et al. The support of human genetic evidence for approved drug indications. *Nature Genetics*, 2015.
2. Minikel EV, et al. Refining the impact of genetic evidence on clinical success. *Nature*, 2024. (≈2.6-fold enrichment; DOI 10.1038/s41586-024-07316-0.)
3. Duffy Á, Forgetta V, et al. A genetics-guided priority score for drug-target pairs. (genetic priority score work)
4. Ochoa D, et al. The Open Targets Platform. *Nucleic Acids Research*. (platform + direction-of-effect documentation)
5. Sabatine MS, et al. Evolocumab and clinical outcomes in patients with cardiovascular disease (FOURIER). *NEJM*, 2017.
6. Armstrong AW, et al. Deucravacitinib in moderate-to-severe plaque psoriasis (POETYK PSO-1). *JAAD*, 2023.
7. Swerdlow DI, et al. The interleukin-6 receptor as a target for prevention of coronary heart disease: a Mendelian randomisation analysis. *Lancet*, 2012.
8. Ferreira RC, et al. Functional IL6R Asp358Ala variant and autoimmunity.
9. Dendrou CA, et al. Resolving TYK2 locus genotype-to-phenotype differences in autoimmunity. *Science Translational Medicine*, 2016.
10. falsifiable-targets — direction-of-effect audit and benchmark. github.com/crisprking/falsifiable-targets
