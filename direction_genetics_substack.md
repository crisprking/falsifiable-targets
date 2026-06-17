# Genetics tells you the target. It won't tell you which way to push.

### I built a tool to predict whether a drug should activate or inhibit its target. On unselected targets it stays silent 98.5% of the time — and the obvious fix makes it worse than a coin. Here's why that's the honest result, and what it says about a field that rewards confident answers.

---

There's a number drug hunters quote like scripture. Targets with human genetic support behind them succeed in the clinic at roughly **2.6 times** the rate of targets without it. It shows up in pitch decks and portfolio reviews; it's the reason "is there genetics behind this?" is now the first question in the room. The number is real. It has been replicated for a decade, it gets stronger when the causal gene is unambiguous, and it barely moves with effect size or how common the variant is.

But sit with what that number is actually about. It tells you a target **matters** — that perturbing this gene moves the disease. It says nothing about **which way to perturb it.** Should your drug turn the target up or turn it down? Agonist or antagonist? That second question is the entire ballgame. Get it backwards and the genetics was right, the biology was right, and you still spent a billion dollars building a drug that pushes the wrong direction.

So I went looking for how much human genetics can answer the second question on its own. The honest answer turned out to be less comfortable — and more interesting — than I expected.

## The hard half of the problem

Direction-of-effect is the hard half, and the field knows it. The current state of the art is a machine-learning framework from Ron Do's lab (Chen et al., *npj Drug Discovery*, 2025). Given a gene and its broad druggability, it predicts the right direction with an AUROC of 0.95 — excellent. Strip it down to *isolated* direction among druggable genes and it's 0.85 — still good. But narrow to the question a drug program actually faces — **the right direction for this gene at this specific disease** — and across 47,822 gene–disease pairs it scores **0.59.**

0.59 is barely off the floor. A coin is 0.50. Even a full embeddings model, trained for exactly this, lands a hair above chance on the disease-specific problem. That isn't a knock on their work; it's the shape of the data. Direction, at the disease level, is genuinely close to unpredictable from genetics alone.

Which left me a choice. I could try to push 0.59 higher — an arms race against a well-resourced ML team, and frankly the wrong job for a public auditor. Or I could ask the complementary question, the one nobody seemed to be answering: **on an unselected set of real drug targets, where does directional genetic signal exist at all — and when it's absent, is the correct move to admit it and stay silent?**

I built the tool to do the second thing. To refuse.

## A tool designed to shut up

The caller reads only the Open Targets evidence that carries a real direction label — rare-variant burden tests (which tell you whether a loss or gain of function raises or lowers disease risk) and curated clinical-variant sources like ClinVar. It deliberately ignores GWAS hits, because on the evidence row a GWAS association carries no direction at all — it tells you *something* is here, not which way it points.

Each directional row gets mapped to a therapeutic call by one rule, **genetic mimicry**: a good drug mimics the protective direction of nature's own experiment. If losing the gene's function protects you from disease, inhibit it. (That's PCSK9 — people born with broken PCSK9 have low LDL and healthy hearts, so we built PCSK9 *inhibitors*, and they work.) If losing function causes disease, you'd want to restore or activate. And so on, four cases.

The tool commits only when the evidence it's reading all points one way. If the rows conflict, it abstains. If there's no directional row at all, it abstains. And — this is the part that keeps the whole thing honest — the approved drug's known mechanism, the answer key, is locked out of the prediction by a guard that fails at runtime if it ever leaks in. The drug can never be used to predict the drug.

Then I ran it on 132 real approved-drug target–indication pairs across eleven common diseases — coronary disease, rheumatoid arthritis, psoriasis, COPD, osteoporosis, and the rest — and read the result straight.

## 98.5% silence

At the actual disease indication, the tool committed on **2 of 132 pairs.** It abstained on the other 130. Not because the evidence conflicted — almost none of those 130 had any directional row to conflict over. **The signal simply wasn't there.**

A 1.5% commit rate looks, at first glance, like a tool that doesn't work. I want to be precise about why it isn't. The tool was built to find directional genetic evidence anchored to a disease and call from it. The finding is that on unselected targets, *that evidence mostly doesn't exist at the disease endpoint.* The silence isn't the tool failing to answer; it's the tool correctly reporting that genetics, here, has nothing directional to say. A bootstrap over the eleven diseases (resampling them two thousand times) put a number even on the silence: in only 1,797 of 2,000 resamples does the tool commit on even a single pair. The accuracy on unselected targets isn't low — it's **unmeasurable.** And saying *unmeasurable*, with the resampling behind it, is more honest than dressing up an n-of-2 as a result.

That could have been the whole post: directional genetics is sparse, abstention is correct, here's the map. But then I did the obvious thing — the thing anyone reaching for coverage would do — and it broke in a way worth the whole rest of this piece.

## The trap: coverage bought below chance

If the signal is too sparse at the disease, why not pool it? Take every trait a target is associated with — not just the drug's indication — and gather all the directional genetics across all of them. More evidence, more coverage.

Coverage jumped from 1.5% to **34%.** Twenty-two times more targets answered. And the accuracy on those answers was **33%.**

Read that again. Not 33% above some baseline — 33% accuracy, full stop. The targets in this set are about 73% inhibitor, so a rock that says "inhibitor" every single time scores 73%. The pooled tool, with all that extra evidence, scored **forty points below the rock.** Pooling didn't just fail to help. It produced a tool that is *actively worse than guessing the majority class.*

And it failed one specific, diagnosable way. When the pooled caller committed, its recall was **9% on true inhibitors and 100% on true activators.** It is, to a first approximation, a machine that says **"activator"** about everything. It nails every target that really is an activator and gets nearly every inhibitor backwards. The faint statistical signal left over (a Matthews correlation of +0.16, for those who track such things) isn't discrimination — it's the fingerprint of a near-constant guess plus a few lucky hits. When I resampled the diseases, that residual's confidence interval ran from **+0.000** upward. It includes zero. There's essentially no real direction-telling in it.

## Why "activator," every time

The mechanism is the actual contribution here, and it's teachable in one sentence: **direction is a property of the trait you're reading, not of the target.**

When you pool a target's genetics across all its traits, what you're mostly pooling is *loss-of-function-causes-a-rare-disease* evidence — annotated at some Mendelian condition, not at the common disease your drug treats. And here's the trap: burden tests, the way Open Targets runs them, **assume every variant is a loss of function.** So any gene where rare variants raise the risk of *some* disease gets mechanically labeled "loss of function raises risk" → which the mimicry rule reads as "activate to fix it." It doesn't matter that the actual drug inhibits the target at a *different*, common disease. The rare-disease frame has already flipped the sign.

The cleanest example is a bone gene called **SOST.** Its only directional genetic evidence sits at *sclerosteosis*, a rare disorder of excessive bone growth. For sclerosteosis, the sign is correct. For osteoporosis — the common disease, the one the drug treats — it's exactly backwards. And the indication-anchored version of the tool, the one that refuses to pool, **correctly says nothing about SOST,** because at osteoporosis there's no directional row to mislead it. The same story repeats across the misses: ACE carries its signal at a rare kidney disorder; RANKL at recessive osteopetrosis; CETP at CETP deficiency; JAK2's drug exploits a gain-of-function mutation while the pooled evidence drags in the loss-of-function arm.

I checked three separate ways that this is the data's structure and not a bug in my code, because that's the difference between a finding and an embarrassment:

The mapping I used **is Open Targets' own documented convention** — it matches their logic four cases out of four, anchored on their own PCSK9 worked example. Every single one of the 80 burden rows in my set is labeled loss-of-function — **80 out of 80** — which is the assumption baked into their method, made visible. And critically, when I pointed the *same mapping* at a place where the genetics and the drug share a trait — LDL cholesterol — it recovered the truth unanimously: PCSK9, ANGPTL3, NPC1L1, and HMGCR all came back **inhibitor, zero conflict,** which is correct for all four. The mapping works perfectly where the frames align. It inverts everywhere they don't. That's not a coding error. That's the boundary.

## Where the signal actually lives

So genetics *can* tell you direction — just not where most people look for it. It works at **biomarkers.** When the genetic evidence and the therapeutic goal point at the same measurable thing — LDL for a lipid drug, a continuous trait a variant directly moves — the sign is trustworthy. A hand-curated panel of lipid and metabolic targets scored **8 out of 8**, and it discriminated *both* directions correctly (LDLR and GCK are activator targets in that set, not inhibitors, and it got them right). But that biomarker-proximal regime is only about **5 to 7%** of target space once you anchor to disease. It's a flashlight, not a floodlight. It illuminates a small, specific, genuinely useful region — and gives you nothing reliable outside it.

That's the map, complete, across all three ways of asking: **1.5%** of targets answerable at the disease, **~7%** answerable at biomarkers, **34%** answerable if you pool — but that 34% is below a coin. There's no fourth setting that buys you both coverage and correctness. I looked. Looking harder for one would have meant hunting for a friendlier way to score the same data until a number looked good, which is the one move a tool built on falsification isn't allowed to make.

## A field that rewards the wrong thing

Here's what I keep coming back to. Almost everyone in this corner of the field publishes **predictors** — models that output a number for every gene–disease pair, all 47,822 of them. That's valuable, and the 0.59 ceiling is honest work. But a predictor that's at chance on most of its inputs and confident-looking on all of them is a particular kind of dangerous, because the number *looks the same* whether it's the trustworthy 7% or the 0.59 noise. The score doesn't wear a label saying "this one's real."

The complementary tool — the one that's harder to publish because it brags about how often it shuts up — is the one that knows the difference. That refuses where the evidence isn't aligned, and tells you plainly that the silence is the answer. On an unselected target you hand this tool, the right output is, overwhelmingly, *"I don't know, and here's the structural reason genetics can't tell you."* That's not a failure mode. For most of target space, that's the only correct thing to say.

I'll be straight about what this is and isn't, because the whole point is the honesty. **This is not a better direction predictor.** It commits on two unselected targets and gets one of them wrong; the IL17RA pair is a real counterexample where the inversion happens *at* the indication, not just at remote rare diseases. The published 0.59 is the state of the art and this doesn't beat it. What this work does is **reproduce that 0.59 from the opposite direction** — falsification instead of prediction — and name a trap that anyone querying these databases directly will fall into: pool the genetics for coverage, and you get a confident "activator" machine that scores below a coin. The map of where the signal lives, the named and quantified failure, the mechanism traced down to specific genes, and a tool disciplined enough to refuse — that's the contribution. It's a boundary, and a boundary you can trust is worth more than a number you can't.

Genetics tells you the target is worth pursuing. On which way to push it, most of the time, the only honest thing it can say is nothing — and a tool that says nothing, loudly and for the right reason, is doing its job.

---

*The auditor, the full audit suite, and the reproducible numbers behind every figure here are public — including the bootstrap, the three-way mechanism attribution, and the circularity guard that locks the answer key out of the prediction. It demotes my own headline target before it will vouch for an FDA-approved one; the credibility is in what it refuses. Code and the full technical write-up: github.com/crisprking/falsifiable-targets.*

*Sources: Nelson et al., Nat Genet 2015; King, Davis & Degner, PLoS Genet 2019; Minikel et al., Nature 2024 (the 2.6× enrichment); Duffy et al., Nat Genet 2024; Chen, Duffy, Park et al. (Do, senior author), npj Drug Discovery 2025 (the 0.59 disease-specific ceiling and the inhibitor-majority imbalance); Buniello et al., Nucleic Acids Research 2025 (the Open Targets direction-of-effect assessment).*
