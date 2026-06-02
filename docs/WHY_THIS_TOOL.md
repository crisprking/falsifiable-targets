# Why this tool exists

A short note on what falsifiable-targets does that other target-validation
checklists don't, and why that matters if you're using it (or thinking about
contributing).

---

## The problem

Published target claims fail in predictable ways:

- **Phantom evidence.** A paper says "we inhibit HDAC4" but the molecule
  binds SAT instead. The mechanism is real; the named target is wrong.
- **Pan-family collapse.** A "selective" inhibitor turns out to hit
  every paralog in the family at comparable potency.
- **Retraction blindness.** A foundational result has been retracted but
  the conclusion still circulates because subsequent work cites the
  pre-retraction paper.
- **Selectivity gap.** A novel target's selectivity argument rests on
  orthology distance from a human paralog, but the orthology DBs
  disagree on whether the human paralog is even an ortholog.
- **Cherry-picked replication.** Three labs tried; one succeeded; the
  successful result becomes the citation.

These failure modes don't get caught by "is the paper in a good journal"
or "does the author have a track record." They get caught by asking
specific structured questions of public data — and asking them in the
*same order, the same way, every time*.

That's what this tool does. Seven rules. One verdict per claim. Cheapest
falsification experiment named. Everything content-addressed by SHA-256 so
two audits with matching hashes must agree byte-for-byte.

---

## What's different from a checklist

The framework looks superficially like a checklist. But three properties
make it qualitatively different:

### 1. SHA-locked ruleset

The seven rules — their applicability conditions, their evaluation logic,
their thresholds, their descriptions — are hashed into a single
SHA-256 digest. Two audits with matching ruleset SHA are guaranteed to
have been run under exactly the same rule definitions.

Why this matters: when someone publishes a falsifiable-targets audit
that says "TYK2 / psoriasis: SURVIVED", you can verify they ran it under
the ruleset version they claim by computing the hash of their checked-in
`smoke_test.py`. If the SHA matches, the rules they applied are
identical to the ones the framework's other audits used. If it doesn't
match, the audit is from a fork and any verdict comparison is invalid.

Checklists don't have this property. "Did you check selectivity?" can
mean fifteen different things across fifteen different reviewers. R7
means exactly one thing under v1.1.0 of the ruleset, and exactly a
*different* one thing under v1.2.0.

### 2. Sentinel-calibrated

The framework ships with 11 sentinel test cases — known-positive and
known-negative claims drawn from real history:

- PCSK9 / LDL (SURVIVED — should obviously survive)
- BCR-ABL / CML (SURVIVED — gold standard)
- TNF / RA (SURVIVED)
- SAT/HDAC4 phantom (FALSIFIED — should obviously die)
- STAP cells (FALSIFIED — Obokata retraction)
- Arsenic life GFAJ-1 (FALSIFIED — wood fume / Erb 2012)
- CETP / LDL clinical-genetics gap (FALSIFIED_WITH_CAVEATS)
- TEX10 orthology threshold case (FALSIFIED_WITH_CAVEATS)
- Cruzain peptidomimetic selectivity gap (FALSIFIED_WITH_CAVEATS)
- Novel-target selectivity gap (FALSIFIED_WITH_CAVEATS)
- Kinase class-collapse on a validated mechanism (FALSIFIED_WITH_CAVEATS)

Every release must produce the right verdict on every sentinel before
it ships. The sentinel suite is the framework's *unit tests* —
calibration that mistakes can't slip past.

A new rule that breaks any sentinel must either be fixed to preserve
existing sentinels OR justified with an updated sentinel + RFC. There
is no third path.

### 3. Falsification-cheapest, not evidence-strongest

Most validation frameworks ask "what's the strongest evidence for this
target?" That's the wrong question. The right question is "what's the
cheapest experiment that could prove this target wrong?"

Every rule, when it flags a claim, names a specific experiment at a
specific cost tier:

- `public_data_lookup` — a database query (free, minutes)
- `cheap_in_silico` — a docking calculation, a sequence alignment
  (cents, hours)
- `targeted_assay` — a single biochemical experiment (\$100s, days)
- `cohort_study` — a clinical cohort analysis (\$10Ks, months)
- `clinical_trial` — what you're trying to avoid (\$10Ms+, years)

The aggregator picks the *cheapest* available falsification. This
biases the framework toward catching errors before they get expensive,
which is the only kind of catching that matters in drug discovery
economics.

---

## What this is not

- **Not a substitute for human judgment.** The verdict is a structured
  output, not a recommendation. `FALSIFIED_WITH_CAVEATS` doesn't mean
  "abandon this target"; it means "the framework finds these specific
  open questions; resolve them or accept their risk."
- **Not a complete failure-mode catalogue.** Seven rules cover the
  failure modes the framework's authors and contributors have
  encoded. Many other failure modes exist. PRs welcome.
- **Not a way to launder questionable claims.** Producing a SURVIVED
  verdict requires the inputs to support SURVIVED. Massaging fixture
  values to engineer a SURVIVED is technically possible but the
  audit report stamps the claim SHA — anyone with the YAML can
  recompute the verdict, and any mismatch is detectable.

---

## Three things this framework has done that prove it works

1. **Demoted its creator's own headline.** The Ipi1 / Madurella audit
   (the in-house novel antifungal claim) returned
   `FALSIFIED_WITH_CAVEATS` with 2 substantive caveats under v1.1.0.
   The framework's first job was to audit the project that built it,
   and it didn't blink. Most validation tools written by a research
   group never produce a verdict the group dislikes about the group's
   own work.

2. **Survived a prediction it should have followed.** Before running
   the v1.3.0 TYK2 re-audit on live ChEMBL data, the framework's
   author predicted R6 would fire a substantive caveat (JAK family
   class-collapse). The data said otherwise — JAK1/2/3 pools were
   1.55× TYK2's pool, below the 2.0× threshold. The threshold stayed
   at 2.0. The verdict stayed `SURVIVED`. Lowering the threshold to
   force the predicted verdict would have been motivated reasoning;
   the sentinel calibration prevents it.

3. **Surfaced its own scope gaps in public.** The v1.2.x TYK2 audit
   ran under a ruleset where R6 didn't apply to validated-mechanism
   claims. The audit itself flagged this as a gap. v1.3.0 closed it.
   The framework's documentation names every known gap (see
   `docs/AUDIT_LIMITATIONS_v1_3.md`), in writing, before anyone
   complains.

---

## How to read an audit report

A falsifiable-targets audit produces a JSON object with these fields:

- `verdict` — one of SURVIVED, FALSIFIED_WITH_CAVEATS, FALSIFIED,
  INSUFFICIENT_DATA. **Always read this with the caveat list.**
- `ruleset_sha256` — which ruleset version was applied.
- `claim_sha256` — content-hash of the claim itself. Re-running
  with the same claim against the same ruleset must produce the
  same verdict.
- `per_rule` — list of seven rule results, each with `status`
  (passed / falsified / not_applicable / abstained), `confidence`,
  and `input_data` (the actual values the rule saw).
- `cheapest_falsification` — if anything was falsified, the cheapest
  named experiment.
- `substantive_caveats` — open scientific questions. Treat these as
  blockers.
- `operational_notes` — tooling/data-availability gaps. Treat these
  as "the framework couldn't fully evaluate this rule" — not as
  scientific findings.

The byte-level reproducibility property means: if you don't trust an
audit's verdict, you can re-run it from the JSON. Check out the repo
at the ruleset SHA, install with `pip install -e .`, run
`ft-audit <claim.yaml>`, and verify the verdict matches.

---

## Citation

See [`CITATION.cff`](CITATION.cff) for the BibTeX-equivalent metadata,
or just point at the GitHub URL.

If you publish an audit produced by this framework, please cite the
specific ruleset SHA you used. That way readers can recompute under
later rulesets and see how the verdict drifts (or doesn't) as the
framework matures.
