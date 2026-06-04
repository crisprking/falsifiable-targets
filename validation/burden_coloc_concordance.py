"""
burden_coloc_concordance.py — orthogonal validation of the recovery engine.

QUESTION. Where a gene carries BOTH, on the SAME Open Targets trait:
  * a rare-variant gene_burden direction   (coding / loss-of-function -> FUNCTION), and
  * a coloc-recovered direction            (cis eQTL/pQTL -> ABUNDANCE-inclusive),
do the two agree?

  AGREE     -> coloc-recovery matched the functional gold standard at that locus.
  DISAGREE  -> abundance != function decoupling (the IL6R soluble-decoy / TYK2 /
               GCK post-translational class). A DISAGREE is to be CLASSIFIED
               (real regulatory decoupling vs spurious trait match), never assumed
               to be engine error -- the engine recovers abundance honestly.

This is the axis ORTHOGONAL to the approved-drug benchmark: the benchmark scores
recovery against drug mechanism; this scores recovery against rare-variant function.

DESIGN (matches the validated pipeline; nothing reinvented):
  * COLOC direction is da.recover(ens, efo)["recovered"] VERBATIM -- the same path,
    same thresholds, and same FIX A/C/D/E hardening the benchmark and gates use.
    A transient API failure RAISES (OpenTargetsError) and is logged-and-skipped per
    gene, so a network blip can never silently zero a gene and inflate the rate.
  * BURDEN direction is harmonised with the engine's OWN da.desired_from_label().
  * Traits are DISCOVERED from each gene's real associatedDiseases (per-datasource
    scores), then burden vs coloc are matched on the IDENTICAL trait id -- no EFO
    guessing and no sign-mapping (both signals sit on the same trait axis). These
    discovery/burden query shapes are the ones that self-tested clean on PCSK9.
  * Reported at BOTH the pair level (gene x trait) AND the gene level. The GENE LEVEL
    is the defensible headline unit: Open Targets nests overlapping trait labels
    (e.g. PCSK9 appears under 4 lipid synonyms), so pair counts over-state n. Wilson
    95% CIs on every rate; with tiny n the interval, not the point estimate, is the
    honest figure.
  * Pins the Open Targets data version and the engine sha; writes a content-addressed
    manifest in the same provenance style as the other validation/results outputs.

RUN as ONE Kaggle cell, internet ON (~15-25 min: it walks each gene's real
associations rather than one guessed trait). "Restart kernel, run this first" is
safest, but it self-heals the engine import (and refuses a test-stubbed module)
regardless of kernel state.
"""
import sys, os, json, time, hashlib, math, glob, re, shutil, subprocess
from collections import defaultdict, Counter

# ---------------------------------------------------------------------------
# 1) acquire the HARDENED engine. Ignore stale /kaggle/input snapshots (they lack
#    cis_qtl_dirs); clone the repo if no good local copy exists. This is the exact
#    guard the benchmark/gates use after a stale dataset copy once shadowed the import.
# ---------------------------------------------------------------------------
REPO    = "https://github.com/crisprking/falsifiable-targets.git"
_FTREPO = "/kaggle/working/ftrepo"

def _ensure_hardened_engine():
    seen = []
    cands = [os.getcwd(), _FTREPO] + [os.path.dirname(p)
             for p in glob.glob("/kaggle/working/**/direction_audit.py", recursive=True)]
    for c in cands:
        f = os.path.join(c, "direction_audit.py")
        if c in seen or not os.path.exists(f):
            continue
        seen.append(c)
        src = open(f, encoding="utf-8", errors="ignore").read()
        if "def cis_qtl_dirs" in src and "class OpenTargetsError" in src:   # hardened markers
            return c
    if os.path.isdir(_FTREPO):
        shutil.rmtree(_FTREPO)
    print("  no hardened engine found locally -> cloning the repo")
    subprocess.run(["git", "clone", "--depth", "1", REPO, _FTREPO],
                   check=True, capture_output=True, text=True)
    return _FTREPO

_d = _ensure_hardened_engine()
if _d not in sys.path:
    sys.path.insert(0, _d)
sys.modules.pop("direction_audit", None)            # drop any test-stubbed copy; import fresh
import direction_audit as da
from direction_audit import post, recover, desired_from_label, OpenTargetsError
if getattr(da.gwas_locus_ids, "__name__", "") == "<lambda>":
    raise RuntimeError("direction_audit is monkeypatched by a test cell (gwas_locus_ids is a stub). "
                       "Run > Restart kernel, then run THIS cell FIRST.")
H4_MIN, MAG, FRAC, MIN_LOCI = da.H4_MIN, da.MAG, da.FRAC, da.MIN_LOCI
ENGINE_SHA = hashlib.sha256(open(da.__file__, "rb").read()).hexdigest()
print(f"  engine: {da.__file__}  (sha {ENGINE_SHA[:12]}; H4>={H4_MIN} |s|>={MAG} frac>={FRAC} min_loci={MIN_LOCI})")

# ---------------------------------------------------------------------------
# 2) helpers
# ---------------------------------------------------------------------------
def wilson(k, n, z=1.96):
    if not n:
        return (None, None)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = (z / d) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (round(max(0.0, c - h), 3), round(min(1.0, c + h), 3))

def ci_str(k, n):
    lo, hi = wilson(k, n)
    base = f"{k}/{n} = {(100 * k / n):.0f}%" if n else "n/a"
    return base + (f"  [95% CI {lo:.0%}-{hi:.0%}]" if n else "")

def ot_version():
    try:
        d = post("query{ meta{ name dataVersion{year month iteration} apiVersion{x y z} } }", label="meta")
        return ((d or {}).get("data") or {}).get("meta") or {}
    except OpenTargetsError as e:
        return {"error": str(e)[:160]}

def resolve(sym):
    d = post('query($q:String!){ search(queryString:$q, entityNames:["target"]){ hits{ id name entity } } }',
             {"q": sym}, "search")
    hits = ((((d or {}).get("data") or {}).get("search") or {}).get("hits")) or []
    exact = [h for h in hits if h.get("entity") == "target" and (h.get("name") or "").upper() == sym.upper()]
    pool = exact or [h for h in hits if h.get("entity") == "target"]
    return pool[0]["id"] if pool else None

def discover(ens):
    """Diseases this gene is associated with, split by which datasource scored them (>0)."""
    q = """query($e:String!){ target(ensemblId:$e){ approvedSymbol
      associatedDiseases(page:{index:0,size:200}){
        rows{ disease{ id name } datasourceScores{ id score } } } } }"""
    d = post(q, {"e": ens}, "assoc")
    tgt = (((d or {}).get("data") or {}).get("target") or {})
    rows = ((tgt.get("associatedDiseases") or {}).get("rows")) or []
    burden, gwas = {}, {}
    for r in rows:
        dis = r.get("disease") or {}; did = dis.get("id")
        for ds in (r.get("datasourceScores") or []):
            if did and (ds.get("score") or 0) > 0:
                if ds.get("id") == "gene_burden":           burden[did] = dis.get("name")
                elif ds.get("id") == "gwas_credible_sets":  gwas[did]   = dis.get("name")
    return {"sym": tgt.get("approvedSymbol"), "burden": burden, "gwas": gwas}

def burden_dir(ens, did):
    """Functional direction from rare-variant burden, harmonised with the engine's desired_from_label."""
    q = """query($e:String!,$f:String!){ target(ensemblId:$e){
      evidences(efoIds:[$f], datasourceIds:["gene_burden"], size:200){
        rows{ directionOnTarget directionOnTrait } } } }"""
    d = post(q, {"e": ens, "f": did}, "burden")
    rows = (((((d or {}).get("data") or {}).get("target") or {}).get("evidences") or {}).get("rows")) or []
    inh = act = 0
    for r in rows:
        e = desired_from_label(r.get("directionOnTarget"), r.get("directionOnTrait"))
        if e == "inhibitor": inh += 1
        elif e == "activator": act += 1
    if inh == act:
        return None, inh + act          # no net functional direction (0/0 or a tie)
    return ("inhibitor" if inh > act else "activator"), inh + act

def gene_level(results):
    """Collapse pair rows to one call per gene (majority burden vs majority coloc).
    Handles ontology duplication: PCSK9's 4 inhibitor/inhibitor labels -> 1 AGREE;
    GCK's 2 activator/inhibitor -> 1 DISAGREE. Flags within-gene inconsistency."""
    by = defaultdict(list)
    for r in results:
        by[r["gene"]].append(r)
    out = []
    for g, rs in sorted(by.items()):
        bmaj = Counter(r["burden"] for r in rs).most_common(1)[0][0]
        cmaj = Counter(r["coloc"]  for r in rs).most_common(1)[0][0]
        mixed = len({(r["burden"], r["coloc"]) for r in rs}) > 1
        out.append({"gene": g, "burden": bmaj, "coloc": cmaj, "agree": bmaj == cmaj,
                    "n_traits": len(rs), "mixed_within_gene": mixed})
    return out

# ---------------------------------------------------------------------------
# 3) panel — cardiometabolic + immune + neuro genes likely to carry burden evidence.
#    Extensible; the scan keeps only genes that actually have BOTH signals.
# ---------------------------------------------------------------------------
GENES = ["PCSK9","LDLR","APOB","NPC1L1","ABCG5","ABCG8","APOE","ANGPTL3","ANGPTL4","APOC3","LPL",
         "APOA5","CETP","LPA","SORT1","PNPLA3","GCKR","SLC30A8","TM6SF2","GCK","HNF1A","HNF4A",
         "KCNJ11","ABCC8","MC4R","GIPR","IRS1","WFS1","IL6R","IL23R","IL12B","TYK2","TNFRSF1A",
         "PTPN22","IL4R","CARD9","NOD2","IFIH1","JAK2","SH2B3","TREM2","GRN","ASGR1","ANGPTL8"]

W = 116
print("=" * W)
print("BURDEN x COLOC CONCORDANCE  —  through the HARDENED engine (recover() == benchmark path)")
print("=" * W)
OTV = ot_version()
print(f"  OpenTargets: {OTV}")

# ---- self-test: discovery must work before spending 20 minutes ----
probe_ens = resolve("PCSK9")
probe = discover(probe_ens) if probe_ens else None
if not probe or not probe.get("burden"):
    raise SystemExit("!! discovery self-test failed on PCSK9 (no burden traits). "
                     "Paste this line and I'll fix the associatedDiseases query — no long run wasted.")
print(f"  self-test PCSK9: burden traits={len(probe['burden'])}  gwas traits={len(probe['gwas'])}  "
      f"shared={len(set(probe['burden']) & set(probe['gwas']))}")

# ---------------------------------------------------------------------------
# 4) scan — per gene, per shared trait: burden(function) vs recover()(abundance)
# ---------------------------------------------------------------------------
print("\nscanning (burden vs coloc on identical traits)...")
results, failed = [], []
for sym in GENES:
    try:
        ens = resolve(sym)
        if not ens:
            continue
        info = discover(ens)
        shared = sorted(set(info["burden"]) & set(info["gwas"]))   # both a burden score AND gwas-coloc potential
        for did in shared:
            bdir, nb = burden_dir(ens, did)
            if not bdir:
                continue
            rec = recover(ens, did)                 # <-- hardened engine; abundance-inclusive coloc direction
            cdir = rec["recovered"]
            if cdir is None:                        # coloc refused (insufficient/conflicted) -> nothing to compare
                continue
            tname = info["burden"].get(did) or did
            row = {"gene": sym, "ens": ens, "trait_id": did, "trait": tname,
                   "burden": bdir, "coloc": cdir, "agree": bdir == cdir,
                   "consensus": rec["consensus"], "n_loci": rec["n_loci_decided"], "n_burden": nb}
            results.append(row)
            tag = "AGREE" if row["agree"] else "DISAGREE  <-- decoupled"
            print(f"  {sym:<9} {tname[:34]:<35} burden={bdir:<10} coloc={cdir:<10} "
                  f"{int(rec['consensus']*100):>3}%/{rec['n_loci_decided']:<3}  {tag}")
    except OpenTargetsError as e:
        failed.append({"gene": sym, "error": str(e)[:160]})
        print(f"  !! {sym:<9} QUERY_FAILED (excluded): {str(e)[:110]}")
        continue

# ---------------------------------------------------------------------------
# 5) headline — pair level AND gene level, with Wilson CIs
# ---------------------------------------------------------------------------
pair_agree = [r for r in results if r["agree"]]
pair_diss  = [r for r in results if not r["agree"]]
gene_rows  = gene_level(results)
gene_agree = [g for g in gene_rows if g["agree"]]
gene_diss  = [g for g in gene_rows if not g["agree"]]

print("\n" + "=" * W)
print("HEADLINE  (coloc-recovery vs the rare-variant functional gold standard)")
print("=" * W)
print(f"  genes attempted                     : {len(GENES)}"
      + (f"   | {len(failed)} QUERY_FAILED (excluded): {[f['gene'] for f in failed]}" if failed else ""))
print(f"  PAIR level  (gene x trait, over-counts ontology):")
print(f"     pairs with BOTH signals          : {len(results)}  across {len({r['gene'] for r in results})} genes")
print(f"     AGREE  coloc == burden(function) : {ci_str(len(pair_agree), len(results))}")
print(f"     DISAGREE abundance != function   : {ci_str(len(pair_diss),  len(results))}")
print(f"  GENE level  (one call per gene -- the defensible unit):")
print(f"     genes with BOTH signals          : {len(gene_rows)}")
print(f"     AGREE                            : {ci_str(len(gene_agree), len(gene_rows))}")
print(f"     DISAGREE (decoupling rate)       : {ci_str(len(gene_diss),  len(gene_rows))}")

print(f"\n  decoupled genes (audit -- classify: regulatory decoupling vs spurious trait):")
for g in gene_diss:
    egs = [r for r in results if r["gene"] == g["gene"] and not r["agree"]]
    traits = "; ".join(sorted({r["trait"][:30] for r in egs}))
    flag = "  [MIXED within gene]" if g["mixed_within_gene"] else ""
    print(f"     {g['gene']:<8} burden={g['burden']:<10} coloc={g['coloc']:<10} on: {traits}{flag}")
if not gene_diss:
    print("     (none)")

# ---------------------------------------------------------------------------
# 6) manifest — provenance-matched, content-addressed
# ---------------------------------------------------------------------------
manifest = {
    "scan": "burden_x_coloc_concordance_hardened",
    "engine_sha": ENGINE_SHA,
    "opentargets_version": OTV,
    "n_genes_attempted": len(GENES),
    "n_query_failed": len(failed),
    "query_failed": failed,
    "pair_level": {
        "n_pairs": len(results), "n_genes": len({r["gene"] for r in results}),
        "agree": len(pair_agree), "agreement_rate": (len(pair_agree) / len(results) if results else None),
        "agreement_ci": wilson(len(pair_agree), len(results)),
        "decoupling": len(pair_diss),  "decoupling_rate": (len(pair_diss) / len(results) if results else None),
    },
    "gene_level": {
        "n_genes": len(gene_rows),
        "agree": len(gene_agree), "agreement_rate": (len(gene_agree) / len(gene_rows) if gene_rows else None),
        "agreement_ci": wilson(len(gene_agree), len(gene_rows)),
        "decoupling": len(gene_diss),  "decoupling_rate": (len(gene_diss) / len(gene_rows) if gene_rows else None),
        "decoupling_ci": wilson(len(gene_diss), len(gene_rows)),
        "decoupled_genes": gene_diss,
        "genes": gene_rows,
    },
    "pairs": results,
}
manifest["manifest_sha"] = hashlib.sha256(
    json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

out_dir = "/kaggle/working/validation/results" if os.path.isdir("/kaggle/working") else "validation/results"
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "burden_coloc_concordance.json")
with open(out_path, "w") as fh:
    json.dump(manifest, fh, indent=2)
print(f"\n  wrote {out_path}   manifest sha {manifest['manifest_sha'][:12]}")
print("\nREAD: AGREE = coloc-recovery tracked rare-variant function at that locus; DISAGREE = abundance!=function,")
print("to be classified (e.g. GCK's post-translational regulation is a REAL decoupling; a poorly-matched trait is")
print("spurious). The GENE-LEVEL agreement (with its Wilson CI) is the figure for the writeup, not the pair count.")
