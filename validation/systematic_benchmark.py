"""
systematic_benchmark.py — pre-specified approved-drug benchmark for direction_audit.
Run as ONE Kaggle cell (internet ON; ~8-12 min). Self-healing: it force-reimports the
REAL engine, so a leftover test stub (the SYM/0-loci monkeypatch) can no longer poison it;
if a stub is somehow still active it refuses with a clear message instead of reporting 0%.

Reports every panel target (no cherry-picking), with Wilson 95% intervals, accuracy split
by confidence tier (pQTL-corroborated vs eQTL-only — where documented decouplings like TYK2
live), QUERY_FAILED rows excluded from denominators, and the Open Targets data version +
engine sha pinned in the manifest so the numbers are reproducible.
"""
import sys, os, json, hashlib, math, glob, time, shutil, subprocess

# ---- acquire the HARDENED engine, ignoring stale copies (e.g. an old /kaggle/input snapshot) ----
REPO = "https://github.com/crisprking/falsifiable-targets.git"
_FTREPO = "/kaggle/working/ftrepo"

def _ensure_hardened_engine():
    """Return a dir whose direction_audit.py actually IS the hardened engine (defines
    cis_qtl_dirs + OpenTargetsError). Prefers a local hardened copy; else clones the repo
    fresh. Deliberately skips stale copies such as an attached /kaggle/input dataset."""
    seen = []
    locals_ = [os.getcwd(), _FTREPO] + [os.path.dirname(p) for p in glob.glob("/kaggle/working/**/direction_audit.py", recursive=True)]
    for c in locals_:
        f = os.path.join(c, "direction_audit.py")
        if c in seen or not os.path.exists(f):
            continue
        seen.append(c)
        src = open(f, encoding="utf-8", errors="ignore").read()
        if "def cis_qtl_dirs" in src and "class OpenTargetsError" in src:
            return c
    if os.path.isdir(_FTREPO):
        shutil.rmtree(_FTREPO)
    print("  no hardened engine found locally -> cloning the repo")
    subprocess.run(["git", "clone", "--depth", "1", REPO, _FTREPO], check=True, capture_output=True, text=True)
    return _FTREPO

_d = _ensure_hardened_engine()
if _d not in sys.path:
    sys.path.insert(0, _d)
sys.modules.pop("direction_audit", None)   # CRUCIAL: drop any stale/stubbed copy, import the hardened one
import direction_audit as _da
from direction_audit import post, verdict, cis_qtl_dirs, confidence_tier, OpenTargetsError
if getattr(_da.gwas_locus_ids, "__name__", "") == "<lambda>":
    raise RuntimeError("direction_audit is monkeypatched by a test cell (gwas_locus_ids is a stub). "
                       "Run > Restart kernel, then run THIS cell first.")
print(f"  engine: {_da.__file__}")
ENGINE_SHA = hashlib.sha256(open(_da.__file__, "rb").read()).hexdigest()

# documented: measured abundance decoupled from function (soluble/decoy isoform). curated, not auto-detected.
DECOY = {"ENSG00000160712"}   # IL6R (sIL6R decoy; rs2228145)

# symbol, ensembl, efo, known_mechanism, drug (citation for the mechanism call)
PANEL = [
  ("PCSK9",   "ENSG00000169174","EFO_0001645",  "inhibitor","evolocumab"),
  ("HMGCR",   "ENSG00000113161","EFO_0001645",  "inhibitor","statins"),
  ("NPC1L1",  "ENSG00000015520","EFO_0004611",  "inhibitor","ezetimibe"),
  ("ANGPTL3", "ENSG00000132855","EFO_0004530",  "inhibitor","evinacumab"),
  ("APOC3",   "ENSG00000110245","EFO_0004530",  "inhibitor","volanesorsen"),
  ("LPL",     "ENSG00000175445","EFO_0004530",  "activator","LPL-pathway (GoF protective)"),
  ("CETP",    "ENSG00000087237","EFO_0001645",  "inhibitor","anacetrapib (borderline outcome)"),
  ("PPARG",   "ENSG00000132170","MONDO_0005148","activator","thiazolidinediones"),
  ("GLP1R",   "ENSG00000112164","MONDO_0005148","activator","semaglutide"),
  ("DPP4",    "ENSG00000197635","MONDO_0005148","inhibitor","sitagliptin"),
  ("SLC22A12","ENSG00000197891","EFO_0004274",  "inhibitor","lesinurad (URAT1)"),
  ("IL6R",    "ENSG00000160712","EFO_0000685",  "inhibitor","tocilizumab [decoy]"),
  ("TNF",     "ENSG00000232810","EFO_0000685",  "inhibitor","adalimumab"),
  ("IL23R",   "ENSG00000162594","EFO_0003767",  "inhibitor","via IL23 (ustekinumab/risankizumab)"),
  ("IL12B",   "ENSG00000113302","EFO_0003767",  "inhibitor","ustekinumab"),
  ("IL4R",    "ENSG00000077238","EFO_0000274",  "inhibitor","dupilumab"),
  ("TYK2",    "ENSG00000105397","EFO_0000676",  "inhibitor","deucravacitinib"),
  ("TSLP",    "ENSG00000145777","EFO_0000270",  "inhibitor","tezepelumab"),
  ("GIPR",    "ENSG00000010310","MONDO_0005148","contested","agonist-vs-antagonist debate"),
]
RESCUED = {"RECOVERED_CONCORDANT","RECOVERED_DISCORDANT"}
REFUSED = {"RECOVERY_CONFLICTED","INSUFFICIENT_DIRECTION"}
LABELED = {"LABEL_CONCORDANT","LABEL_DISCORDANT"}

def wilson(k, n, z=1.96):
    if not n: return (None, None)
    p = k / n; d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d
    h = (z/d) * math.sqrt(p*(1-p)/n + z*z/(4*n*n))
    return (round(max(0.0, c-h), 3), round(min(1.0, c+h), 3))

def line(k, n):
    if not n: return f"{k}/0 = n/a (no calls in this bucket)"
    lo, hi = wilson(k, n)
    return f"{k}/{n} = {(100*k/n):.0f}%  [95% CI {lo:.0%}-{hi:.0%}]"

def ot_version():
    for q in ("query{ meta{ dataVersion{ year month iteration } apiVersion{ x y z } } }",
              "query{ meta{ name dataVersion{ year month } } }"):
        try:
            d = post(q, {}, "meta"); m = ((d or {}).get("data") or {}).get("meta")
            if m: return m
        except Exception:
            continue
    return None

def efo_name(efo):
    d = post('query($f:String!){ disease(efoId:$f){ name } }', {"f": efo}, "dn")
    return (((d or {}).get("data") or {}).get("disease") or {}).get("name") or "?? EFO NOT FOUND"

def call_target(ens, efo, claim, attempts=2):
    last = None
    for _ in range(attempts):
        try:
            return verdict(ens, efo, claim), None
        except OpenTargetsError as e:
            last = str(e); time.sleep(2)
    return None, last

def tier_of(ens, efo, direction):
    """Confidence tier — extra locus pass, run ONLY for rescued targets to keep API load down."""
    try:
        return confidence_tier(direction, cis_qtl_dirs(ens, efo)["pqtl_dir"])
    except OpenTargetsError:
        return None

def run(panel):
    out = []
    print(f"{'target':<9}{'disease':<28}{'known':<10}{'verdict':<24}{'dir':<10}{'tier':<10}{'c%':>4}{'loc':>5}  flags")
    print("-" * 126)
    for sym, ens, efo, known, drug in panel:
        time.sleep(0.2)
        try:
            nm = efo_name(efo)
        except OpenTargetsError:
            nm = "?? (efo lookup failed)"
        claim = known if known in ("inhibitor","activator") else "inhibitor"
        v, err = call_target(ens, efo, claim)
        if v is None:
            out.append({"sym":sym,"ens":ens,"efo":efo,"efo_name":nm,"known":known,"drug":drug,
                        "verdict":"QUERY_FAILED","direction":None,"vouches":False,"tier":None,
                        "consensus":None,"n_loci":None,"source":"error","decoy":(ens in DECOY),
                        "sha":None,"error":err})
            print(f"{sym:<9}{nm[:27]:<28}{known:<10}{'QUERY_FAILED':<24}{'-':<10}{'-':<10}{'':>4}{'':>5}  !! {err[:38] if err else ''}")
            continue
        rec = v.get("recovery", {}) or {}
        resolved = rec.get("symbol")
        tier = tier_of(ens, efo, v.get("direction")) if v["verdict"] in RESCUED else None
        flags = []
        if ens in DECOY and v["verdict"] in RESCUED:
            v = dict(v); v["verdict"] = "RECOVERED_CAVEAT_DECOY"; v["vouches"] = False
            flags.append("DECOY")
        symflag = "" if (resolved is None or resolved == sym) else f"!{resolved}"
        efoflag = "EFO?" if nm == "?? EFO NOT FOUND" else ""
        c = int(round(rec.get("consensus",0)*100)) if rec else ""
        nl = rec.get("n_loci_decided","") if rec else ""
        out.append({"sym":sym,"ens":ens,"efo":efo,"efo_name":nm,"known":known,"drug":drug,
                    "verdict":v["verdict"],"direction":v.get("direction"),"vouches":v["vouches"],
                    "tier":tier,"consensus":rec.get("consensus"),"n_loci":rec.get("n_loci_decided"),
                    "source":v["source"],"decoy":(ens in DECOY),"sha":v["sha256"]})
        print(f"{sym:<9}{nm[:27]:<28}{known:<10}{v['verdict']:<24}{str(v.get('direction')):<10}"
              f"{str(tier):<10}{str(c):>4}{str(nl):>5}  {symflag} {efoflag} {'/'.join(flags)}")
    return out

print("="*126); print("SYSTEMATIC VALIDATION  -  pre-specified approved-drug benchmark (hardened engine)"); print("="*126)
OTV = ot_version()
print(f"  engine sha {ENGINE_SHA[:12]} | OpenTargets {OTV if OTV else '(version unavailable)'}")
R = run(PANEL)

known    = [r for r in R if r["known"] in ("inhibitor","activator")]
failed   = [r for r in R if r["verdict"] == "QUERY_FAILED"]
rescued  = [r for r in known if r["verdict"] in RESCUED]
refused  = [r for r in known if r["verdict"] in REFUSED]
labeled  = [r for r in known if r["verdict"] in LABELED]
decoyfl  = [r for r in known if r["verdict"] == "RECOVERED_CAVEAT_DECOY"]
correct  = [r for r in rescued if r["direction"] == r["known"]]
wrong    = [r for r in rescued if r["direction"] != r["known"]]
contested = [r for r in R if r["known"] == "contested"]
contested_ok = [r for r in contested if r["verdict"] in REFUSED]
labcorrect = [r for r in labeled if r["verdict"] == "LABEL_CONCORDANT"]

def tier_acc(t):
    g = [r for r in rescued if r["tier"] == t]
    return len([r for r in g if r["direction"] == r["known"]]), len(g)
hi_c, hi_n = tier_acc("HIGH"); st_c, st_n = tier_acc("STANDARD"); cav_c, cav_n = tier_acc("CAVEAT")
pqtl_wrong = [r["sym"] for r in wrong if r["tier"] == "HIGH"]
cov_den = len(rescued) + len(refused) + len(decoyfl)

print("\n" + "="*126); print("HEADLINE  (every panel target counted; failures shown, not hidden)"); print("="*126)
print(f"  panel: {len(known)} known-mechanism + {len(contested)} contested control" + (f"  | {len(failed)} QUERY_FAILED (excluded): {[r['sym'] for r in failed]}" if failed else ""))
print(f"  COVERAGE  rescued / GWAS-only total:    {line(len(rescued), cov_den)}   [refused {len(refused)}, decoy-flagged {len(decoyfl)}]")
print(f"  ACCURACY  rescued matching mechanism:   {line(len(correct), len(rescued))}")
print(f"    - HIGH  (pQTL-corroborated):          {line(hi_c, hi_n)}")
print(f"    - STANDARD (eQTL-only):               {line(st_c, st_n)}   <- documented decouplings (e.g. TYK2) live here")
if cav_n: print(f"    - CAVEAT (eQTL/pQTL conflict):        {line(cav_c, cav_n)}")
print(f"  RELIABILITY  wrong among pQTL-corroborated:  {len(pqtl_wrong)}  {pqtl_wrong if pqtl_wrong else '(none — no protein-backed call missed)'}")
print(f"  CONTROLS  contested refused (correct):  {len(contested_ok)}/{len(contested)}")
print(f"  LABELS    curated-label targets match:  {len(labcorrect)}/{len(labeled)}  (recovery not needed)")

print("\nFAILURE / EDGE CATALOGUE:")
print(f"  wrong rescues (recovered != drug):     {[(r['sym'],r['direction'],r['known'],r['tier']) for r in wrong]}")
print(f"  decoy-flagged (abundance!=function):   {[r['sym'] for r in decoyfl]}")
print(f"  refused (insufficient/conflicted):     {[(r['sym'],r['verdict']) for r in refused]}")
print(f"  query failures (transient/schema):     {[(r['sym'],(r.get('error') or '')[:60]) for r in failed]}")

manifest = {
  "engine_sha": ENGINE_SHA, "opentargets_version": OTV,
  "panel_sha": hashlib.sha256(json.dumps(PANEL, sort_keys=True).encode()).hexdigest()[:16],
  "n_known": len(known), "n_query_failed": len(failed),
  "n_rescued": len(rescued), "n_refused": len(refused), "n_decoy_flagged": len(decoyfl),
  "coverage": (len(rescued)/cov_den if cov_den else None), "coverage_ci": wilson(len(rescued), cov_den),
  "n_correct": len(correct), "accuracy_rescued": (len(correct)/len(rescued) if rescued else None),
  "accuracy_ci": wilson(len(correct), len(rescued)),
  "accuracy_by_tier": {"HIGH": [hi_c, hi_n], "STANDARD": [st_c, st_n], "CAVEAT": [cav_c, cav_n]},
  "wrong_among_pqtl_corroborated": pqtl_wrong,
  "wrong": [r["sym"] for r in wrong], "contested_refused": f"{len(contested_ok)}/{len(contested)}",
  "label_correct": f"{len(labcorrect)}/{len(labeled)}",
  "per_target": [{"sym":r["sym"],"efo":r["efo"],"verdict":r["verdict"],"dir":r["direction"],
                  "tier":r["tier"],"cons":r["consensus"],"sha":(r["sha"][:12] if r["sha"] else None)} for r in R],
}
manifest["manifest_sha"] = hashlib.sha256(json.dumps(manifest, sort_keys=True, separators=(",",":")).encode()).hexdigest()
for outdir in ("validation/results", "/kaggle/working/validation/results", "/kaggle/working/ftrepo/validation/results", "."):
    try:
        os.makedirs(outdir, exist_ok=True)
        open(os.path.join(outdir, "systematic_run.json"), "w").write(json.dumps(manifest, indent=2))
        print(f"\nwrote {outdir}/systematic_run.json   manifest sha {manifest['manifest_sha'][:12]}")
        break
    except Exception:
        continue

print("\nINTERPRET: coverage = how much of the direction-mute GWAS bulk the coloc layer reaches on a pre-committed")
print("set; accuracy = of those, how many match the drug, Wilson CI making small-n explicit; the tier split is the")
print("reliability axis (pQTL-corroborated = trustworthy; eQTL-only = where decoupling like TYK2 lives). Commit this")
print("manifest beside the figure; the writeup cites these numbers, not prose estimates.")
