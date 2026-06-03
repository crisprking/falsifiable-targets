"""
Pre-specified approved-drug benchmark for direction_audit (no cherry-picking).
Run:  python validation/systematic_benchmark.py   (needs internet; ~5-8 min)
Every panel target is reported; mechanism direction is cited from the approved drug;
documented abundance!=function decoys (IL6R) are FLAGGED, not scored as plain wrong.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from direction_audit import post, verdict
import json, hashlib

# documented: measured abundance decoupled from function (soluble/decoy isoform). curated, not auto-detected.
DECOY = {"ENSG00000160712"}   # IL6R (sIL6R decoy; rs2228145)

# symbol, ensembl, efo, known_mechanism, drug (citation for the mechanism call)
PANEL = [
  # --- cardiometabolic ---
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
  # --- immunology ---
  ("IL6R",    "ENSG00000160712","EFO_0000685",  "inhibitor","tocilizumab [decoy]"),
  ("TNF",     "ENSG00000232810","EFO_0000685",  "inhibitor","adalimumab"),
  ("IL23R",   "ENSG00000162594","EFO_0003767",  "inhibitor","via IL23 (ustekinumab/risankizumab)"),
  ("IL12B",   "ENSG00000113302","EFO_0003767",  "inhibitor","ustekinumab"),
  ("IL4R",    "ENSG00000077238","EFO_0000274",  "inhibitor","dupilumab"),
  ("TYK2",    "ENSG00000105397","EFO_0000676",  "inhibitor","deucravacitinib"),
  ("TSLP",    "ENSG00000145777","EFO_0000270",  "inhibitor","tezepelumab"),
  # --- contested control (should refuse) ---
  ("GIPR",    "ENSG00000010310","MONDO_0005148","contested","agonist-vs-antagonist debate"),
]

RESCUED = {"RECOVERED_CONCORDANT","RECOVERED_DISCORDANT"}
REFUSED = {"RECOVERY_CONFLICTED","INSUFFICIENT_DIRECTION"}
LABELED = {"LABEL_CONCORDANT","LABEL_DISCORDANT"}

def efo_name(efo):
    d = post('query($f:String!){ disease(efoId:$f){ name } }', {"f": efo}, "dn")
    return (((d or {}).get("data") or {}).get("disease") or {}).get("name") or "?? EFO NOT FOUND"

def run(panel):
    out = []
    print(f"{'target':<9}{'disease':<30}{'known':<10}{'verdict':<24}{'dir':<10}{'c%':>4}{'loc':>5}  src/flags")
    print("-" * 122)
    for sym, ens, efo, known, drug in panel:
        nm = efo_name(efo)
        claim = known if known in ("inhibitor","activator") else "inhibitor"
        v = verdict(ens, efo, claim)
        rec = v.get("recovery", {})
        resolved = rec.get("symbol")  # approvedSymbol from gwas_locus_ids (None if no gwas evidence path)
        # post-hoc decoy flag (curated): a coloc rescue on a documented decoy is NOT a mechanism call
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
                    "consensus":rec.get("consensus"),"n_loci":rec.get("n_loci_decided"),
                    "source":v["source"],"decoy":(ens in DECOY),"sha":v["sha256"]})
        print(f"{sym:<9}{nm[:29]:<30}{known:<10}{v['verdict']:<24}{str(v.get('direction')):<10}"
              f"{str(c):>4}{str(nl):>5}  {v['source']} {symflag} {efoflag} {'/'.join(flags)}")
    return out

print("="*122); print("SYSTEMATIC VALIDATION  -  pre-specified approved-drug benchmark"); print("="*122)
R = run(PANEL)

known   = [r for r in R if r["known"] in ("inhibitor","activator")]
rescued = [r for r in known if r["verdict"] in RESCUED]
refused = [r for r in known if r["verdict"] in REFUSED]
labeled = [r for r in known if r["verdict"] in LABELED]
decoyfl = [r for r in known if r["verdict"] == "RECOVERED_CAVEAT_DECOY"]
correct = [r for r in rescued if r["direction"] == r["known"]]
wrong   = [r for r in rescued if r["direction"] != r["known"]]
contested = [r for r in R if r["known"] == "contested"]
contested_ok = [r for r in contested if r["verdict"] in REFUSED]
labcorrect = [r for r in labeled if r["verdict"] == "LABEL_CONCORDANT"]

cov_den = len(rescued) + len(refused) + len(decoyfl)
def pct(a,b): return f"{(100*a/b):.0f}%" if b else "n/a"

print("\n" + "="*122); print("HEADLINE  (every panel target counted; nothing dropped)"); print("="*122)
print(f"  panel: {len(known)} known-mechanism + {len(contested)} contested control")
print(f"  COVERAGE  (GWAS-only rescued / GWAS-only total):   {len(rescued)}/{cov_den} = {pct(len(rescued),cov_den)}"
      f"   [refused honestly {len(refused)}, decoy-flagged {len(decoyfl)}]")
print(f"  ACCURACY  rescued matching approved mechanism:     {len(correct)}/{len(rescued)} = {pct(len(correct),len(rescued))}")
print(f"  CONTROLS  contested refused (correct):             {len(contested_ok)}/{len(contested)}")
print(f"  LABELS    curated-label targets matching mech:     {len(labcorrect)}/{len(labeled)}  (recovery not needed)")

print("\nFAILURE / EDGE CATALOGUE:")
print(f"  wrong rescues (recovered != drug):     {[(r['sym'],r['direction'],r['known']) for r in wrong]}")
print(f"  decoy-flagged (abundance!=function):   {[r['sym'] for r in decoyfl]}")
print(f"  refused-insufficient/thin/conflicted:  {[(r['sym'],r['verdict']) for r in refused]}")
print(f"  symbol/EFO sanity flags above ('!','EFO?') -> inspect any that fired.")

manifest = {
  "panel_sha": hashlib.sha256(json.dumps(PANEL, sort_keys=True).encode()).hexdigest()[:16],
  "n_known": len(known), "n_rescued": len(rescued), "n_refused": len(refused),
  "n_decoy_flagged": len(decoyfl), "coverage": (len(rescued)/cov_den if cov_den else None),
  "n_correct": len(correct), "accuracy_rescued": (len(correct)/len(rescued) if rescued else None),
  "wrong": [r["sym"] for r in wrong], "contested_refused": f"{len(contested_ok)}/{len(contested)}",
  "label_correct": f"{len(labcorrect)}/{len(labeled)}",
  "per_target": [{"sym":r["sym"],"efo":r["efo"],"verdict":r["verdict"],"dir":r["direction"],
                  "cons":r["consensus"],"sha":r["sha"][:12]} for r in R],
}
print("\n" + "="*122); print("PROVENANCE (commit beside the figure):"); print(json.dumps(manifest, indent=2))
print("\nINTERPRET: coverage = how much of the direction-mute GWAS bulk the coloc layer reaches on a")
print("pre-specified set; accuracy = of those, how many match the drug; wrong+decoy = the honest error")
print("budget. A pre-committed denominator with failures shown is the cherry-pick-proof version of the claim.")
