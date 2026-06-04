import hashlib
import json
import time

import requests

OT = "https://api.platform.opentargets.org/api/v4/graphql"
H4_MIN  = 0.8     # colocalisation posterior floor
MAG     = 0.8     # a coloc votes only if |betaRatioSignAverage| >= MAG (drops the ~0 ambiguous ones)
FRAC    = 0.8     # >= this fraction of DECIDED loci must agree to call a direction
MIN_LOCI = 3      # need at least this many decided loci to vouch
GWAS = "gwas_credible_sets"
GENETIC_SOURCES = {"gwas_credible_sets","gene_burden","eva","eva_somatic","clingen",
                   "genomics_england","orphanet","gene2phenotype","uniprot_variants",
                   "uniprot_literature","ot_genetics_portal"}   # clinical_precedence excluded (circular)
MECH = {"inhibitor":"inhibitor","antagonist":"inhibitor","blocker":"inhibitor","degrader":"inhibitor",
        "activator":"activator","agonist":"activator","potentiator":"activator"}

def post(q, variables=None, label="", retries=3):
    for a in range(retries):
        try:
            r = requests.post(OT, json={"query": q, "variables": variables or {}}, timeout=60)
        except requests.RequestException as e:
            print(f"  !! net {label}: {e}"); return None
        if r.status_code == 200:
            d = r.json()
            if "errors" in d: print(f"  !! gql {label}: {json.dumps(d['errors'])[:240]}")
            return d
        if r.status_code in (429, 502, 503): time.sleep(2 ** a); continue
        print(f"  !! HTTP {r.status_code} {label}"); return None
    return None

def gwas_locus_ids(ensembl, efo, cap=40):   # cap for the validation sweep
    q = """query($e:String!,$f:String!){ target(ensemblId:$e){ approvedSymbol
      evidences(efoIds:[$f], datasourceIds:["gwas_credible_sets"], size:300){
        rows{ credibleSet { studyLocusId } } } } }"""
    d = post(q, {"e": ensembl, "f": efo}, label="gwas-loci")
    tgt = ((d or {}).get("data") or {}).get("target") or {}
    rows = ((tgt.get("evidences") or {}).get("rows")) or []
    ids = [(r.get("credibleSet") or {}).get("studyLocusId") for r in rows]
    return tgt.get("approvedSymbol"), [i for i in dict.fromkeys(ids) if i][:cap]

COLOC_Q = """query($id:String!){ credibleSet(studyLocusId:$id){ studyLocusId
  colocalisation(studyTypes:[eqtl, pqtl]){ rows {
    betaRatioSignAverage h4 rightStudyType
    otherStudyLocus { qtlGeneId study { target { id approvedSymbol } } } } } } }"""

def cis_signs(studyLocusId, ensembl, symbol, h4_min):
    """Strong directional cis colocs as [(sign,h4,qtltype)]: gene IS the target (qtlGeneId or
    study.target, Ensembl/symbol), h4>=floor, sign present. Excludes trans (PLA2G7-type)."""
    d = post(COLOC_Q, {"id": studyLocusId}, label=f"coloc:{studyLocusId[:10]}")
    cs = ((d or {}).get("data") or {}).get("credibleSet") or {}
    rows = ((cs.get("colocalisation") or {}).get("rows")) or []
    sym_u = symbol.upper() if symbol else None
    out = []
    for r in rows:
        o = r.get("otherStudyLocus") or {}
        tgt = ((o.get("study") or {}).get("target") or {})
        is_cis = (o.get("qtlGeneId") == ensembl) or (tgt.get("id") == ensembl) \
                 or (sym_u and (tgt.get("approvedSymbol") or "").upper() == sym_u)
        h4 = r.get("h4"); sign = r.get("betaRatioSignAverage")
        if not is_cis or sign is None: continue
        if isinstance(h4, (int, float)) and h4 < h4_min: continue
        qt = "pqtl" if "pqtl" in str(r.get("rightStudyType")).lower() else \
             ("eqtl" if "eqtl" in str(r.get("rightStudyType")).lower() else "other")
        out.append((sign, h4, qt))
    return out

def recover(ensembl, efo, h4_min=H4_MIN):
    """LOCUS-CONSENSUS recovery. Each GWAS locus -> +1/-1 from its strong cis colocs
    (|sign|>=MAG); a direction is called iff >=FRAC of >=MIN_LOCI decided loci agree."""
    sym, ids = gwas_locus_ids(ensembl, efo)
    loci_dir = {}; cpos = cneg = 0
    for lid in ids:
        strong = [s for (s, _, _) in cis_signs(lid, ensembl, sym, h4_min) if abs(s) >= MAG]
        p = sum(1 for s in strong if s > 0); n = sum(1 for s in strong if s < 0)
        cpos += p; cneg += n
        loci_dir[lid] = (1 if p > n else (-1 if n > p else 0))
    decided = [d for d in loci_dir.values() if d != 0]
    lpos = sum(1 for d in decided if d > 0); lneg = len(decided) - lpos
    frac = (max(lpos, lneg) / len(decided)) if decided else 0.0
    drug = None
    if len(decided) >= MIN_LOCI and frac >= FRAC:
        drug = "inhibitor" if lpos >= lneg else "activator"
    # strong-coloc-majority as a cross-check (the simpler, pseudo-replicating rule)
    cx = None
    if cpos + cneg:
        cf = max(cpos, cneg) / (cpos + cneg)
        cx = ("inhibitor" if cpos >= cneg else "activator") if cf >= FRAC else "conflicted"
    return {"symbol": sym, "recovered": drug, "consensus": round(frac, 2),
            "n_loci_total": len(loci_dir), "n_loci_decided": len(decided),
            "loci_pos": lpos, "loci_neg": lneg, "coloc_pos": cpos, "coloc_neg": cneg,
            "crosscheck_strong_majority": cx, "loci": sorted([l for l, d in loci_dir.items() if d != 0])}

# ---- curated label (genetic only) so labelled targets bypass recovery ------
DT = {"lof":"down","loss_of_function":"down","gof":"up","gain_of_function":"up","decrease":"down","increase":"up"}
def desired_from_label(dot, dotr):
    t = DT.get(str(dot).strip().lower()) if dot else None
    tr = str(dotr).strip().lower() if dotr else None
    if not t or tr not in ("risk", "protect"): return None
    if (t, tr) in (("up","risk"),("down","protect")): return "inhibitor"
    if (t, tr) in (("down","risk"),("up","protect")): return "activator"
    return None

def labelled_direction(ensembl, efo):
    q = """query($e:String!,$f:String!){ target(ensemblId:$e){
      evidences(efoIds:[$f], size:1000){ rows{ datasourceId directionOnTarget directionOnTrait } } } }"""
    d = post(q, {"e": ensembl, "f": efo}, label="labels")
    rows = (((((d or {}).get("data") or {}).get("target") or {}).get("evidences") or {}).get("rows")) or []
    votes = {}
    for r in rows:
        ds = r.get("datasourceId")
        if ds not in GENETIC_SOURCES or ds == GWAS: continue
        des = desired_from_label(r.get("directionOnTarget"), r.get("directionOnTrait"))
        if des: votes[des] = votes.get(des, 0) + 1
    if not votes: return None, 0
    return max(votes, key=votes.get), sum(votes.values())

def verdict(ensembl, efo, mechanism, h4_min=H4_MIN):
    claim = MECH.get(str(mechanism).strip().lower())
    if claim is None: raise ValueError(f"bad mechanism {mechanism!r}")
    lab, nlab = labelled_direction(ensembl, efo)
    if lab:
        res = {"verdict": "LABEL_CONCORDANT" if lab == claim else "LABEL_DISCORDANT",
               "vouches": lab == claim, "source": "curated-label", "confidence": "high",
               "direction": lab, "detail": f"{nlab} curated genetic label rows -> {lab}"}
    else:
        rec = recover(ensembl, efo, h4_min)
        pct = int(round(rec["consensus"] * 100))
        if rec["n_loci_decided"] < MIN_LOCI:
            res = {"verdict": "INSUFFICIENT_DIRECTION", "vouches": False, "source": "none",
                   "confidence": "n/a", "direction": None, "recovery": rec,
                   "detail": f"only {rec['n_loci_decided']} decided cis-coloc locus/loci (<{MIN_LOCI})",
                   "falsifier": "too few cis molQTL coloc loci over H4>=%.2f/|sign|>=%.1f; an MR or rare-variant burden test is the entry ticket." % (h4_min, MAG)}
        elif rec["recovered"] is None:
            res = {"verdict": "RECOVERY_CONFLICTED", "vouches": False, "source": "coloc",
                   "confidence": "n/a", "direction": None, "recovery": rec,
                   "detail": f"locus consensus {pct}% over {rec['n_loci_decided']} loci is below {int(FRAC*100)}%",
                   "falsifier": "loci disagree on direction; a trait-specific MR would adjudicate."}
        else:
            concord = rec["recovered"] == claim
            res = {"verdict": "RECOVERED_CONCORDANT" if concord else "RECOVERED_DISCORDANT",
                   "vouches": concord, "source": "coloc-recovered",
                   "confidence": f"moderate (coloc-derived; {pct}% locus consensus over {rec['n_loci_decided']} loci)",
                   "direction": rec["recovered"], "recovery": rec,
                   "falsifier": "coloc-inferred direction; a rare-variant burden test or MR would upgrade it to direct."}
    payload = {"rule": "R3+recovery/locus-consensus", "ensembl": ensembl, "efo": efo, "mechanism": claim,
               "verdict": res["verdict"], "vouches": res["vouches"], "direction": res.get("direction"),
               "loci": res.get("recovery", {}).get("loci", [])}
    res["sha256"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    res["ensembl"], res["efo"], res["claim"] = ensembl, efo, claim
    return res


# ============================================================================
# CHOLESTEROL DEEP-DIVE  -  dual-endpoint (CAD direct + LDL) on the lipid pathway
#   (1) CONCORDANCE: where a target resolves on BOTH coronary-artery-disease and
#       LDL, do the two directions agree? (the faithful-bridge test)
#   (2) ACTIVATOR-direction live check: an unlabelled, LDL-LOWERING target that
#       refuses on CAD should rescue via LDL as ACTIVATOR (the gap the demo left).
# Run on Kaggle, internet ON. ~6-10 min (two endpoints x a small lipid panel).
# ============================================================================
LDL = "EFO_0004611"; CAD = "EFO_0001645"
LDL_RISK_SIGN = +1          # higher LDL -> higher CAD risk: identity map (HDL would flip)
def flip(d): return {"inhibitor":"activator","activator":"inhibitor"}.get(d, d)

def ldl_endorsement(ens):
    r = recover(ens, LDL)                      # direction treating LDL as if it were risk
    raw = r.get("recovered")
    endo = None if raw is None else (raw if LDL_RISK_SIGN > 0 else flip(raw))
    return endo, r

# symbol, ensembl, mechanism-from-drug-or-genetics, note   (more gene product -> ? LDL)
LIPID = [
  ("PCSK9",  "ENSG00000169174","inhibitor","evolocumab; more->higher LDL"),
  ("HMGCR",  "ENSG00000113161","inhibitor","statins; more->higher LDL"),
  ("NPC1L1", "ENSG00000015520","inhibitor","ezetimibe; more->higher LDL"),
  ("APOB",   "ENSG00000084674","inhibitor","mipomersen; more->higher LDL"),
  ("MTTP",   "ENSG00000138823","inhibitor","lomitapide; more->higher LDL"),
  ("CETP",   "ENSG00000087237","inhibitor","anacetrapib (borderline outcome)"),
  ("LPA",    "ENSG00000198670","inhibitor","pelacarsen (trials); more->higher risk"),
  ("LDLR",   "ENSG00000130164","activator","FH gene; more->LOWER LDL (often labelled)"),
  ("LDLRAP1","ENSG00000157978","activator","LDLR adaptor; more->LOWER LDL"),
  ("ABCG5",  "ENSG00000138075","activator","sterol efflux; more->LOWER LDL"),
  ("ABCG8",  "ENSG00000143921","activator","sterol efflux; more->LOWER LDL"),
]
REFUSED = {"RECOVERY_CONFLICTED","INSUFFICIENT_DIRECTION"}

def run_lipid(panel):
    rows=[]
    print(f"{'target':<9}{'mech':<10}{'CAD verdict':<24}{'CADdir':<10}{'LDLdir':<10}{'LDLloci':>7}  status")
    print("-"*96)
    for sym,ens,mech,note in panel:
        vc = verdict(ens, CAD, mech)                 # disease endpoint
        cad_dir = vc.get("direction")
        ldl_dir, r = ldl_endorsement(ens)            # biomarker endpoint (identity-mapped)
        nloci = r.get("n_loci_decided")
        # classify
        if cad_dir and ldl_dir:
            status = "CONCORDANT" if cad_dir==ldl_dir else "DISCORDANT(bridge!)"
        elif vc["verdict"] in REFUSED and ldl_dir:
            status = f"RESCUED->{ldl_dir}"
        elif vc["verdict"] in REFUSED:
            status = "refused both"
        else:
            status = "direct only"
        rsym = r.get("symbol")
        sflag = "" if (rsym is None or rsym==sym) else f" !{rsym}"
        rows.append({"sym":sym,"ens":ens,"mech":mech,"cad_verdict":vc["verdict"],"cad_dir":cad_dir,
                     "ldl_dir":ldl_dir,"nloci":nloci,"status":status,
                     "sha":hashlib.sha256(f"{ens}|{CAD}|{LDL}|{ldl_dir}".encode()).hexdigest()[:12]})
        print(f"{sym:<9}{mech:<10}{vc['verdict']:<24}{str(cad_dir):<10}{str(ldl_dir):<10}{str(nloci):>7}  {status}{sflag}")
    return rows

print("="*96); print("CHOLESTEROL DEEP-DIVE  -  CAD endpoint vs LDL endpoint on the lipid pathway"); print("="*96)
R = run_lipid(LIPID)

both   = [r for r in R if r["cad_dir"] and r["ldl_dir"]]
conc   = [r for r in both if r["cad_dir"]==r["ldl_dir"]]
disc   = [r for r in both if r["cad_dir"]!=r["ldl_dir"]]
rescue = [r for r in R if r["cad_verdict"] in REFUSED and r["ldl_dir"]]
resc_ok= [r for r in rescue if r["ldl_dir"]==r["mech"]]
act_resc = [r for r in rescue if r["ldl_dir"]=="activator"]
inh_resc = [r for r in rescue if r["ldl_dir"]=="inhibitor"]
def pct(a,b): return f"{100*a/b:.0f}%" if b else "n/a"

print("\n"+"="*96); print("HEADLINE"); print("="*96)
print(f"  (1) BRIDGE CONCORDANCE  (resolves on BOTH CAD and LDL): {len(conc)}/{len(both)} agree = {pct(len(conc),len(both))}")
if disc: print(f"      DISCORDANT (bridge problem, inspect): {[(d['sym'],d['cad_dir'],d['ldl_dir']) for d in disc]}")
print(f"  (2) CAD refusals rescued via LDL: {len(rescue)}   matching mechanism: {len(resc_ok)}/{len(rescue)} = {pct(len(resc_ok),len(rescue))}")
print(f"      inhibitor-direction rescues: {[r['sym'] for r in inh_resc]}")
print(f"      ACTIVATOR-direction rescues: {[r['sym'] for r in act_resc]}   <- closes the both-ways gap if non-empty")
print("\nPROVENANCE (LDL-route verdicts):")
for r in R:
    if r["ldl_dir"]: print(f"  {r['sym']:<9}{r['status']:<22} LDL-dir={r['ldl_dir']:<10} loci={r['nloci']}  sha={r['sha']}")
print("\nREAD: concordance where both endpoints speak tests that the LDL bridge is FAITHFUL (not just")
print("a coverage trick). An activator-direction rescue proves the sign map works both ways on real")
print("data -- the thing the LDLR control couldn't show because it had a label and bypassed recovery.")
