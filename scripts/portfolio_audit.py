#!/usr/bin/env python3
"""Production portfolio audit. Run: python scripts/portfolio_audit.py <claims_dir> <reports_dir>"""
from __future__ import annotations
import sys, csv, json, subprocess
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter

TIER_MAP = {
    "SURVIVED":               "TIER_3_PROCEED",
    "FALSIFIED_WITH_CAVEATS": "TIER_2_INVESTIGATE",
    "INSUFFICIENT_DATA":      "TIER_2_INVESTIGATE",
    "FALSIFIED":              "TIER_1_DROP",
    "ERROR":                  "TIER_0_FIX_INPUT",
}

def extract_actions(report):
    fired, substantive, abstained = [], [], []
    for r in report.get("per_rule", []):
        rid = r.get("rule_id", "?")
        if r.get("status") == "falsified":
            fired.append((rid, r.get("falsification_tier") or "?",
                          (r.get("falsification_experiment") or "")[:160]))
        elif r.get("status") == "abstained":
            abstained.append(rid)
        for c in r.get("caveats", []):
            if c.get("kind") == "substantive":
                substantive.append((rid, c.get("text", "")[:160]))
    return {"fired": fired, "substantive": substantive, "abstained": abstained}

def format_actions(a):
    if a["fired"]:
        return "; ".join(f"{rid} FALSIFIED @ {t}: {e}" if e else f"{rid} FALSIFIED @ {t}"
                         for rid, t, e in a["fired"])
    if a["substantive"]:
        return "; ".join(f"{rid}: {txt}" for rid, txt in a["substantive"])
    if a["abstained"]:
        return f"data gap: {', '.join(a['abstained'])} abstained (no input)"
    return "all rules cleared"

def audit_one(claim_path, reports_dir, ft_root):
    out_json = reports_dir / f"{claim_path.stem}.json"
    r = subprocess.run([sys.executable, str(ft_root / "run_audit.py"),
                        str(claim_path), "--no-live", "--json-out", str(out_json)],
                       cwd=ft_root, capture_output=True, text=True)
    if r.returncode == 5 or not out_json.exists():
        err = next((l for l in r.stderr.strip().splitlines() if l.startswith("ERROR:")),
                   r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "unknown")
        return {"claim_file": claim_path.name, "target": "(unparsed)",
                "verdict": "ERROR", "risk_tier": TIER_MAP["ERROR"], "score": "",
                "caveats": 0, "fired_rules": "", "cheapest_fix": "-",
                "action_items": err[:200], "report_json": ""}
    d = json.loads(out_json.read_text())
    a = extract_actions(d)
    c = d.get("cheapest_falsification")
    return {"claim_file": claim_path.name, "target": d["claim"]["target_symbol"],
            "verdict": d["verdict"], "risk_tier": TIER_MAP.get(d["verdict"], "TIER_2_INVESTIGATE"),
            "score": d["score"], "caveats": d["substantive_caveat_count"],
            "fired_rules": ",".join(rid for rid, _, _ in a["fired"]) or "-",
            "cheapest_fix": f"{c['rule_id']} @ {c['tier']}" if c else "-",
            "action_items": format_actions(a), "report_json": str(out_json)}

def main():
    if len(sys.argv) < 3:
        print("usage: portfolio_audit.py <claims_dir> <reports_dir>"); sys.exit(2)
    claims_dir, reports_dir = Path(sys.argv[1]), Path(sys.argv[2])
    ft_root = Path(__file__).parent.parent
    reports_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(p for p in claims_dir.glob("*.yaml")
                   if p.name != "paralog_map.yaml")
    if not files: print(f"no claims in {claims_dir}"); sys.exit(2)
    rows = [audit_one(f, reports_dir, ft_root) for f in files]
    glyphs = {"SURVIVED":"PASS","FALSIFIED_WITH_CAVEATS":"WARN",
              "FALSIFIED":"FAIL","INSUFFICIENT_DATA":"????","ERROR":"ERR!"}
    for r in rows:
        print(f"  [{glyphs.get(r['verdict'],'?'):4}]  {r['target']:<28}  "
              f"{r['verdict']:<25} tier={r['risk_tier']:<20} caveats={r['caveats']}")
        print(f"           {r['action_items'][:140]}")
    with open(reports_dir / "portfolio_decision.csv", "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    manifest = {"generated_utc": datetime.now(timezone.utc).isoformat(),
                "claim_count": len(rows),
                "verdict_summary": dict(Counter(r["verdict"] for r in rows)),
                "tier_summary":    dict(Counter(r["risk_tier"] for r in rows)),
                "rows": rows}
    (reports_dir / "portfolio_manifest.json").write_text(json.dumps(manifest, indent=2))
    has_drop  = any(r["risk_tier"] == "TIER_1_DROP" for r in rows)
    has_error = any(r["risk_tier"] == "TIER_0_FIX_INPUT" for r in rows)
    has_warn  = any(r["risk_tier"] == "TIER_2_INVESTIGATE" for r in rows)
    sys.exit(2 if has_drop else (5 if has_error else (1 if has_warn else 0)))

if __name__ == "__main__":
    main()
