#!/usr/bin/env python3
"""Production portfolio audit.

Usage:
    python scripts/portfolio_audit.py <claims_dir> <reports_dir>

Audits every ``*.yaml`` in ``claims_dir`` (excluding ``paralog_map.yaml``)
and writes three output files into ``reports_dir``:

- ``portfolio_decision.csv``  for spreadsheets / program managers
- ``portfolio_decision.md``   for PR descriptions / wiki pages
- ``portfolio_manifest.json`` for downstream automation

Each row carries an ``action_items`` column built from (in priority order)
fired rules, substantive caveats, abstained rules, or a clean signal.

Exit code reflects worst tier:
    0  all TIER_3_PROCEED
    1  TIER_2_INVESTIGATE present
    2  TIER_1_DROP present
    5  TIER_0_FIX_INPUT present (broken YAML)
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

TIER_MAP = {
    "SURVIVED":               "TIER_3_PROCEED",
    "FALSIFIED_WITH_CAVEATS": "TIER_2_INVESTIGATE",
    "INSUFFICIENT_DATA":      "TIER_2_INVESTIGATE",
    "FALSIFIED":              "TIER_1_DROP",
    "ERROR":                  "TIER_0_FIX_INPUT",
}

GLYPHS = {
    "SURVIVED":               "PASS",
    "FALSIFIED_WITH_CAVEATS": "WARN",
    "FALSIFIED":              "FAIL",
    "INSUFFICIENT_DATA":      "????",
    "ERROR":                  "ERR!",
}

TIER_ORDER = {"TIER_1_DROP": 0, "TIER_0_FIX_INPUT": 1,
              "TIER_2_INVESTIGATE": 2, "TIER_3_PROCEED": 3}


def extract_actions(report: dict) -> dict:
    """Bucket per_rule entries by category for action_items rendering."""
    fired, substantive, abstained = [], [], []
    for r in report.get("per_rule", []):
        rid = r.get("rule_id", "?")
        status = r.get("status")
        if status == "falsified":
            fired.append((
                rid,
                r.get("falsification_tier") or "?",
                (r.get("falsification_experiment") or "")[:160],
            ))
        elif status == "abstained":
            abstained.append(rid)
        for c in r.get("caveats", []):
            if c.get("kind") == "substantive":
                substantive.append((rid, c.get("text", "")[:160]))
    return {"fired": fired, "substantive": substantive, "abstained": abstained}


def format_actions(a: dict) -> str:
    """Render the action_items string in priority order."""
    if a["fired"]:
        return "; ".join(
            f"{rid} FALSIFIED @ {t}: {e}" if e else f"{rid} FALSIFIED @ {t}"
            for rid, t, e in a["fired"]
        )
    if a["substantive"]:
        return "; ".join(f"{rid}: {txt}" for rid, txt in a["substantive"])
    if a["abstained"]:
        return f"data gap: {', '.join(a['abstained'])} abstained (no input)"
    return "all rules cleared"


def audit_one(claim_path: Path, reports_dir: Path, ft_root: Path) -> dict:
    out_json = reports_dir / f"{claim_path.stem}.json"
    cmd = [sys.executable, str(ft_root / "run_audit.py"),
           str(claim_path), "--no-live", "--json-out", str(out_json)]
    r = subprocess.run(cmd, cwd=ft_root, capture_output=True, text=True)

    if r.returncode == 5 or not out_json.exists():
        err_lines = r.stderr.strip().splitlines() if r.stderr.strip() else ["unknown error"]
        err = next((l for l in err_lines if l.startswith("ERROR:")), err_lines[-1])
        return {
            "claim_file":   claim_path.name,
            "target":       "(unparsed)",
            "verdict":      "ERROR",
            "risk_tier":    TIER_MAP["ERROR"],
            "score":        "",
            "caveats":      0,
            "fired_rules":  "",
            "cheapest_fix": "-",
            "action_items": err[:200],
            "report_json":  "",
        }

    d = json.loads(out_json.read_text())
    actions = extract_actions(d)
    cheapest = d.get("cheapest_falsification")
    return {
        "claim_file":   claim_path.name,
        "target":       d["claim"]["target_symbol"],
        "verdict":      d["verdict"],
        "risk_tier":    TIER_MAP.get(d["verdict"], "TIER_2_INVESTIGATE"),
        "score":        d["score"],
        "caveats":      d["substantive_caveat_count"],
        "fired_rules":  ",".join(rid for rid, _, _ in actions["fired"]) or "-",
        "cheapest_fix": (f"{cheapest['rule_id']} @ {cheapest['tier']}"
                        if cheapest else "-"),
        "action_items": format_actions(actions),
        "report_json":  str(out_json),
    }


def write_csv(rows: list[dict], path: Path) -> None:
    with open(path, "w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def write_markdown(rows: list[dict], path: Path) -> None:
    tier_counts = Counter(r["risk_tier"] for r in rows)
    with open(path, "w") as fp:
        fp.write("# Portfolio audit decision sheet\n\n")
        fp.write(f"_Generated {datetime.now(timezone.utc).isoformat()}_\n\n")
        fp.write("## Summary by risk tier\n\n")
        for tier in ["TIER_1_DROP", "TIER_0_FIX_INPUT",
                     "TIER_2_INVESTIGATE", "TIER_3_PROCEED"]:
            if tier_counts.get(tier):
                fp.write(f"- **{tier}**: {tier_counts[tier]}\n")
        fp.write("\n## Per-claim detail (highest risk first)\n\n")
        fp.write("| Tier | Target | Verdict | Caveats | Fired | Action items |\n")
        fp.write("|---|---|---|---|---|---|\n")
        for r in sorted(rows, key=lambda x: TIER_ORDER.get(x["risk_tier"], 9)):
            fp.write(f"| `{r['risk_tier']}` | **{r['target']}** | {r['verdict']} "
                     f"| {r['caveats']} | `{r['fired_rules']}` "
                     f"| {r['action_items'][:200]} |\n")


def main():
    if len(sys.argv) < 3:
        print("usage: portfolio_audit.py <claims_dir> <reports_dir>",
              file=sys.stderr)
        sys.exit(2)

    claims_dir = Path(sys.argv[1])
    reports_dir = Path(sys.argv[2])
    ft_root = Path(__file__).resolve().parent.parent
    reports_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(p for p in claims_dir.glob("*.yaml")
                   if p.name != "paralog_map.yaml")
    if not files:
        print(f"no claim YAMLs in {claims_dir}", file=sys.stderr)
        sys.exit(2)

    rows = []
    for f in files:
        row = audit_one(f, reports_dir, ft_root)
        glyph = GLYPHS.get(row["verdict"], "?")
        print(f"  [{glyph:4}]  {row['target']:<28}  "
              f"{row['verdict']:<25} tier={row['risk_tier']:<20} "
              f"caveats={row['caveats']}")
        print(f"           {row['action_items'][:140]}")
        rows.append(row)

    write_csv(rows, reports_dir / "portfolio_decision.csv")
    write_markdown(rows, reports_dir / "portfolio_decision.md")
    (reports_dir / "portfolio_manifest.json").write_text(json.dumps({
        "generated_utc":    datetime.now(timezone.utc).isoformat(),
        "claim_count":      len(rows),
        "verdict_summary":  dict(Counter(r["verdict"] for r in rows)),
        "tier_summary":     dict(Counter(r["risk_tier"] for r in rows)),
        "rows":             rows,
    }, indent=2))

    print(f"\nCSV:  {reports_dir / 'portfolio_decision.csv'}")
    print(f"MD:   {reports_dir / 'portfolio_decision.md'}")
    print(f"JSON: {reports_dir / 'portfolio_manifest.json'}")

    has_drop  = any(r["risk_tier"] == "TIER_1_DROP" for r in rows)
    has_error = any(r["risk_tier"] == "TIER_0_FIX_INPUT" for r in rows)
    has_warn  = any(r["risk_tier"] == "TIER_2_INVESTIGATE" for r in rows)
    sys.exit(2 if has_drop else (5 if has_error else (1 if has_warn else 0)))


if __name__ == "__main__":
    main()
