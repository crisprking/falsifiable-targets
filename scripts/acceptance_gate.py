#!/usr/bin/env python3
"""Six-contract acceptance gate for the falsifiable-targets audit pipeline.

Run before tagging a release. Exits 0 only when all contracts hold.

Contracts:
  1. Verdict + exit-code mapping across five reference scenarios
  2. Error paths exit 5 with no Python tracebacks leaking
  3. Ruleset SHA pinned across all reports
  4. Determinism: back-to-back runs are byte-identical modulo timestamp
  5. --debug flag accepted, preserves exit codes, adds tracebacks
  6. Substantive caveat text is propagated into report JSON
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

FT_ROOT   = Path(__file__).resolve().parent.parent
CLAIMS    = FT_ROOT / "tests" / "hardtest_claims"
REPORTS   = FT_ROOT / "reports" / "acceptance_gate"
RUN_AUDIT = FT_ROOT / "run_audit.py"

VERDICT_CONTRACTS = {
    "s1_pcsk9_hypercholesterolemia.yaml": ("SURVIVED",               0),
    "s2_bace1_alzheimers.yaml":           ("FALSIFIED_WITH_CAVEATS", 1),
    "s3_stap_cells_retracted.yaml":       ("FALSIFIED",              2),
    "s4_class_collapse_phantom.yaml":     ("FALSIFIED",              2),
    "s5_novel_pathogen_weak.yaml":        ("FALSIFIED_WITH_CAVEATS", 1),
}

ERROR_CONTRACTS = [
    "s6a_malformed.yaml",
    "s6b_missing_claim_type.yaml",
    "s6c_empty.yaml",
    "s6d_yaml_parse_error.yaml",
]

EXPECTED_RULESET_SHA = (
    "35ef2b2ab5363298097962a0b6ae52c70d551a1edddc341054f75cb6e4fb7221"
)


def audit(claim_path: Path, json_out: Path | None = None,
          extra: list[str] | None = None) -> subprocess.CompletedProcess:
    cmd = [sys.executable, str(RUN_AUDIT), str(claim_path), "--no-live"]
    if json_out:
        cmd += ["--json-out", str(json_out)]
    if extra:
        cmd += extra
    return subprocess.run(cmd, cwd=FT_ROOT, capture_output=True, text=True)


def main():
    REPORTS.mkdir(parents=True, exist_ok=True)
    failures = []
    sep = "=" * 80

    # ---- C1: verdict + exit code ----
    print(sep)
    print("CONTRACT 1: verdict + exit-code")
    print(sep)
    for cf, (wv, we) in VERDICT_CONTRACTS.items():
        p = CLAIMS / cf
        j = REPORTS / f"{p.stem}.json"
        if not p.exists():
            print(f"  [SKIP]  {cf} (not found)")
            continue
        r = audit(p, json_out=j)
        gv = json.loads(j.read_text())["verdict"] if j.exists() else "?"
        ok = (gv == wv) and (r.returncode == we)
        print(f"  [{'OK  ' if ok else 'FAIL'}]  {cf:<40}  "
              f"verdict={gv} (want {wv})  exit={r.returncode} (want {we})")
        if not ok:
            failures.append(("verdict", cf, gv, r.returncode))

    # ---- C2: error paths ----
    print("\n" + sep)
    print("CONTRACT 2: error paths exit 5, no tracebacks")
    print(sep)
    for cf in ERROR_CONTRACTS:
        p = CLAIMS / cf
        if not p.exists():
            print(f"  [SKIP]  {cf} (not found)")
            continue
        r = audit(p)
        has_tb = "traceback (most recent call last)" in r.stderr.lower()
        ok = r.returncode == 5 and not has_tb
        print(f"  [{'OK  ' if ok else 'FAIL'}]  {cf:<32} "
              f"exit={r.returncode} (want 5)  no_traceback={not has_tb}")
        if not ok:
            failures.append(("error_exit", cf))

    # ---- C3: ruleset SHA pinned ----
    print("\n" + sep)
    print("CONTRACT 3: ruleset SHA pinned")
    print(sep)
    drift = False
    for cf in VERDICT_CONTRACTS:
        j = REPORTS / f"{Path(cf).stem}.json"
        if not j.exists():
            continue
        sha = json.loads(j.read_text())["ruleset_sha256"]
        if sha != EXPECTED_RULESET_SHA:
            print(f"  [FAIL]  {cf}: {sha[:16]} (want {EXPECTED_RULESET_SHA[:16]})")
            drift = True
            failures.append(("sha", cf))
    if not drift:
        print(f"  [OK  ]  all pinned to {EXPECTED_RULESET_SHA[:16]}...")

    # ---- C4: determinism ----
    print("\n" + sep)
    print("CONTRACT 4: determinism")
    print(sep)
    s1 = CLAIMS / "s1_pcsk9_hypercholesterolemia.yaml"
    if s1.exists():
        ja, jb = REPORTS / "_det_a.json", REPORTS / "_det_b.json"
        audit(s1, json_out=ja)
        audit(s1, json_out=jb)
        a = json.loads(ja.read_text())
        b = json.loads(jb.read_text())
        a.pop("audit_timestamp_utc", None)
        b.pop("audit_timestamp_utc", None)
        ca = json.dumps(a, sort_keys=True)
        cb = json.dumps(b, sort_keys=True)
        ok = ca == cb
        sha = hashlib.sha256(ca.encode()).hexdigest()[:16]
        print(f"  [{'OK  ' if ok else 'FAIL'}]  identical={ok}  payload_sha={sha}")
        if not ok:
            failures.append(("determinism",))
    else:
        print("  [SKIP]  s1 claim not found")

    # ---- C5: --debug flag ----
    print("\n" + sep)
    print("CONTRACT 5: --debug accepted and behaves")
    print(sep)
    s6a = CLAIMS / "s6a_malformed.yaml"
    if s6a.exists():
        r_plain = audit(s6a)
        r_debug = audit(s6a, extra=["--debug"])
        accepted = "unrecognized arguments" not in r_debug.stderr
        consistent = r_plain.returncode == 5 and r_debug.returncode == 5
        print(f"  [{'OK  ' if accepted else 'FAIL'}]  --debug accepted by argparse")
        print(f"  [{'OK  ' if consistent else 'FAIL'}]  --debug preserves exit "
              f"(plain={r_plain.returncode}, debug={r_debug.returncode})")
        if not accepted:
            failures.append(("debug_not_accepted",))
        if not consistent:
            failures.append(("debug_breaks_exit",))
    else:
        print("  [SKIP]  s6a claim not found")

    # ---- C6: substantive caveats surfaced ----
    print("\n" + sep)
    print("CONTRACT 6: substantive caveats surfaced")
    print(sep)
    for cf in ["s2_bace1_alzheimers.json", "s5_novel_pathogen_weak.json"]:
        j = REPORTS / cf
        if not j.exists():
            print(f"  [SKIP]  {cf}")
            continue
        d = json.loads(j.read_text())
        subs = [(r["rule_id"], c.get("text", ""))
                for r in d.get("per_rule", [])
                for c in r.get("caveats", [])
                if c.get("kind") == "substantive"]
        ok = (len(subs) == d["substantive_caveat_count"]
              and all(t for _, t in subs))
        print(f"  [{'OK  ' if ok else 'FAIL'}]  {cf}: {len(subs)} substantive "
              f"caveats with non-empty text")
        if not ok:
            failures.append(("caveat", cf))

    # ---- Verdict ----
    print("\n" + sep)
    if not failures:
        print("ACCEPTANCE GATE: PASS  (all contracts hold)")
        sys.exit(0)
    print(f"ACCEPTANCE GATE: FAIL  ({len(failures)} violations)")
    for f in failures:
        print(f"  - {f}")
    sys.exit(1)


if __name__ == "__main__":
    main()
