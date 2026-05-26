"""
Run an audit on a claim YAML using the live + fixture composite adapter.

Usage:
    python -m run_audit claims/tyk2_psoriasis.yaml
    python -m run_audit claims/tyk2_psoriasis.yaml --offline

The --offline flag is for hermetic testing: forces the live adapters to
serve from .ae_cache/ only, matching what AE_OFFLINE=1 would do globally.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml


def _ruleset_sha(rules, version):
    descriptor = {
        "version": version,
        "rules": [
            {"rule_id": r.rule_id, "version": r.version, "description": r.description}
            for r in sorted(rules, key=lambda x: x.rule_id)
        ],
    }
    canon = json.dumps(descriptor, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canon.encode()).hexdigest()


def _claim_sha(claim):
    canon = json.dumps(
        {
            "target_symbol": claim.target_symbol,
            "uniprot_id": claim.uniprot_id,
            "indication": claim.indication,
            "mechanism": claim.mechanism,
            "claim_type": claim.claim_type.value,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(canon.encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("claim_path", help="Path to claim YAML")
    parser.add_argument(
        "--offline", action="store_true",
        help="Force adapters to serve from cache only (hermetic mode)",
    )
    parser.add_argument(
        "--no-live", action="store_true",
        help="Skip live adapters entirely (fixture only)",
    )
    parser.add_argument(
        "--json-out", default=None,
        help="Write structured JSON report to this path",
    )
    args = parser.parse_args()

    if args.offline:
        os.environ["AE_OFFLINE"] = "1"

    # Local imports so the env var is set before adapter import
    from smoke_test import TargetClaim, ClaimType, RULES, RuleStatus, run_audit, aggregate
    from adapters import default_composite, FixtureAdapter

    claim_path = Path(args.claim_path)
    with open(claim_path) as f:
        spec = yaml.safe_load(f)

    c = spec["claim"]
    claim = TargetClaim(
        target_symbol=c["target_symbol"],
        indication=c["indication"],
        mechanism=c["mechanism"],
        claim_type=ClaimType(c["claim_type"]),
        uniprot_id=c.get("uniprot_id"),
    )

    fixture = spec.get("fixture", {})

    # Build adapter
    if args.no_live:
        adapter = FixtureAdapter(fixture)
        adapter_mode = "fixture-only"
    else:
        adapter = default_composite(fixture, use_live=True)
        adapter_mode = "live+fixture composite" + (" (offline cache)" if args.offline else "")

    # Walk rules manually so we can capture per-rule input data
    results = []
    rule_inputs = {}
    for rule in RULES:
        if not rule.applies_to(claim):
            from smoke_test import RuleResult
            results.append(RuleResult(rule.rule_id, RuleStatus.NOT_APPLICABLE, 0.0))
            continue
        # Snapshot the data each rule sees, for the report
        if rule.rule_id == "R1_orthology":
            rule_inputs[rule.rule_id] = adapter.get("orthology", claim)
        elif rule.rule_id in ("R2_chemistry_support", "R6_chemistry_class_collapse"):
            rule_inputs[rule.rule_id] = adapter.get("chemistry", claim)
        elif rule.rule_id == "R3_genetics_support":
            rule_inputs[rule.rule_id] = adapter.get("genetics", claim)
        elif rule.rule_id == "R4_expression":
            rule_inputs[rule.rule_id] = adapter.get("expression", claim)
        elif rule.rule_id == "R5_replication":
            rule_inputs[rule.rule_id] = adapter.get("reproducibility", claim)
        elif rule.rule_id == "R7_selectivity_counterscreen":
            rule_inputs[rule.rule_id] = adapter.get("selectivity", claim)
        results.append(rule.evaluate(claim, adapter))

    verdict, score, cheapest, substantive, operational = aggregate(results)

    # Compute deterministic stamps
    ruleset_sha = _ruleset_sha(RULES, "1.1.0")
    claim_sha = _claim_sha(claim)
    audit_ts = datetime.now(timezone.utc).isoformat()

    # Console output
    print("=" * 72)
    print(f"Audit: {claim.target_symbol}  [{claim.indication}]")
    print(f"Claim type: {claim.claim_type.value}")
    print(f"Adapter:    {adapter_mode}")
    print(f"Ruleset:    v1.1.0  ({ruleset_sha[:16]}...)")
    print(f"Claim SHA:  {claim_sha[:16]}...")
    print(f"Timestamp:  {audit_ts}")
    print("=" * 72)
    print(f"\nVERDICT: {verdict.value}")
    print(f"Score:   {score:.2f}  ({len([r for r in results if r.status == RuleStatus.FALSIFIED])} falsified / {len([r for r in results if r.status != RuleStatus.NOT_APPLICABLE])} applicable)")
    if cheapest:
        print(f"Cheapest falsification: {cheapest[2]} @ {cheapest[1].value}")
        print(f"  Experiment: {cheapest[0]}")

    print(f"\nPer-rule results:")
    for r in results:
        status_str = r.status.value.upper().ljust(15)
        print(f"  {r.rule_id:35} {status_str} conf={r.confidence}")
        if r.rule_id in rule_inputs and rule_inputs[r.rule_id]:
            data = rule_inputs[r.rule_id]
            data_str = ", ".join(f"{k}={v!r}" for k, v in list(data.items())[:4])
            if len(data) > 4:
                data_str += f", ... +{len(data) - 4} more"
            print(f"  {'':35}   input: {data_str}")

    print(f"\nSubstantive caveats ({len(substantive)}):")
    if not substantive:
        print("  (none)")
    for cav in substantive:
        print(f"  - [{cav.rule_id}] {cav.text}")

    print(f"\nOperational notes ({len(operational)}):")
    if not operational:
        print("  (none)")
    for cav in operational:
        print(f"  - [{cav.rule_id}] {cav.text}")

    # JSON report
    report = {
        "schema_version": "1.0",
        "audit_timestamp_utc": audit_ts,
        "ruleset_version": "1.1.0",
        "ruleset_sha256": ruleset_sha,
        "claim": {
            "target_symbol": claim.target_symbol,
            "uniprot_id": claim.uniprot_id,
            "indication": claim.indication,
            "claim_type": claim.claim_type.value,
            "claim_sha256": claim_sha,
        },
        "verdict": verdict.value,
        "score": score,
        "cheapest_falsification": (
            {
                "rule_id": cheapest[2],
                "tier": cheapest[1].value,
                "experiment": cheapest[0],
            } if cheapest else None
        ),
        "per_rule": [
            {
                "rule_id": r.rule_id,
                "status": r.status.value,
                "confidence": r.confidence,
                "input_data": rule_inputs.get(r.rule_id, {}),
                "evidence": r.evidence,
                "falsification_tier": r.falsification_tier.value if r.falsification_tier else None,
                "falsification_experiment": r.falsification_experiment,
                "caveats": [
                    {"kind": c.kind.value, "text": c.text}
                    for c in r.caveats
                ],
            }
            for r in results
        ],
        "substantive_caveat_count": len(substantive),
        "operational_note_count": len(operational),
        "adapter_mode": adapter_mode,
    }

    if args.json_out:
        Path(args.json_out).write_text(json.dumps(report, indent=2))
        print(f"\nJSON report written to: {args.json_out}")

    # Exit code reflects verdict for CI use:
    # 0 = SURVIVED, 1 = FALSIFIED_WITH_CAVEATS, 2 = FALSIFIED, 3 = INSUFFICIENT_DATA
    from smoke_test import Verdict
    exit_codes = {
        Verdict.SURVIVED: 0,
        Verdict.FALSIFIED_WITH_CAVEATS: 1,
        Verdict.FALSIFIED: 2,
        Verdict.INSUFFICIENT_DATA: 3,
    }
    sys.exit(exit_codes.get(verdict, 99))


if __name__ == "__main__":
    main()
