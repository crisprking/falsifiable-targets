"""Version stamp for falsifiable-targets.

Used by:
  - pyproject.toml (dynamic versioning)
  - smoke_test.py (banner)
  - run_audit.py (audit report metadata)
  - tests/ (sanity checks)

When releasing, update this file ONLY. CI verifies the git tag matches.
"""

__version__ = "1.4.0"

# Ruleset version is INDEPENDENT of tool version. Same rules can ship in
# multiple tool releases (e.g., v1.4.0 packaging release keeps the
# v1.2.0 ruleset). Update RULESET_VERSION only when smoke_test.RULES
# changes in a way that bumps the SHA.
RULESET_VERSION = "1.2.0"
RULESET_SHA256 = "35ef2b2ab5363298097962a0b6ae52c70d551a1edddc341054f75cb6e4fb7221"
