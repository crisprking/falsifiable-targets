import os
import sys

try:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
except NameError:
    pass
import pytest

import direction_audit as da


@pytest.mark.live
def test_live_smoke_pcsk9_recovers_loci():
    """Engine must pull >=1 decided locus for PCSK9/hypercholesterolemia off the live API.
    Guards against a stubbed data layer (returns 0) passing CI. Marked live; skipped in CI."""
    if os.environ.get("SKIP_LIVE"):
        pytest.skip("SKIP_LIVE set")
    rec = da.recover("ENSG00000169174", "HP_0003124")
    assert rec["n_loci_decided"] >= 1, f"engine returned 0 decided loci: {rec}"
    assert rec["recovered"] == "inhibitor", f"PCSK9 should recover inhibitor, got {rec['recovered']}"
