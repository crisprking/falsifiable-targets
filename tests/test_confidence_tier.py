import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import direction_audit as da

def test_confidence_tier_logic():
    assert da.confidence_tier("inhibitor", "inhibitor") == "HIGH"
    assert da.confidence_tier("activator", "inhibitor") == "CAVEAT"
    assert da.confidence_tier("inhibitor", None) == "STANDARD"
    assert da.confidence_tier("inhibitor", "conflicted") == "STANDARD"
    assert da.confidence_tier(None, "inhibitor") == "REFUSED"
