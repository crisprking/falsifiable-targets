"""falsifiable-targets adapters."""
from .io import (
    ChEMBLAdapter,
    CompositeAdapter,
    FixtureAdapter,
    UniProtAdapter,
    default_composite,
    load_paralog_map,
)
from .opentargets import OpenTargetsAdapter
from .protocol import SECTIONS, Adapter

__all__ = [
    "Adapter",
    "SECTIONS",
    "FixtureAdapter",
    "CompositeAdapter",
    "UniProtAdapter",
    "ChEMBLAdapter",
    "OpenTargetsAdapter",
    "default_composite",
    "load_paralog_map",
]
