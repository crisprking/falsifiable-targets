"""falsifiable-targets adapters."""
from .io import (
    FixtureAdapter,
    CompositeAdapter,
    UniProtAdapter,
    ChEMBLAdapter,
    default_composite,
    load_paralog_map,
)

__all__ = [
    "FixtureAdapter",
    "CompositeAdapter",
    "UniProtAdapter",
    "ChEMBLAdapter",
    "default_composite",
    "load_paralog_map",
]
