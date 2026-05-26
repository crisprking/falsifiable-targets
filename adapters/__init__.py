"""falsifiable-targets adapters."""
from .io import (
    FixtureAdapter,
    CompositeAdapter,
    UniProtAdapter,
    ChEMBLAdapter,
    default_composite,
)

__all__ = [
    "FixtureAdapter",
    "CompositeAdapter",
    "UniProtAdapter",
    "ChEMBLAdapter",
    "default_composite",
]
