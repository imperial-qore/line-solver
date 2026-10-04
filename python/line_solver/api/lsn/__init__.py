"""
Layered Stochastic Network (LSN) utilities.

Native Python implementations for layered queueing
network analysis.

Key functions:
    lsn_max_multiplicity: Compute maximum multiplicity for tasks in a layered network.

Key classes:
    LayeredNetworkStruct: Structure representing a layered network.
    LayeredNetworkElement: Enumeration of layered network element types.
"""

from .max_multiplicity import (
    LayeredNetworkElement,
    LayeredNetworkStruct,
    kahn_topological_sort,
    lsn_max_multiplicity,
)

__all__ = [
    'LayeredNetworkElement',
    'LayeredNetworkStruct',
    'kahn_topological_sort',
    'lsn_max_multiplicity',
]
