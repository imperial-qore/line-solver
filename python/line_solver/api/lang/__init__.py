"""
Native Python implementations for LINE lang module.

This module provides pure Python implementations of network generation
utilities that do not require the Java backend.

The generator itself lives in line_solver.gen; it is re-exported here so
that both import paths resolve to the same class.
"""

from ...gen.network_generator import (
    NetworkGenerator,
    rand_graph,
    cyclic_graph,
    rand_spanning_tree,
)

__all__ = [
    'NetworkGenerator',
    'rand_graph',
    'cyclic_graph',
    'rand_spanning_tree',
]
