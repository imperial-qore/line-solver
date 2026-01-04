"""
Native Python implementations for LINE lang module.

This module provides pure Python implementations of network generation
utilities that do not require the Java backend.
"""

from .network_generator import (
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
