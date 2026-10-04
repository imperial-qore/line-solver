"""
Native Python implementation of SSA (Stochastic Simulation Algorithm) solver.

This module provides discrete event simulation using Gillespie's algorithm.
"""

from .solver_ssa import SolverSSA, SolverSSAOptions, SamplePath, SampleEvent

__all__ = [
    'SolverSSA',
    'SolverSSAOptions',
    'SamplePath',
    'SampleEvent',
]
