"""
Native Python implementation of AUTO solver.

This module provides automatic solver selection based on model characteristics.
"""

from .solver_auto import SolverAuto, SolverAutoOptions, ModelAnalyzer

__all__ = [
    'SolverAuto',
    'SolverAutoOptions',
    'ModelAnalyzer',
]
