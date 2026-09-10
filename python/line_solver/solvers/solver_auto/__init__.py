"""
Native Python implementation of AUTO solver.

This module provides automatic solver selection based on model characteristics.
"""

from .solver_auto import SolverAUTO, SolverAUTOOptions, ModelAnalyzer

# Pre-rename names, kept so existing scripts and notebooks keep importing.
SolverAuto = SolverAUTO
SolverAutoOptions = SolverAUTOOptions

__all__ = [
    'SolverAUTO',
    'SolverAUTOOptions',
    'ModelAnalyzer',
    'SolverAuto',
    'SolverAutoOptions',
]
