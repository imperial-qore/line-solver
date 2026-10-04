"""
FLD Solver handlers for native Python implementation.

Provides the core analysis handlers for the Fluid/Mean-Field Approximation (FLD) solver,
supporting both matrix method and closing/state-dependent methods.

Key handlers:
    solver_fld: Main FLD solver handler
    solver_fld_matrix: Matrix method analyzer

Analyzers:
    solver_fld_analyzer: Main analyzer with method selection
"""

from .handler import (
    solver_fld,
    SolverFLDReturn,
    SolverFLDOptions,
)

from .analyzers import (
    solver_fld_analyzer,
    FLDResult,
    matrix_method_analyzer,
    mfq_analyzer,
)

__all__ = [
    'solver_fld',
    'SolverFLDReturn',
    'SolverFLDOptions',
    'solver_fld_analyzer',
    'FLDResult',
    'matrix_method_analyzer',
    'mfq_analyzer',
]
