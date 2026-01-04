"""
MAM Solver handlers for native Python implementation.

Provides the core analysis handlers for the Matrix-Analytic Methods (MAM) solver,
supporting analysis of queueing networks with Markovian Arrival Processes (MAP),
phase-type service distributions, and other non-product-form characteristics.

Key handlers:
    solver_mam: Main MAM solver handler
    solver_mam_basic: Basic decomposition method

Analyzers:
    solver_mam_analyzer: Main analyzer with method selection
"""

from .handler import (
    solver_mam,
    solver_mam_basic,
    SolverMAMReturn,
    SolverMAMOptions,
)

from .analyzers import (
    solver_mam_analyzer,
    MAMResult,
)

__all__ = [
    'solver_mam',
    'solver_mam_basic',
    'SolverMAMReturn',
    'SolverMAMOptions',
    'solver_mam_analyzer',
    'MAMResult',
]
