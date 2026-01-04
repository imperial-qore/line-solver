"""
Polling System Analysis Algorithms.

Native Python implementations for analyzing polling/vacation
queue systems with various disciplines.

Key algorithms:
    polling_qsys_exhaustive: Exhaustive polling discipline
    polling_qsys_gated: Gated polling discipline
    polling_qsys_1limited: 1-Limited polling discipline
"""

from .polling import (
    polling_qsys_exhaustive,
    polling_qsys_gated,
    polling_qsys_1limited,
)

__all__ = [
    'polling_qsys_exhaustive',
    'polling_qsys_gated',
    'polling_qsys_1limited',
]
