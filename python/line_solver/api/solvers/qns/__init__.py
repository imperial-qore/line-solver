"""
QNS Solver - Native Python handler for qnsolver external tool.

This module provides a native Python implementation that calls the
external qnsolver tool via subprocess, bypassing the Java wrapper.

Requires:
    qnsolver command-line tool in system PATH
"""

from .handler import solver_qns, SolverQNSOptions, SolverQNSReturn, is_qns_available

__all__ = [
    'solver_qns',
    'SolverQNSOptions',
    'SolverQNSReturn',
    'is_qns_available',
]
