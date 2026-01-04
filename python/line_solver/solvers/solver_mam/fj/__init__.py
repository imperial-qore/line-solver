"""
Fork-Join topology analysis for SolverMAM.

This module provides utilities for detecting and analyzing Fork-Join queueing
networks, extracting parameters, and computing percentile response times using
the FJ_codes toolkit.
"""

from .validator import (
    fj_isfj,
    FJValidationResult,
)

from .extractor import (
    fj_extract_params,
    FJParams,
)

__all__ = [
    'fj_isfj',
    'FJValidationResult',
    'fj_extract_params',
    'FJParams',
]
