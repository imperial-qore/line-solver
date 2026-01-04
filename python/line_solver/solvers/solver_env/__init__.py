"""
SolverENV: Native Python implementation for Random Environment analysis.

Provides Numba JIT-optimized solvers for queueing networks immersed in
random environments (CTMC-modulated systems).
"""

from .env_jit import (
    HAS_NUMBA,
    exponential_cdf_jit,
    exponential_survival_jit,
    sojourn_survival_jit,
    sojourn_cdf_jit,
    transition_cdf_jit,
    compute_dtmc_transition_matrix_jit,
    simpson_integrate_survival_jit,
    compute_hold_times_jit,
    blend_stage_results_jit,
    compute_embedding_weights_jit,
    check_convergence_jit,
)

__all__ = [
    'HAS_NUMBA',
    'exponential_cdf_jit',
    'exponential_survival_jit',
    'sojourn_survival_jit',
    'sojourn_cdf_jit',
    'transition_cdf_jit',
    'compute_dtmc_transition_matrix_jit',
    'simpson_integrate_survival_jit',
    'compute_hold_times_jit',
    'blend_stage_results_jit',
    'compute_embedding_weights_jit',
    'check_convergence_jit',
]
