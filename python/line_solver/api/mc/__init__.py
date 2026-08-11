"""
Markov Chain analysis algorithms.

Native Python implementations for continuous-time and discrete-time
Markov chain analysis.

Key algorithms:
    ctmc_solve: CTMC steady-state distribution
    ctmc_sens: Steady-state sensitivity to a scalar parameter
    ctmc_transient: CTMC transient analysis
    ctmc_uniformization: Uniformization method for transient analysis
    ctmc_foxglynn: Fox-Glynn uniformization for transient analysis
    ctmc_gmres: Restarted GMRES with ILUT preconditioning
    ctmc_stochcomp: Stochastic complementation
    dtmc_solve: DTMC steady-state distribution
"""

from .ctmc import (
    ctmc_solve,
    ctmc_solve_reducible,
    ctmc_sens,
    ctmc_makeinfgen,
    ctmc_transient,
    ctmc_timeaverage,
    ctmc_uniformization,
    ctmc_randomization,
    ctmc_stochcomp,
    ctmc_timereverse,
    ctmc_rand,
    ctmc_simulate,
    ctmc_isfeasible,
    ctmc_ssg,
    ctmc_ssg_reachability,
    ctmc_memory_gate,
    CtmcSsgResult,
)

from .foxglynn import (
    ctmc_foxglynn,
    ctmc_foxglynn_weights,
)

from .gmres import (
    ctmc_gmres,
    ctmc_gmres_multi,
)

from .dtmc import (
    dtmc_solve,
    dtmc_solve_reducible,
    dtmc_makestochastic,
    dtmc_isfeasible,
    dtmc_simulate,
    dtmc_rand,
    dtmc_timereverse,
    dtmc_stochcomp,
    dtmc_transient,
    dtmc_hitting_time,
)

from .aggregation import (
    CourtoisResult,
    KMSResult,
    TakahashiResult,
    ctmc_courtois,
    ctmc_kms,
    ctmc_takahashi,
    ctmc_multi,
)

__all__ = [
    # CTMC functions
    'ctmc_solve',
    'ctmc_solve_reducible',
    'ctmc_sens',
    'ctmc_makeinfgen',
    'ctmc_transient',
    'ctmc_timeaverage',
    'ctmc_uniformization',
    'ctmc_foxglynn',
    'ctmc_foxglynn_weights',
    'ctmc_gmres',
    'ctmc_gmres_multi',
    'ctmc_randomization',
    'ctmc_stochcomp',
    'ctmc_timereverse',
    'ctmc_rand',
    'ctmc_simulate',
    'ctmc_isfeasible',
    'ctmc_ssg',
    'ctmc_ssg_reachability',
    'ctmc_memory_gate',
    'CtmcSsgResult',
    # DTMC functions
    'dtmc_solve',
    'dtmc_solve_reducible',
    'dtmc_makestochastic',
    'dtmc_isfeasible',
    'dtmc_simulate',
    'dtmc_rand',
    'dtmc_timereverse',
    'dtmc_stochcomp',
    'dtmc_transient',
    'dtmc_hitting_time',
    # Aggregation methods
    'CourtoisResult',
    'KMSResult',
    'TakahashiResult',
    'ctmc_courtois',
    'ctmc_kms',
    'ctmc_takahashi',
    'ctmc_multi',
]
