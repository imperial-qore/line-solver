"""
Markov Chain analysis algorithms.

Native Python implementations for continuous-time and discrete-time
Markov chain analysis.

Key algorithms:
    ctmc_solve: CTMC steady-state distribution
    ctmc_sens: Steady-state sensitivity to a scalar parameter
    ctmc_transient_sens: Transient sensitivity to a scalar parameter
    ctmc_transient: CTMC transient analysis
    ctmc_uniformization: Transient distribution by Jensen uniformization
    ctmc_randomization: Randomized (uniformized) DTMC P = I + Q/q
    ctmc_foxglynn: Fox-Glynn uniformization for transient analysis
    ctmc_fau: Fast adaptive uniformization for transient analysis
    ctmc_gmres: Restarted GMRES with ILUT preconditioning
    ctmc_bicgstab: BiCGSTAB with the same ILUT preconditioner
    ctmc_saddlepoint: Pr{N(t)=k} of a MAP counting process by saddlepoint
    ctmc_stochcomp: Stochastic complementation
    dtmc_solve: DTMC steady-state distribution
"""

from .ctmc import (
    ctmc_solve,
    ctmc_solve_reducible,
    ctmc_solve_reducible_blkdecomp,
    ctmc_sens,
    ctmc_transient_sens,
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

from .fau import (
    ctmc_fau,
    CtmcFauInfo,
)

from .gmres import (
    ctmc_gmres,
    ctmc_gmres_multi,
)

from .bicgstab import (
    ctmc_bicgstab,
    ctmc_bicgstab_multi,
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
    dtmc_stochcomp_full,
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

from .passage import (
    ctmc_passage_ph,
    ctmc_passage_lst,
    ctmc_passage_moments,
    ctmc_passage_time,
    ctmc_hitting_time,
    smp_passage_lst,
    smp_passage_moments,
    smp_passage_time,
)

from .saddlepoint import ctmc_saddlepoint, K2_MIN

__all__ = [
    # CTMC functions
    'ctmc_solve',
    'ctmc_solve_reducible',
    'ctmc_solve_reducible_blkdecomp',
    'ctmc_sens',
    'ctmc_transient_sens',
    'ctmc_makeinfgen',
    'ctmc_transient',
    'ctmc_timeaverage',
    'ctmc_uniformization',
    'ctmc_foxglynn',
    'ctmc_foxglynn_weights',
    'ctmc_fau',
    'CtmcFauInfo',
    'ctmc_gmres',
    'ctmc_gmres_multi',
    'ctmc_bicgstab',
    'ctmc_bicgstab_multi',
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
    'dtmc_stochcomp_full',
    'dtmc_transient',
    'dtmc_hitting_time',
    # First passage times (Harrison and Knottenbelt 2002)
    'ctmc_passage_ph',
    'ctmc_passage_lst',
    'ctmc_passage_moments',
    'ctmc_passage_time',
    'ctmc_hitting_time',
    'smp_passage_lst',
    'smp_passage_moments',
    'smp_passage_time',
    # Aggregation methods
    'CourtoisResult',
    'KMSResult',
    'TakahashiResult',
    'ctmc_courtois',
    'ctmc_kms',
    'ctmc_takahashi',
    'ctmc_multi',
    # Saddlepoint approximation of the MAP counting process (Daniels 1954)
    'ctmc_saddlepoint',
    'K2_MIN',
]
