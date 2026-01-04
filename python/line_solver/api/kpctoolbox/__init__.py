"""
KPC-Toolbox: Markov Chain and Phase-Type Distribution Analysis.

Native Python implementations of the KPC-Toolbox functions for analyzing
Markov chains, phase-type distributions, and Markovian arrival processes.

Key modules:
    basic: Basic utility functions (minpos, maxpos, logspacei, spectd)
    mc: CTMC and DTMC analysis (ctmc_solve, dtmc_solve, etc.)
    aph: Acyclic Phase-Type distributions
    mmpp: Markov Modulated Poisson Processes
"""

from .basic import (
    minpos,
    maxpos,
    logspacei,
    spectd,
    ones,
    e,
    eye,
    zeros,
)

from .mc import (
    # CTMC functions
    ctmc_makeinfgen,
    ctmc_solve,
    ctmc_solve_full,
    ctmc_rand,
    ctmc_timereverse,
    ctmc_randomization,
    ctmc_uniformization,
    ctmc_relsolve,
    weaklyconncomp,
    # DTMC functions
    dtmc_makestochastic,
    dtmc_isfeasible,
    dtmc_solve,
    dtmc_rand,
    dtmc_simulate,
    dtmc_stochcomp,
    dtmc_timereverse,
    dtmc_uniformization,
)

from .aph import (
    aph_simplify,
    aph_convpara,
    aph_convseq,
    aph_rand,
    aph_fit,
    ph2hyper,
    hyper_rand,
    ConvolutionPattern,
)

from .mmpp import (
    mmpp2_fit,
    mmpp2_fit1,
    mmpp2_fit2,
    mmpp2_fit3,
    mmpp2_fit4,
    mmpp2_fitc,
    mmpp2_fitc_approx,
    mmpp2_fitc_theoretical,
    mmpp_rand,
)

from .kpcfit import (
    KPCFIT_TOL,
    KpcfitTraceData,
    KpcfitPhOptions,
    KpcfitResult,
    kpcfit_tol,
    logspacei,
    kpcfit_init,
    kpcfit_sub_eval_acfit,
    kpcfit_sub_bic,
    kpcfit_sub_compose,
    kpcfit_hyper_charpoly,
    kpcfit_ph_prony,
    kpcfit_ph_options,
    kpcfit_ph_exact,
    kpcfit_ph_auto,
)

__all__ = [
    # Basic utilities
    'minpos',
    'maxpos',
    'logspacei',
    'spectd',
    'ones',
    'e',
    'eye',
    'zeros',
    # CTMC functions
    'ctmc_makeinfgen',
    'ctmc_solve',
    'ctmc_solve_full',
    'ctmc_rand',
    'ctmc_timereverse',
    'ctmc_randomization',
    'ctmc_uniformization',
    'ctmc_relsolve',
    'weaklyconncomp',
    # DTMC functions
    'dtmc_makestochastic',
    'dtmc_isfeasible',
    'dtmc_solve',
    'dtmc_rand',
    'dtmc_simulate',
    'dtmc_stochcomp',
    'dtmc_timereverse',
    'dtmc_uniformization',
    # APH functions
    'aph_simplify',
    'aph_convpara',
    'aph_convseq',
    'aph_rand',
    'aph_fit',
    'ph2hyper',
    'hyper_rand',
    'ConvolutionPattern',
    # MMPP functions
    'mmpp2_fit',
    'mmpp2_fit1',
    'mmpp2_fit2',
    'mmpp2_fit3',
    'mmpp2_fit4',
    'mmpp2_fitc',
    'mmpp2_fitc_approx',
    'mmpp2_fitc_theoretical',
    'mmpp_rand',
    # KPC fitting functions
    'KPCFIT_TOL',
    'KpcfitTraceData',
    'KpcfitPhOptions',
    'KpcfitResult',
    'kpcfit_tol',
    'kpcfit_init',
    'kpcfit_sub_eval_acfit',
    'kpcfit_sub_bic',
    'kpcfit_sub_compose',
    'kpcfit_hyper_charpoly',
    'kpcfit_ph_prony',
    'kpcfit_ph_options',
    'kpcfit_ph_exact',
    'kpcfit_ph_auto',
]
