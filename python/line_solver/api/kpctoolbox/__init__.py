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
    ctmc_transient,
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

from .mvph import (
    mvph_joint,
    mvph_mean_x,
    mvph_mean_y,
    mvph_cov,
    mvph_corr,
)

from .erchmm import (
    erchmm_emfit,
)

from .trace import (
    trace_mean,
    trace_var,
    trace_scv,
    trace_acf,
    trace_skew,
    trace_joint,
    trace_bicov,
    trace_idi,
    trace_idc,
    trace_gamma,
    trace_shuffle,
    trace_iat2counts,
    trace_iat2bins,
    trace_pmf,
    trace_summary,
    mtrace_mean,
    autocov,
)

from .det import (
    det_embedded,
    det_moment,
    det_scv,
    det_acf,
    det_sample,
    det_sum,
)

from .map import (
    map2ph,
    map_mmpp2,
    map_bernstein,
    me_sample,
    rap_sample,
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
    'ctmc_transient',
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
    # Trace analysis functions
    'trace_mean',
    'trace_var',
    'trace_scv',
    'trace_acf',
    'trace_skew',
    'trace_joint',
    'trace_bicov',
    'trace_idi',
    'trace_idc',
    'trace_gamma',
    'trace_shuffle',
    'trace_iat2counts',
    'trace_iat2bins',
    'trace_pmf',
    'trace_summary',
    'mtrace_mean',
    'autocov',
    # MVPH (Multivariate Phase-Type) functions
    'mvph_joint',
    'mvph_mean_x',
    'mvph_mean_y',
    'mvph_cov',
    'mvph_corr',
    # ER-CHMM functions
    'erchmm_emfit',
    # DET (deterministic process) functions
    'det_embedded',
    'det_moment',
    'det_scv',
    'det_acf',
    'det_sample',
    'det_sum',
    # MAP functions
    'map2ph',
    'map_mmpp2',
    'map_bernstein',
    'me_sample',
    'rap_sample',
]
