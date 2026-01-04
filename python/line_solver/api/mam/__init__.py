"""
Matrix-Analytic Methods (MAM) for MAP/PH distributions.

Native Python implementations for analyzing Markovian Arrival Processes (MAPs),
Phase-Type (PH) distributions, and related matrix-analytic methods.

Key algorithms:
    map_piq: CTMC steady-state of MAP
    map_pie: Embedded DTMC steady-state
    map_lambda: Arrival rate computation
    map_mean, map_var, map_scv: Moment computations
    solver_mam_map_bmap_1: MAP/BMAP/1 queue solver using GI/M/1-type ETAQA
    solver_mam_bmap_map_1: BMAP/MAP/1 queue solver using M/G/1-type ETAQA
"""

from .map_analysis import (
    map_infgen,
    map_piq,
    map_pie,
    map_lambda,
    map_mean,
    map_var,
    map_scv,
    map_moment,
    map_scale,
    map_normalize,
    map_isfeasible,
    # Rate-based constructors
    exp_map,
    erlang_map,
    hyperexp_map,
    # MATLAB-style mean-based constructors
    map_exponential,
    map_erlang,
    map_hyperexp,
    map_gamma,
    # MAP operations
    map_sumind,
    map_cdf,
    map_pdf,
    map_sample,
    # Statistics functions
    map_skew,
    map_kurt,
    map_acf,
    map_acfc,
    map_idc,
    # Count process functions
    map_count_mean,
    map_count_var,
    map_varcount,
    map_count_moment,
    # Constructor functions
    map_mmpp2,
    map_gamma2,
    map_rand,
    map_randn,
    map_renewal,
    # Operation functions
    map_embedded,
    map_sum,
    map_super,
    map_mixture,
    map_max,
    map_timereverse,
    map_mark,
    map_stochcomp,
    # Fitting functions
    map_kpc,
    map_bernstein,
    map_pntiter,
    map_pntquad,
    map2_fit,
    # Utility functions
    map_joint,
    map_issym,
    map_feastol,
    map_feasblock,
    map_block,
    map_largemap,
)

from .mmap_ops import (
    mmap_infgen,
    mmap_normalize,
    mmap_super_safe,
    mmap_mark,
    mmap_scale,
    mmap_hide,
    mmap_compress,
    mmap_exponential,
    # New MMAP functions
    mmap_issym,
    mmap_isfeasible,
    mmap_lambda,
    mmap_count_lambda,
    mmap_pie,
    mmap_pc,
    mmap_embedded,
    mmap_sample,
    mmap_rand,
    mmap_timereverse,
    mmap_sum,
    mmap_super,
    mmap_mixture,
    mmap_max,
    mmap_maps,
    mmap_count_mean,
    mmap_count_var,
    mmap_count_idc,
    mmap_count_mcov,
    mmap_idc,
    mmap_sigma,
    mmap_sigma2,
    mmap_forward_moment,
    mmap_backward_moment,
    mmap_cross_moment,
    mmap_modulate,
)

from .ldqbd import (
    ldqbd,
    LdqbdResult,
    LdqbdOptions,
)

from .map_bmap_1 import (
    MAPBMAP1Result,
    solver_mam_map_bmap_1,
)

from .bmap_map_1 import (
    BMAPMAP1Result,
    solver_mam_bmap_map_1,
)

from .qbd import (
    QBDResult,
    qbd_R,
    qbd_R_logred,
    qbd_rg,
    qbd_blocks_mapmap1,
    qbd_bmapbmap1,
    qbd_mapmap1,
    qbd_raprap1,
    qbd_setupdelayoff,
)

from .map_derivatives import (
    map_ccdf_derivative,
    map_jointpdf_derivative,
    map_factorial_moment,
    map_joint_moment,
)

from .mapm1ps import (
    map_m1ps_cdf_respt,
    map_compute_R,
    map_m1ps_h_recursive,
    map_m1ps_sojourn,
)

from ..mmdp import mmdp_isfeasible

__all__ = [
    # MAP analysis
    'map_infgen',
    'map_piq',
    'map_pie',
    'map_lambda',
    'map_mean',
    'map_var',
    'map_scv',
    'map_moment',
    'map_scale',
    'map_normalize',
    'map_isfeasible',
    # MAP constructors (rate-based)
    'exp_map',
    'erlang_map',
    'hyperexp_map',
    # MAP constructors (MATLAB-style mean-based)
    'map_exponential',
    'map_erlang',
    'map_hyperexp',
    'map_gamma',
    # MAP operations
    'map_sumind',
    'map_cdf',
    'map_pdf',
    'map_sample',
    # Statistics functions
    'map_skew',
    'map_kurt',
    'map_acf',
    'map_acfc',
    'map_idc',
    # Count process functions
    'map_count_mean',
    'map_count_var',
    'map_varcount',
    'map_count_moment',
    # Constructor functions
    'map_mmpp2',
    'map_gamma2',
    'map_rand',
    'map_randn',
    'map_renewal',
    # Operation functions
    'map_embedded',
    'map_sum',
    'map_super',
    'map_mixture',
    'map_max',
    'map_timereverse',
    'map_mark',
    'map_stochcomp',
    # Fitting functions
    'map_kpc',
    'map_bernstein',
    'map_pntiter',
    'map_pntquad',
    'map2_fit',
    # Utility functions
    'map_joint',
    'map_issym',
    'map_feastol',
    'map_feasblock',
    'map_block',
    'map_largemap',
    # MMAP operations (marked arrivals)
    'mmap_infgen',
    'mmap_normalize',
    'mmap_super_safe',
    'mmap_mark',
    'mmap_scale',
    'mmap_hide',
    'mmap_compress',
    'mmap_exponential',
    # New MMAP functions
    'mmap_issym',
    'mmap_isfeasible',
    'mmap_lambda',
    'mmap_count_lambda',
    'mmap_pie',
    'mmap_pc',
    'mmap_embedded',
    'mmap_sample',
    'mmap_rand',
    'mmap_timereverse',
    'mmap_sum',
    'mmap_super',
    'mmap_mixture',
    'mmap_max',
    'mmap_maps',
    'mmap_count_mean',
    'mmap_count_var',
    'mmap_count_idc',
    'mmap_count_mcov',
    'mmap_idc',
    'mmap_sigma',
    'mmap_sigma2',
    'mmap_forward_moment',
    'mmap_backward_moment',
    'mmap_cross_moment',
    'mmap_modulate',
    # LDQBD solver
    'ldqbd',
    'LdqbdResult',
    'LdqbdOptions',
    # MAP/BMAP/1 queue solver (GI/M/1-type)
    'MAPBMAP1Result',
    'solver_mam_map_bmap_1',
    # BMAP/MAP/1 queue solver (M/G/1-type)
    'BMAPMAP1Result',
    'solver_mam_bmap_map_1',
    # QBD utilities
    'QBDResult',
    'qbd_R',
    'qbd_R_logred',
    'qbd_rg',
    'qbd_blocks_mapmap1',
    'qbd_bmapbmap1',
    'qbd_mapmap1',
    'qbd_raprap1',
    'qbd_setupdelayoff',
    # MAP derivatives
    'map_ccdf_derivative',
    'map_jointpdf_derivative',
    'map_factorial_moment',
    'map_joint_moment',
    # MAP/M/1-PS sojourn time functions
    'map_m1ps_cdf_respt',
    'map_compute_R',
    'map_m1ps_h_recursive',
    'map_m1ps_sojourn',
    # MMDP
    'mmdp_isfeasible',
]
