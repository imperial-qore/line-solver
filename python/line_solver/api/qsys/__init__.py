"""
Native Python implementations for queueing system analysis.

This module provides pure Python/NumPy implementations for analyzing
single queueing systems, including basic queues (M/M/1, M/M/k, M/G/1),
G/G/1 approximations, MAP-based queues, and scheduling disciplines.

Key algorithms:
    Basic queues: qsys_mm1, qsys_mmk, qsys_mg1, qsys_gm1, qsys_mminf, qsys_mginf
    G/G/1 approximations: Allen-Cunneen, Kingman, Marchal, Whitt, Heyman, etc.
    G/G/k approximations: qsys_gigk_approx
    MAP/D queues: qsys_mapdc, qsys_mapd1
    MAP/PH queues: qsys_phph1, qsys_mapph1, qsys_mapm1, qsys_mapmc, qsys_mapmap1
    Scheduling: qsys_mg1_prio, qsys_mg1_srpt, qsys_mg1_fb, etc.
    Loss systems: qsys_mm1k_loss, qsys_mg1k_loss, qsys_mxm1
    Exact finite buffer: qsys_mapg1k, qsys_mmapg1k, qsys_mapg1k_perflow
    Conditional Lindley: qsys_mm1_lindley, qsys_hh1_lindley, qsys_tandem_lindley
    Tandem tail bounds: qsys_tandem_ub_ciucu
    Abandonment: qsys_mgisrgi_whitt (M/GI/s/r+GI), qsys_erlanga (M/M/s/r+M),
    qsys_ggisgi_fluid (G/GI/s+GI fluid limit)
    QED regime: qsys_mmk_qed, qsys_mmk_qed_alpha, qsys_mmk_qed_staffing
    Time-varying: qsys_mtginf (Mt/G/inf, exact)
    Extremal bounds: qsys_gig1_bnds_extremal
    Time-varying fluid: qsys_gtmtst_fluid (Gt/Mt/st+GI)
    Diffusion: qsys_ggnm_diffusion (G/GI/n/m), qsys_ggingi_tga (G/GI/n+GI)
"""

from .mapdc import qsys_mapdc, qsys_mapd1
from .mapphc import qsys_mapphc, MapPhcResult
from .mmapgk1 import qsys_mmapgk1, MmapGk1Result
from .mdc_crommelin import qsys_mdc_crommelin
from .dmc import qsys_dmc
from .phm1 import qsys_phm1
from .bmapm1 import qsys_bmapm1
from .phmc import qsys_phmc
from .ps import qsys_mm1_ps
from .mg1ps import qsys_mg1_ps
from .lindley import (
    qsys_lindley_moment,
    qsys_mm1_lindley,
    qsys_hh1_lindley,
    qsys_mm1_tandem_lindley,
    qsys_tandem_lindley,
)

from .tandem_bounds import qsys_tandem_ub_ciucu

from .abandonment import (
    qsys_mgisrgi_whitt,
    qsys_erlanga,
)

from .fluid_abandonment import qsys_ggisgi_fluid

from .mtginf import qsys_mtginf

from .extremal import qsys_gig1_bnds_extremal

from .tvfluid import qsys_gtmtst_fluid

from .diffusion import qsys_ggnm_diffusion

from .gaussian_ed import qsys_ggingi_tga

from .mol import qsys_mtgs0_mol, erlang_b, erlang_c

from .maxima import qsys_maxima_twomoment

from .qed import (
    qsys_mmk_qed,
    qsys_mmk_qed_alpha,
    qsys_mmk_qed_staffing,
)

from .rqt import (
    qsys_gigk_rqt,
    qsys_gig1_rqt,
    qsys_gigk_rqt_gamma,
)

from .basic import (
    qsys_mm1,
    qsys_mmk,
    qsys_mmck,
    qsys_mg1,
    qsys_gm1,
    qsys_mminf,
    qsys_mginf,
    qsys_mmcc_retrial_fp,
    qsys_gig1_rq,
)

from .approximations import (
    qsys_gig1_approx_allencunneen,
    qsys_gig1_approx_kingman,
    qsys_gig1_approx_marchal,
    qsys_gig1_approx_whitt,
    qsys_gig1_approx_heyman,
    qsys_gig1_approx_kobayashi,
    qsys_gig1_approx_klb,
    qsys_gig1_approx_gelenbe,
    qsys_gig1_approx_kimura,
    qsys_gig1_approx_myskja,
    qsys_gig1_approx_myskja2,
    qsys_gigk_approx,
    qsys_gig1_ubnd_kingman,
    qsys_gigk_approx_kingman,
    qsys_gg1,
    qsys_gig1_lbnd,
    qsys_gigk_approx_cosmetatos,
    qsys_gigk_approx_whitt,
)

from .scheduling import (
    qsys_mg1_prio,
    qsys_mm1_dps,
    qsys_mg1_srpt,
    qsys_mg1_fb,
    qsys_mg1_lrpt,
    qsys_mg1_psjf,
    qsys_mg1_setf,
)

from .loss import (
    qsys_mm1k_loss,
    qsys_mg1k_loss,
    qsys_mg1k_loss_mgs,
    qsys_mxm1,
)

from .workload import (
    qsys_ldps_workload,
)

from .mapg1k import (
    qsys_mapg1k,
    qsys_mmapg1k,
    qsys_mapg1k_perflow,
)

from .map_queues import (
    QueueResult,
    ph_to_map,
    qsys_phph1,
    qsys_mapph1,
    qsys_mapm1,
    qsys_mapmc,
    qsys_mapmap1,
    qsys_mapg1,
)

from .retrial import (
    QueueType,
    BmapMatrix,
    PhDistribution,
    QbdStatespace,
    RetrialQueueResult,
    RetrialQueueAnalyzer,
    qsys_bmapphnn_retrial,
    qsys_is_retrial,
    RetrialInfo,
    RenegingInfo,
    detect_reneging_topology,
    has_reneging_patience,
    extract_bmap_matrices,
    extract_ph_params,
    convert_patience_to_regimes,
    solver_mam_retrial,
)

__all__ = [
    # MAP/D queues
    'qsys_mapdc',
    'qsys_mapphc',
    'qsys_mmapgk1',
    'MmapGk1Result',
    'MapPhcResult',
    'qsys_mapd1',
    'qsys_mdc_crommelin',
    'qsys_dmc',
    'qsys_phm1',
    'qsys_phmc',
    'qsys_bmapm1',
    # Discrete-time (slotted) queues
    # Multiclass processor sharing
    'qsys_mm1_ps',
    'qsys_mg1_ps',
    # Basic queues
    'qsys_mm1',
    'qsys_mmk',
    'qsys_mmck',
    'qsys_mg1',
    'qsys_gm1',
    'qsys_mminf',
    'qsys_mginf',
    'qsys_mmcc_retrial_fp',
    # G/G/1 approximations
    'qsys_gig1_rq',
    'qsys_gigk_rqt',
    'qsys_gig1_rqt',
    'qsys_gigk_rqt_gamma',
    'qsys_gig1_approx_allencunneen',
    'qsys_gig1_approx_kingman',
    'qsys_gig1_approx_marchal',
    'qsys_gig1_approx_whitt',
    'qsys_gig1_approx_heyman',
    'qsys_gig1_approx_kobayashi',
    'qsys_gig1_approx_klb',
    'qsys_gig1_approx_gelenbe',
    'qsys_gig1_approx_kimura',
    'qsys_gig1_approx_myskja',
    'qsys_gig1_approx_myskja2',
    # G/G/1 analysis
    'qsys_gg1',
    # G/G/1 lower bound
    'qsys_gig1_lbnd',
    # G/G/k approximations
    'qsys_gigk_approx',
    'qsys_gigk_approx_cosmetatos',
    'qsys_gigk_approx_whitt',
    # Upper bounds and Kingman multi-server
    'qsys_gig1_ubnd_kingman',
    'qsys_gigk_approx_kingman',
    # Scheduling disciplines
    'qsys_mg1_prio',
    'qsys_mm1_dps',
    'qsys_mg1_srpt',
    'qsys_mg1_fb',
    'qsys_mg1_lrpt',
    'qsys_mg1_psjf',
    'qsys_mg1_setf',
    # Loss systems
    'qsys_mm1k_loss',
    'qsys_mg1k_loss',
    'qsys_mg1k_loss_mgs',
    'qsys_mxm1',
    # Exact MAP/G/1/K finite-buffer family
    'qsys_mapg1k',
    'qsys_mmapg1k',
    'qsys_mapg1k_perflow',
    # MAP/PH queues
    'QueueResult',
    'ph_to_map',
    'qsys_phph1',
    'qsys_mapph1',
    'qsys_mapm1',
    'qsys_mapmc',
    'qsys_mapmap1',
    'qsys_mapg1',
    # Retrial queueing framework
    'QueueType',
    'BmapMatrix',
    'PhDistribution',
    'QbdStatespace',
    'RetrialQueueResult',
    'RetrialQueueAnalyzer',
    'qsys_bmapphnn_retrial',
    'qsys_is_retrial',
    'RetrialInfo',
    'RenegingInfo',
    'detect_reneging_topology',
    'has_reneging_patience',
    'extract_bmap_matrices',
    'extract_ph_params',
    'convert_patience_to_regimes',
    'solver_mam_retrial',
    # Workload distribution
    'qsys_ldps_workload',
    # Conditional Lindley recursion
    'qsys_lindley_moment',
    'qsys_mm1_lindley',
    'qsys_hh1_lindley',
    'qsys_mm1_tandem_lindley',
    'qsys_tandem_lindley',
    'qsys_tandem_ub_ciucu',
    # Multiserver queues with customer abandonment
    'qsys_mgisrgi_whitt',
    'qsys_erlanga',
    'qsys_ggisgi_fluid',
    # Halfin-Whitt QED regime and square-root staffing
    'qsys_mmk_qed',
    'qsys_mmk_qed_alpha',
    'qsys_mmk_qed_staffing',
    # Time-varying infinite-server queue
    'qsys_mtginf',
    # Extremal two-moment bounds
    'qsys_gig1_bnds_extremal',
    # Time-varying many-server fluid queue
    'qsys_gtmtst_fluid',
    # G/GI/n/m diffusion approximation
    'qsys_ggnm_diffusion',
    # Heavily-loaded G/GI/n+GI Gaussian approximation
    'qsys_ggingi_tga',
    # Modified offered load for time-varying systems
    'qsys_mtgs0_mol',
    'erlang_b',
    'erlang_c',
    # Two-moment approximation for maxima
    'qsys_maxima_twomoment',
]
