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
    Discrete time (slotted): qsys_geogeo1, qsys_geoxgeo1
"""

from .mapdc import qsys_mapdc, qsys_mapd1
from .mdc_crommelin import qsys_mdc_crommelin
from .dmc import qsys_dmc
from .phm1 import qsys_phm1
from .bmapm1 import qsys_bmapm1
from .phmc import qsys_phmc
from .discrete import qsys_geogeo1, qsys_geoxgeo1, qsys_geoxgeo1_moments
from .ps import qsys_mm1_ps

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
    'qsys_mapd1',
    'qsys_mdc_crommelin',
    'qsys_dmc',
    'qsys_phm1',
    'qsys_phmc',
    'qsys_bmapm1',
    # Discrete-time (slotted) queues
    'qsys_geogeo1',
    'qsys_geoxgeo1',
    'qsys_geoxgeo1_moments',
    # Multiclass processor sharing
    'qsys_mm1_ps',
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
]
