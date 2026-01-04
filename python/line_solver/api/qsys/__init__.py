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
"""

from .mapdc import qsys_mapdc, qsys_mapd1

from .basic import (
    qsys_mm1,
    qsys_mmk,
    qsys_mg1,
    qsys_gm1,
    qsys_mminf,
    qsys_mginf,
)

from .approximations import (
    qsys_gig1_approx_allencunneen,
    qsys_gig1_approx_kingman,
    qsys_gig1_approx_marchal,
    qsys_gig1_approx_whitt,
    qsys_gig1_approx_heyman,
    qsys_gig1_approx_kobayashi,
    qsys_gig1_approx_gelenbe,
    qsys_gig1_approx_kimura,
    qsys_gigk_approx,
    qsys_gig1_ubnd_kingman,
    qsys_gigk_approx_kingman,
)

from .scheduling import (
    qsys_mg1_prio,
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
)

__all__ = [
    # MAP/D queues
    'qsys_mapdc',
    'qsys_mapd1',
    # Basic queues
    'qsys_mm1',
    'qsys_mmk',
    'qsys_mg1',
    'qsys_gm1',
    'qsys_mminf',
    'qsys_mginf',
    # G/G/1 approximations
    'qsys_gig1_approx_allencunneen',
    'qsys_gig1_approx_kingman',
    'qsys_gig1_approx_marchal',
    'qsys_gig1_approx_whitt',
    'qsys_gig1_approx_heyman',
    'qsys_gig1_approx_kobayashi',
    'qsys_gig1_approx_gelenbe',
    'qsys_gig1_approx_kimura',
    # G/G/k approximations
    'qsys_gigk_approx',
    # Upper bounds and Kingman multi-server
    'qsys_gig1_ubnd_kingman',
    'qsys_gigk_approx_kingman',
    # Scheduling disciplines
    'qsys_mg1_prio',
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
]
