"""
Non-Product-Form Queueing Network (NPFQN) algorithms.

Native Python implementations for approximating performance
of non-product-form queueing networks.

Key algorithms:
    npfqn_nonexp_approx: Non-exponential distribution approximation
    npfqn_traffic_merge: Merge multiple MMAP traffic flows
    npfqn_traffic_merge_cs: Merge traffic flows with class switching
    npfqn_traffic_split_cs: Split traffic flows with class switching
"""

from .nonexp import (
    npfqn_nonexp_approx,
    NpfqnNonexpApproxResult,
)

from .traffic import (
    npfqn_traffic_merge,
    npfqn_traffic_merge_cs,
    npfqn_traffic_split_cs,
)

from .sqd import (
    npfqn_sqd,
    NpfqnSqdResult,
)

from .rqt import npfqn_traffic_rqt
from .bpt import npfqn_bnd_bpt, NpfqnBndBptResult
from .bgt import npfqn_bnd_bgt, NpfqnBndBgtResult
from .rqna import (
    npfqn_rqna_weight,
    npfqn_traffic_idc,
)

from .split_rr import (
    npfqn_traffic_split_rr,
)

from ..qsys.tvfluid import npfqn_gtmtst_fluid
from .feedback_elim import npfqn_feedback_elim
from .dps_morrison import npfqn_dps_morrison, NpfqnDpsMorrisonResult

__all__ = [
    # Non-exponential approximation
    'npfqn_nonexp_approx',
    'NpfqnNonexpApproxResult',
    # Traffic operations
    'npfqn_traffic_merge',
    'npfqn_traffic_merge_cs',
    'npfqn_traffic_split_cs',
    # Blocking-after-service approximation
    'npfqn_sqd',
    'NpfqnSqdResult',
    # Robust queueing network analyzer (RQNA) traffic equations
    'npfqn_rqna_weight',
    'npfqn_traffic_rqt',
    # Achievable-region LP relaxation (Bertsimas-Paschalidis-Tsitsiklis 1994)
    'npfqn_bnd_bpt',
    'NpfqnBndBptResult',
    # Piecewise-linear Lyapunov bound (Bertsimas-Gamarnik-Tsitsiklis 2001)
    'npfqn_bnd_bgt',
    'NpfqnBndBgtResult',
    'npfqn_traffic_idc',
    # Deterministic (round-robin) traffic split degrees
    'npfqn_traffic_split_rr',
    # Time-varying network of many-server fluid queues
    'npfqn_gtmtst_fluid',
    # Near-immediate feedback elimination for RQNA
    'npfqn_feedback_elim',
    # Morrison's heavy-usage expansion for a closed think+DPS network
    'npfqn_dps_morrison',
    'NpfqnDpsMorrisonResult',
]
