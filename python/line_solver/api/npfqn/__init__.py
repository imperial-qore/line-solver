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

from .rqna import (
    npfqn_rqna_weight,
    npfqn_traffic_idc,
)

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
    'npfqn_traffic_idc',
]
