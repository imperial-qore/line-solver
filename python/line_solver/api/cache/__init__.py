"""
Cache Analysis Algorithms.

Native Python implementations for analyzing cache systems, including
exact recursive methods, singular perturbation methods, importance
sampling, TTL-based caches, and various approximation techniques.

Key algorithms:
    cache_erec: Exact recursive normalizing constant
    cache_prob_erec: Exact recursive state probabilities
    cache_spm: Singular perturbation method
    cache_spm_size: Ray (WKB) asymptotics of the cost-capped constant
    cache_miss: Miss rate computation
    cache_is: Importance sampling
    cache_ttl_*: TTL-based cache analysis
    cache_rrm_*: Random replacement model
"""

from .erec import (
    cache_erec,
    cache_erec_aux,
    cache_prob_erec,
    cache_cost,
    cache_cost_pathcheck,
    cache_mva,
)

from .spm_size import (
    cache_spm_size,
)

from .spm import (
    cache_xi_iter,
    cache_spm,
    cache_prob_spm,
    cache_prob_fpi,
)

from .miss import (
    cache_miss,
    cache_xi_fp,
    cache_miss_fpi,
    cache_miss_spm,
    cache_mva_miss,
    cache_miss_asy,
)

from .sampling import (
    cache_is,
    cache_prob_is,
    cache_miss_is,
    logmeanexp,
)

from .ttl import (
    cache_lrum_map_levelstats,
    cache_t_lrum_map,
    cache_ttl_lrum_map,
    cache_t_hlru,
    cache_ttl_hlru,
    cache_ttl_lrua,
)

from .rmf import (
    cache_miss_rmf,
    cache_rmf_lna,
)
from .rmf_sfifo import (
    cache_miss_sfifo_rmf,
)
from .rmf_fifo import (
    cache_miss_fifo_rmf,
)

from .rrm import (
    cache_rrm_meanfield_ode,
    cache_rrm_meanfield,
    cache_gamma_lp,
    cache_gamma,
)

__all__ = [
    'cache_lrum_map_levelstats',
    'cache_t_lrum_map',
    'cache_ttl_lrum_map',
    # Exact recursive methods
    'cache_erec',
    'cache_erec_aux',
    'cache_prob_erec',
    'cache_cost',
    'cache_cost_pathcheck',
    'cache_mva',
    # Ray (WKB) asymptotics of the cost-capped constant
    'cache_spm_size',
    # SPM methods
    'cache_xi_iter',
    'cache_spm',
    'cache_prob_spm',
    'cache_prob_fpi',
    # Miss rate methods
    'cache_miss',
    'cache_xi_fp',
    'cache_miss_fpi',
    'cache_miss_spm',
    'cache_mva_miss',
    'cache_miss_asy',
    'cache_miss_rmf',
    'cache_rmf_lna',
    'cache_miss_sfifo_rmf',
    'cache_miss_fifo_rmf',
    # Importance sampling
    'cache_is',
    'cache_prob_is',
    'cache_miss_is',
    'logmeanexp',
    # TTL-based caches
    'cache_t_hlru',
    'cache_ttl_hlru',
    'cache_ttl_lrua',
    # RRM methods
    'cache_rrm_meanfield_ode',
    'cache_rrm_meanfield',
    'cache_gamma_lp',
    'cache_gamma',
]
