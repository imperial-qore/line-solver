"""
Cache Analysis Algorithms.

Native Python implementations for analyzing cache systems, including
exact recursive methods, singular perturbation methods, importance
sampling, TTL-based caches, and various approximation techniques.

Key algorithms:
    cache_erec: Exact recursive normalizing constant
    cache_prob_erec: Exact recursive state probabilities
    cache_spm: Singular perturbation method
    cache_miss: Miss rate computation
    cache_is: Importance sampling
    cache_ttl_*: TTL-based cache analysis
    cache_rrm_*: Random replacement model
"""

from .erec import (
    cache_erec,
    cache_erec_aux,
    cache_prob_erec,
    cache_mva,
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
)

from .sampling import (
    cache_is,
    cache_prob_is,
    cache_miss_is,
    logmeanexp,
)

from .ttl import (
    cache_t_lrum,
    cache_t_hlru,
    cache_ttl_lrum,
    cache_ttl_hlru,
    cache_ttl_lrua,
)

from .rrm import (
    cache_rrm_meanfield_ode,
    cache_rrm_meanfield,
    cache_gamma_lp,
)

__all__ = [
    # Exact recursive methods
    'cache_erec',
    'cache_erec_aux',
    'cache_prob_erec',
    'cache_mva',
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
    # Importance sampling
    'cache_is',
    'cache_prob_is',
    'cache_miss_is',
    'logmeanexp',
    # TTL-based caches
    'cache_t_lrum',
    'cache_t_hlru',
    'cache_ttl_lrum',
    'cache_ttl_hlru',
    'cache_ttl_lrua',
    # RRM methods
    'cache_rrm_meanfield_ode',
    'cache_rrm_meanfield',
    'cache_gamma_lp',
]
