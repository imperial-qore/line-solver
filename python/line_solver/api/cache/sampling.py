"""
Importance Sampling Methods for Cache Analysis.

Native Python implementations of Monte Carlo importance sampling methods
for cache performance analysis.

Key functions:
    cache_is: Importance sampling estimation of normalizing constant
    cache_prob_is: Importance sampling estimation of hit probabilities
    cache_miss_is: Importance sampling estimation of miss rates

References:
    Original MATLAB: matlab/src/api/cache/cache_is.m, cache_prob_is.m
"""

import numpy as np
from typing import Tuple, Optional
from scipy.special import gammaln
from .erec import cache_erec, cache_prob_erec


def factln(n: float) -> float:
    """Compute log(n!) using log-gamma function."""
    if n <= 0:
        return 0.0
    return gammaln(n + 1)


def logmeanexp(x: np.ndarray) -> float:
    """
    Compute log(mean(exp(x))) in a numerically stable way.

    Uses the log-sum-exp trick for numerical stability.
    """
    x = np.asarray(x).flatten()
    x_max = np.max(x)
    if not np.isfinite(x_max):
        return float('-inf')
    return x_max + np.log(np.mean(np.exp(x - x_max)))


def _assign_items_to_levels(m: np.ndarray, selected: np.ndarray) -> list:
    """
    Assign selected items to cache levels uniformly at random.

    Args:
        m: Cache capacity vector (h,)
        selected: Selected item indices

    Returns:
        List of arrays, one per level, containing assigned item indices
    """
    h = len(m)
    m = m.astype(int)

    # Shuffle selected items
    perm = np.random.permutation(len(selected))
    shuffled = selected[perm]

    # Assign to levels
    assignment = []
    idx = 0
    for j in range(h):
        assignment.append(shuffled[idx:idx + m[j]])
        idx += m[j]

    return assignment


def cache_is(gamma: np.ndarray, m: np.ndarray,
             samples: int = 100000) -> Tuple[float, float]:
    """
    Importance sampling estimation of cache normalizing constant.

    Estimates the normalizing constant for cache models using Monte Carlo
    importance sampling.

    Args:
        gamma: Item popularity probabilities (n x h matrix)
        m: Cache capacity vector (h,)
        samples: Number of Monte Carlo samples (default: 100000)

    Returns:
        Tuple of (E, lE) where:
            - E: Normalizing constant estimate
            - lE: Log of normalizing constant

    References:
        Original MATLAB: matlab/src/api/cache/cache_is.m
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    # Remove items with zero gamma
    row_sums = np.sum(gamma, axis=1)
    gamma = gamma[row_sums > 0, :]

    n, h = gamma.shape
    mt = int(np.sum(m))  # total cache capacity

    # Edge cases
    if n == 0 or mt == 0:
        return 1.0, 0.0

    if n < mt:
        return 0.0, float('-inf')

    if n == mt:
        # All items must be in cache - use exact method
        E = cache_erec(gamma, m)
        lE = np.log(E) if E > 0 else float('-inf')
        return E, lE

    # Pre-compute log(gamma)
    log_gamma = np.log(gamma + 1e-300)

    # Pre-compute log(m(j)!)
    log_m_fact = sum(factln(mj) for mj in m)

    # Log of number of ways to choose mt items from n
    log_combinations = factln(n) - factln(mt) - factln(n - mt)

    lZ_samples = np.zeros(samples)

    for s in range(samples):
        # Sample mt items uniformly without replacement
        selected = np.random.choice(n, mt, replace=False)

        # Assign items to levels uniformly at random
        assignment = _assign_items_to_levels(m, selected)

        # Compute log of unnormalized state probability
        log_state_prob = log_m_fact
        for j in range(h):
            items_in_level = assignment[j]
            for i in items_in_level:
                log_state_prob += log_gamma[i, j]

        # Proposal probability is 1/(C(n,mt) * multinomial(mt; m))
        log_multinomial = factln(mt) - log_m_fact
        log_proposal = -log_combinations - log_multinomial

        lZ_samples[s] = log_state_prob - log_proposal

    lE = logmeanexp(lZ_samples)
    E = np.exp(lE) if np.isfinite(lE) else 0.0

    return E, lE


def cache_prob_is(gamma: np.ndarray, m: np.ndarray,
                  samples: int = 100000) -> np.ndarray:
    """
    Importance sampling estimation of cache hit probabilities.

    Estimates cache hit probability distribution using Monte Carlo
    importance sampling.

    Args:
        gamma: Item popularity probabilities (n x h matrix)
        m: Cache capacity vector (h,)
        samples: Number of Monte Carlo samples (default: 100000)

    Returns:
        Cache hit probability matrix (n x h+1):
            prob[i, 0] = miss probability for item i
            prob[i, 1+j] = hit probability for item i at level j

    References:
        Original MATLAB: matlab/src/api/cache/cache_prob_is.m
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    n, h = gamma.shape
    mt = int(np.sum(m))

    # Initialize probability matrix
    prob = np.zeros((n, h + 1))

    # Edge cases
    if n == 0 or mt == 0:
        prob[:, 0] = 1.0  # all items miss
        return prob

    if n < mt:
        prob[:, 0] = 1.0
        return prob

    if n == mt:
        # All items in cache - use exact method
        return cache_prob_erec(gamma, m)

    log_gamma = np.log(gamma + 1e-300)
    log_m_fact = sum(factln(mj) for mj in m)
    log_combinations = factln(n) - factln(mt) - factln(n - mt)
    log_multinomial = factln(mt) - log_m_fact

    # Accumulators
    item_level_weight = np.zeros((n, h))
    total_weight = 0.0

    # Scale factor to avoid overflow
    scale = 50.0

    for s in range(samples):
        # Sample mt items uniformly without replacement
        selected = np.random.choice(n, mt, replace=False)

        # Assign items to levels
        assignment = _assign_items_to_levels(m, selected)

        # Compute unnormalized state probability
        log_state_prob = log_m_fact
        for j in range(h):
            for i in assignment[j]:
                log_state_prob += log_gamma[i, j]

        # Compute proposal probability
        log_proposal = -log_combinations - log_multinomial

        # Importance weight
        log_is_weight = log_state_prob - log_proposal
        is_weight = np.exp(log_is_weight - scale)

        # Update accumulators
        total_weight += is_weight
        for j in range(h):
            for i in assignment[j]:
                item_level_weight[i, j] += is_weight

    # Compute probabilities
    if total_weight > 0:
        for i in range(n):
            for j in range(h):
                prob[i, 1 + j] = item_level_weight[i, j] / total_weight
            prob[i, 0] = max(0, 1 - np.sum(prob[i, 1:]))
    else:
        prob[:, 0] = 1.0

    return prob


def cache_miss_is(gamma: np.ndarray, m: np.ndarray,
                  lambd: Optional[np.ndarray] = None,
                  samples: int = 100000
                  ) -> Tuple[float, Optional[np.ndarray], Optional[np.ndarray],
                             Optional[np.ndarray], float]:
    """
    Importance sampling estimation of cache miss rates.

    Computes global, per-user, and per-item miss rates using Monte Carlo
    importance sampling.

    Args:
        gamma: Item popularity probabilities (n x h matrix)
        m: Cache capacity vector (h,)
        lambd: Optional arrival rates per user per item (u x n x h+1)
        samples: Number of Monte Carlo samples (default: 100000)

    Returns:
        Tuple of (M, MU, MI, pi0, lE) where:
            - M: Global miss rate
            - MU: Per-user miss rate or None
            - MI: Per-item miss rate or None
            - pi0: Per-item miss probability or None
            - lE: Log of normalizing constant

    References:
        Original MATLAB: matlab/src/api/cache/cache_miss_is.m
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    n = gamma.shape[0]

    # Compute normalizing constant via importance sampling
    _, lE = cache_is(gamma, m, samples)

    # Compute hit probabilities via importance sampling
    pij = cache_prob_is(gamma, m, samples)

    # Extract miss probabilities (first column)
    pi0 = pij[:, 0]

    if lambd is None:
        M = np.mean(pi0)
        return M, None, None, pi0, lE

    lambd = np.asarray(lambd, dtype=np.float64)
    u = lambd.shape[0]

    # Per-user miss rate
    MU = np.zeros(u)
    for v in range(u):
        for k in range(n):
            MU[v] += lambd[v, k, 0] * pi0[k]

    # Per-item miss rate
    MI = np.zeros(n)
    for k in range(n):
        MI[k] = np.sum(lambd[:, k, 0]) * pi0[k]

    # Global miss rate
    M = np.sum(MI)

    return M, MU, MI, pi0, lE


__all__ = [
    'cache_is',
    'cache_prob_is',
    'cache_miss_is',
    'logmeanexp',
]
