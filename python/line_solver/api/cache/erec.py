"""
Exact Recursive (EREC) Cache Analysis Algorithms.

Native Python implementations of exact recursive methods for cache
analysis, providing precise solutions for small to medium-sized
cache systems.

References:
    Che et al., "Hierarchical Web Caching Systems: Modeling, Design and
    Experimental Results." IEEE JSAC, 2002.
"""

import numpy as np
from scipy import special
from typing import Union, Tuple


def cache_erec(gamma: np.ndarray, m: np.ndarray,
               sigma: np.ndarray = None, k: np.ndarray = None) -> float:
    """
    Compute the cache normalizing constant using exact recursive method.

    With item sizes and per-list storage cost caps this evaluates the
    constrained normalizing constant E(m,k) of Casale-Gast (IEEE/ACM Trans.
    Networking 29(2), 2021), Sec. IX, from the recursion
    E(m,k) = E_i(m,k) + sum_j m_j gamma_ij E_i(m-1_j, k-sigma_i 1_j), with the
    extra boundary E(m,k)=0 whenever some residual cap is negative.

    Args:
        gamma: Cache access factors matrix (n x h), where n is number of
               items and h is number of cache levels.
        m: Cache capacity vector (1 x h or h,).
        sigma: Optional item storage costs (sizes), positive integers (n,).
        k: Optional per-list storage cost caps, non-negative integers (h,).

    Returns:
        Normalizing constant E.
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()
    if sigma is None or k is None or len(np.asarray(sigma).ravel()) == 0 \
            or len(np.asarray(k).ravel()) == 0:
        return cache_erec_aux(gamma, m, gamma.shape[0])
    return _cache_erec_cost(gamma, m,
                            np.asarray(sigma, dtype=np.float64).ravel(),
                            np.asarray(k, dtype=np.float64).ravel())


def _cache_erec_cost(gamma: np.ndarray, m: np.ndarray,
                     sigma: np.ndarray, k: np.ndarray) -> float:
    """Dynamic program over the (residual capacity, residual cost cap) lattice."""
    n, h = gamma.shape
    if len(sigma) != n:
        raise ValueError('The item size vector must have one entry per item.')
    if len(k) != h:
        raise ValueError('The cost cap vector must have one entry per cache list.')
    if np.any(sigma <= 0) or np.any(sigma != np.round(sigma)):
        raise ValueError('Item sizes must be positive integers.')
    if np.any(k < 0) or np.any(k != np.round(k)):
        raise ValueError('Storage cost caps must be non-negative integers.')
    if np.min(m) < 0 or np.sum(m) > n:
        return 0.0
    if np.sum(m) == 0:
        return 1.0
    mv = m.astype(int)
    kv = k.astype(int)
    sv = sigma.astype(int)
    dims = np.concatenate([mv + 1, kv + 1])
    lattice = int(np.prod(dims))
    if lattice > 10000000:
        raise ValueError('The cost-constrained normalizing constant lattice has %d '
                         'states, which exceeds the exact method limit; use the '
                         'sampling method.' % lattice)
    # F(mm,kk) over items 1..t; at t=0 only the empty cache contributes
    F = np.zeros(dims, dtype=np.float64)
    empty = tuple([0] * h + [slice(None)] * h)
    F[empty] = 1.0
    for t in range(n):
        Fprev = F
        F = Fprev.copy()
        for j in range(h):
            g = gamma[t, j]
            if g == 0.0:
                continue
            # shift by one unit of capacity in list j and sigma_t of its cap
            src = [slice(None)] * (2 * h)
            dst = [slice(None)] * (2 * h)
            src[j] = slice(0, dims[j] - 1)
            dst[j] = slice(1, dims[j])
            if sv[t] > dims[h + j] - 1:
                continue
            src[h + j] = slice(0, dims[h + j] - sv[t])
            dst[h + j] = slice(sv[t], dims[h + j])
            mult = np.arange(1, dims[j]).reshape(
                [-1 if d == j else 1 for d in range(2 * h)])
            F[tuple(dst)] += g * mult * Fprev[tuple(src)]
    return float(F[tuple(dims - 1)])


def cache_erec_aux(gamma: np.ndarray, m: np.ndarray, k: int) -> float:
    """
    Auxiliary method for computing cache normalizing constant using
    exact recursive method.

    This method performs the core computation recursively, adjusting
    the size of the input matrix.

    Args:
        gamma: Cache access factors matrix (n x h).
        m: Cache capacity vector (h,).
        k: Current number of rows in the recursive step.

    Returns:
        Normalizing constant for the given configuration.
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()
    h = len(m)  # Number of cache levels

    # Base cases
    m_sum = np.sum(m)
    m_min = np.min(m)

    if m_sum == 0:
        return 1.0

    if m_sum > k or m_min < 0:
        return 0.0

    if k == 1 and m_sum == 1.0:
        # Find the index of the non-zero element in m
        j = np.argmax(m)
        return gamma[0, j]

    # Recursive case
    E = cache_erec_aux(gamma, m, k - 1)

    for j in range(h):
        if m[j] > 0:
            # Create m with one element reduced at position j
            m_oner = m.copy()
            m_oner[j] -= 1
            term = cache_erec_aux(gamma, m_oner, k - 1) * gamma[k - 1, j] * m[j]
            E += term

    return E


def cache_prob_erec(gamma: np.ndarray, m: np.ndarray,
                    sigma: np.ndarray = None, k: np.ndarray = None) -> np.ndarray:
    """
    Compute cache state probabilities using exact recursive method.

    This method calculates the probabilities of the cache being in
    different states based on the cache access factors and capacity. With
    item sizes and per-list storage cost caps it evaluates
    pi_ij = m_j gamma_ij E_i(m-1_j, k-sigma_i 1_j) / E(m,k).

    Args:
        gamma: Cache access factors matrix (n x h), where n is number of
               items and h is number of cache levels.
        m: Cache capacity vector (1 x h or h,).
        sigma: Optional item storage costs (sizes) (n,).
        k: Optional per-list storage cost caps (h,).

    Returns:
        Matrix (n x h+1) containing cache state probabilities.
        Column 0 is miss probability, columns 1..h are hit probabilities
        at each cache level.
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()
    capped = not (sigma is None or k is None
                  or len(np.asarray(sigma).ravel()) == 0
                  or len(np.asarray(k).ravel()) == 0)
    if capped:
        sigma = np.asarray(sigma, dtype=np.float64).ravel()
        k = np.asarray(k, dtype=np.float64).ravel()

    n = gamma.shape[0]  # Number of items
    h = gamma.shape[1]  # Number of cache levels

    E = cache_erec(gamma, m, sigma if capped else None, k if capped else None)
    prob = np.zeros((n, h + 1))

    for i in range(n):
        for j in range(h):
            # Create sub-gamma with row i removed
            sub_gamma = np.delete(gamma, i, axis=0)

            # Create m with one element reduced at position j
            m_oner = m.copy()
            m_oner[j] -= 1

            if capped:
                kij = k.copy()
                kij[j] -= sigma[i]
                if kij[j] < 0:
                    Ei = 0.0
                else:
                    Ei = cache_erec(sub_gamma, m_oner, np.delete(sigma, i), kij)
            else:
                Ei = cache_erec_aux(sub_gamma, m_oner, n - 1)

            if E != 0:
                value = m[j] * gamma[i, j] * Ei / E
            else:
                value = 0.0
            prob[i, j + 1] = value

        # Miss probability is 1 - sum of hit probabilities
        row_sum = np.sum(prob[i, 1:h + 1])
        prob[i, 0] = abs(1 - row_sum)

    return prob


def cache_mva(gamma: np.ndarray, m: np.ndarray
              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                         np.ndarray, np.ndarray, float]:
    """
    Mean Value Analysis for cache systems.

    Computes cache performance metrics using MVA approach.

    Args:
        gamma: Request rate matrix (n x h).
        m: Cache size parameters (h,).

    Returns:
        Tuple of (pi, pi0, pij, x, u, E) containing:
            - pi: Steady-state probabilities
            - pi0: Miss probabilities per item
            - pij: Hit probabilities per item per level (n x h)
            - x: Throughput vector
            - u: Utilization vector
            - E: Normalizing constant
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    n = gamma.shape[0]  # Number of items
    h = gamma.shape[1]  # Number of cache levels

    # Compute normalizing constant
    E = cache_erec(gamma, m)

    # Compute state probabilities
    prob = cache_prob_erec(gamma, m)

    pi0 = prob[:, 0]  # Miss probabilities
    pij = prob[:, 1:]  # Hit probabilities

    # Compute throughputs and utilizations
    total_rates = np.sum(gamma, axis=1)
    x = total_rates * (1 - pi0)  # Throughput (hit rate)
    u = x / (np.sum(x) + 1e-10) if np.sum(x) > 0 else np.zeros(n)

    # pi is the full probability vector
    pi = prob

    return pi, pi0, pij, x, u, E


def cache_cost(gamma: np.ndarray, m: np.ndarray, sigma: np.ndarray,
               k: np.ndarray = None, pij: np.ndarray = None) -> np.ndarray:
    """
    Mean storage cost held by each cache list.

    Evaluates K_j = sum_i sigma_i pi_ij, the expected storage cost of the items
    resident in list j at steady state, as defined in Casale-Gast (IEEE/ACM
    Trans. Networking 29(2), 2021), Sec. IX.

    Args:
        gamma: Cache access factors matrix (n x h).
        m: Cache capacity vector (h,).
        sigma: Item storage costs (sizes) (n,).
        k: Optional per-list storage cost caps (h,).
        pij: Optional precomputed occupancy matrix (n x h+1), column 0 the miss
             probability.

    Returns:
        Mean storage cost of each list (h,).
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    sigma = np.asarray(sigma, dtype=np.float64).ravel()
    n, h = gamma.shape
    if len(sigma) != n:
        raise ValueError('The item size vector must have one entry per item.')
    if pij is None:
        pij = cache_prob_erec(gamma, m, sigma, k)
    pij = np.asarray(pij, dtype=np.float64)
    return np.array([float(sigma @ pij[:, 1 + j]) for j in range(h)])


def cache_cost_pathcheck(gamma: np.ndarray, sigma: np.ndarray, k: np.ndarray,
                         parent: np.ndarray) -> np.ndarray:
    """
    Detect promotion paths blocked by storage cost caps.

    The constrained normalizing constant E(m,k) sums the product form over
    every size-feasible cache state, while under RR-C(m) an item only reaches
    list j by being promoted one list at a time along the path from the miss
    list. A cap on an intermediate list therefore makes size-feasible states
    unreachable and E(m,k) normalizes over states the cache never visits. An
    empty report is a necessary, not sufficient, condition for the two sets to
    agree.

    Args:
        gamma: Cache access factors matrix (n x h).
        sigma: Item storage costs (sizes) (n,).
        k: Per-list storage cost caps (h,).
        parent: Parent list of each list, 0-based with -1 for the lists rooted
                in the miss list.

    Returns:
        Array of blocked triples [item, list, blocking list], 0-based.
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    sigma = np.asarray(sigma, dtype=np.float64).ravel()
    k = np.asarray(k, dtype=np.float64).ravel()
    parent = np.asarray(parent, dtype=int).ravel()
    n, h = gamma.shape
    viol = []
    for i in range(n):
        for j in range(h):
            if gamma[i, j] == 0.0 or sigma[i] > k[j]:
                continue  # item i never resides in list j anyway
            l = parent[j]
            while l >= 0:
                if sigma[i] > k[l]:
                    viol.append([i, j, l])
                    break
                l = parent[l]
    return np.array(viol, dtype=int).reshape(-1, 3)


__all__ = [
    'cache_erec',
    'cache_erec_aux',
    'cache_prob_erec',
    'cache_cost',
    'cache_cost_pathcheck',
    'cache_mva',
]
