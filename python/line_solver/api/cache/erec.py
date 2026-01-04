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


def cache_erec(gamma: np.ndarray, m: np.ndarray) -> float:
    """
    Compute the cache normalizing constant using exact recursive method.

    This method serves as a wrapper that calls the auxiliary function
    to perform the actual computation.

    Args:
        gamma: Cache access factors matrix (n x h), where n is number of
               items and h is number of cache levels.
        m: Cache capacity vector (1 x h or h,).

    Returns:
        Normalizing constant E.
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()
    k = gamma.shape[0]  # Number of items
    return cache_erec_aux(gamma, m, k)


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


def cache_prob_erec(gamma: np.ndarray, m: np.ndarray) -> np.ndarray:
    """
    Compute cache state probabilities using exact recursive method.

    This method calculates the probabilities of the cache being in
    different states based on the cache access factors and capacity.

    Args:
        gamma: Cache access factors matrix (n x h), where n is number of
               items and h is number of cache levels.
        m: Cache capacity vector (1 x h or h,).

    Returns:
        Matrix (n x h+1) containing cache state probabilities.
        Column 0 is miss probability, columns 1..h are hit probabilities
        at each cache level.
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    n = gamma.shape[0]  # Number of items
    h = gamma.shape[1]  # Number of cache levels

    E = cache_erec(gamma, m)
    prob = np.zeros((n, h + 1))

    for i in range(n):
        for j in range(h):
            # Create sub-gamma with row i removed
            sub_gamma = np.delete(gamma, i, axis=0)

            # Create m with one element reduced at position j
            m_oner = m.copy()
            m_oner[j] -= 1

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


__all__ = [
    'cache_erec',
    'cache_erec_aux',
    'cache_prob_erec',
    'cache_mva',
]
