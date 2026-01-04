"""
TTL-Based Cache Analysis Methods.

Native Python implementations of Time-To-Live (TTL) based cache analysis
methods, including LRU(m) and hierarchical LRU caches.

Key functions:
    cache_t_lrum: Characteristic times for LRU(m) cache levels
    cache_t_hlru: Characteristic times for hierarchical LRU cache
    cache_ttl_lrum: Steady-state probabilities for TTL-LRU(m) cache
    cache_ttl_hlru: Steady-state probabilities for TTL hierarchical LRU
    cache_ttl_lrua: Steady-state probabilities with arrival-based routing

References:
    Original MATLAB: matlab/src/api/cache/cache_ttl*.m, cache_t*.m
"""

import numpy as np
from scipy.optimize import fsolve
from typing import Tuple, Optional


def cache_t_lrum(gamma: np.ndarray, m: np.ndarray) -> np.ndarray:
    """
    Compute characteristic times for LRU(m) cache levels.

    Uses fixed-point iteration (fsolve) to compute the characteristic time
    for each level of an LRU(m) cache.

    Args:
        gamma: Item popularity probabilities or arrival rates (n x h)
        m: Cache capacity vector (h,)

    Returns:
        Characteristic time for each cache level (h,)

    References:
        Original MATLAB: matlab/src/api/cache/cache_t_lrum.m
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    n, h = gamma.shape

    def time_equations(x):
        """System of equations for characteristic times."""
        F = np.zeros(h)

        # Compute transition probabilities and capacities
        trans = np.zeros((n, h))
        logtrans = np.zeros((n, h))
        denom = np.zeros(n)
        stablecapa = np.zeros(h)

        for k in range(n):
            for j in range(h):
                val = np.exp(gamma[k, j] * x[j]) - 1
                trans[k, j] = max(val, 1e-300)
                logtrans[k, j] = np.log(trans[k, j])
            denom[k] = np.sum(np.exp(np.cumsum(logtrans[k, :])))

        for l in range(h):
            for k in range(n):
                log_prod = np.sum(np.log(trans[k, :l+1]))
                stablecapa[l] += np.exp(log_prod - np.log(1 + denom[k]))
            F[l] = m[l] - stablecapa[l]

        return F

    # Initial guess
    x0 = np.ones(h)

    # Solve
    t, info, ier, mesg = fsolve(time_equations, x0, full_output=True)

    if ier != 1:
        # Try different initial conditions
        for scale in [0.1, 0.5, 2.0, 5.0]:
            t, info, ier, mesg = fsolve(time_equations, x0 * scale, full_output=True)
            if ier == 1:
                break

    return t


def cache_t_hlru(gamma: np.ndarray, m: np.ndarray) -> np.ndarray:
    """
    Compute characteristic times for hierarchical LRU cache levels.

    Uses fixed-point iteration (fsolve) to compute the characteristic time
    for each level of a hierarchical LRU cache.

    Args:
        gamma: Item popularity probabilities or arrival rates (n x h)
        m: Cache capacity vector (h,)

    Returns:
        Characteristic time for each cache level (h,)

    References:
        Original MATLAB: matlab/src/api/cache/cache_t_hlru.m
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    n, h = gamma.shape

    def time_equations(x):
        """System of equations for characteristic times."""
        F = np.zeros(h)

        for a in range(h):
            temp1 = np.ones(n)
            temp2 = np.zeros(n)
            probh = np.zeros(n)

            for k in range(n):
                # Compute temp1: product of (1 - exp(-gamma * t)) for levels 0..a
                for s in range(a + 1):
                    temp1[k] *= (1 - np.exp(-gamma[k, s] * x[s]))

                # Compute temp2
                middtemp2 = 0.0
                for l in range(a):
                    middtemp = 1.0
                    for s in range(l + 1):
                        middtemp *= (1 - np.exp(-gamma[k, s] * x[s]))
                    middtemp2 += middtemp

                temp2[k] = np.exp(-gamma[k, a] * x[a]) * (1 + middtemp2)

                # Hit probability at level a
                denom = temp1[k] + temp2[k]
                probh[k] = temp1[k] / denom if denom > 0 else 0.0

            F[a] = m[a] - np.sum(probh)

        return F

    # Initial guess
    x0 = np.ones(h)

    # Solve
    t, info, ier, mesg = fsolve(time_equations, x0, full_output=True)

    if ier != 1:
        for scale in [0.1, 0.5, 2.0, 5.0]:
            t, info, ier, mesg = fsolve(time_equations, x0 * scale, full_output=True)
            if ier == 1:
                break

    return t


def cache_ttl_lrum(gamma: np.ndarray, m: np.ndarray) -> np.ndarray:
    """
    Compute steady-state probabilities for TTL-based LRU(m) cache.

    Args:
        gamma: Arrival rates per item per list. Can be (n x h) or (1 x n x h+1).
               If 3D, uses gamma[0, :, 1:] as the (n x h) input.
        m: Cache capacity vector (h,)

    Returns:
        Steady-state probability distribution (n x h+1):
            prob[:, 0] = miss probability
            prob[:, 1:] = hit probability at each level

    References:
        Original MATLAB: matlab/src/api/cache/cache_ttl_lrum.m
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    # Handle 3D input (MATLAB convention)
    if gamma.ndim == 3:
        gamma = gamma[0, :, 1:]

    n, h = gamma.shape

    # Compute characteristic times
    t = cache_t_lrum(gamma, m)

    # Compute steady-state probabilities
    probh = np.zeros((n, h))
    prob0 = np.zeros(n)
    trans = np.zeros((n, h))
    denom = np.zeros(n)

    for k in range(n):
        for j in range(h):
            trans[k, j] = np.exp(gamma[k, j] * t[j]) - 1
        denom[k] = np.sum(np.cumprod(trans[k, :]))

    for k in range(n):
        for l in range(h):
            probh[k, l] = np.prod(trans[k, :l+1]) / (1 + denom[k])
        prob0[k] = 1 - np.sum(probh[k, :])

    prob = np.column_stack([prob0, probh])
    return prob


def cache_ttl_hlru(gamma: np.ndarray, m: np.ndarray) -> np.ndarray:
    """
    Compute steady-state probabilities for TTL-based hierarchical LRU cache.

    Args:
        gamma: Arrival rates per item per list. Can be (n x h) or (1 x n x h+1).
               If 3D, uses gamma[0, :, 1:] as the (n x h) input.
        m: Cache capacity vector (h,)

    Returns:
        Steady-state probability distribution (n x 2):
            prob[:, 0] = miss probability
            prob[:, 1] = hit probability at level h

    References:
        Original MATLAB: matlab/src/api/cache/cache_ttl_hlru.m
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    # Handle 3D input
    if gamma.ndim == 3:
        gamma = gamma[0, :, 1:]

    n, h = gamma.shape

    # Compute characteristic times
    t = cache_t_hlru(gamma, m)

    # Compute steady-state probabilities
    probh = np.zeros(n)
    prob0 = np.zeros(n)
    temp1 = np.ones(n)
    temp2 = np.zeros(n)

    for k in range(n):
        # Product of (1 - exp(-gamma * t)) for all levels
        for s in range(h):
            temp1[k] *= (1 - np.exp(-gamma[k, s] * t[s]))

        # Compute temp2
        middtemp2 = 0.0
        for l in range(h - 1):
            middtemp = 1.0
            for s in range(l + 1):
                middtemp *= (1 - np.exp(-gamma[k, s] * t[s]))
            middtemp2 += middtemp

        temp2[k] = np.exp(-gamma[k, h-1] * t[h-1]) * (1 + middtemp2)

        denom = temp1[k] + temp2[k]
        probh[k] = temp1[k] / denom if denom > 0 else 0.0
        prob0[k] = 1 / (denom * np.exp(gamma[k, h-1] * t[h-1])) if denom > 0 else 1.0

    prob = np.column_stack([prob0, probh])
    return prob


def cache_ttl_lrua(lambd: np.ndarray, R: list, m: np.ndarray,
                   seed: int = 23000) -> np.ndarray:
    """
    Compute steady-state probabilities for TTL-LRU cache with arrival routing.

    Uses fixed-point iteration with DTMC solving for cache systems with
    multiple users, items, and levels with routing.

    Args:
        lambd: Arrival rates per user per item per list (u x n x h+1)
        R: Routing probability structure (list of n matrices, each (h+1 x h+1))
        m: Cache capacity vector (h,)
        seed: Random seed for initialization (default: 23000)

    Returns:
        Steady-state probability distribution (n x h+1)

    References:
        Original MATLAB: matlab/src/api/cache/cache_ttl_lrua.m
    """
    np.random.seed(seed)

    lambd = np.asarray(lambd, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    u = lambd.shape[0]  # number of users
    n = lambd.shape[1]  # number of items
    h = lambd.shape[2] - 1  # number of lists

    def ttl_tree_time(x):
        """Compute capacity difference for given characteristic times."""
        steadystateprob = np.zeros((n, h + 1))
        randprob = np.zeros((n, h + 1))
        avgtime = np.zeros((n, h + 1))
        capa = np.zeros(h)
        rpdenominator = np.zeros(n)

        for i in range(n):
            # Build transition matrix
            Ri = R[i]
            transmatrix = np.zeros((h + 1, h + 1))

            for j in range(h + 1):
                leafnodes = np.where(Ri[j, :] > 0)[0]
                for k in leafnodes:
                    if j == 0:
                        transmatrix[j, k] = Ri[j, k]
                    else:
                        transmatrix[j, k] = (1 - np.exp(-lambd[0, i, j] * x[j-1])) * Ri[j, k]
                    if j != k and k > 0:
                        transmatrix[k, j] = np.exp(-lambd[0, i, k] * x[k-1])

            # Remove disconnected nodes
            connected = []
            for j in range(h + 1):
                if np.any(transmatrix[j, :] > 0) or np.any(transmatrix[:, j] > 0):
                    connected.append(j)

            if len(connected) == 0:
                continue

            # Extract submatrix for connected nodes
            submatrix = transmatrix[np.ix_(connected, connected)]

            # Solve DTMC
            dtmcprob = _dtmc_solve(submatrix)

            for idx, node in enumerate(connected):
                steadystateprob[i, node] = dtmcprob[idx]
                if node > 0:
                    rate = lambd[0, i, node]
                    if rate > 0:
                        avgtime[i, node] = (1 - np.exp(-rate * x[node-1])) / rate
                    else:
                        avgtime[i, node] = 0
                else:
                    rate = lambd[0, i, node]
                    avgtime[i, node] = 1 / rate if rate > 0 else 1.0
                rpdenominator[i] += steadystateprob[i, node] * avgtime[i, node]

            for idx, node in enumerate(connected):
                if rpdenominator[i] > 0:
                    randprob[i, node] = (steadystateprob[i, node] *
                                         avgtime[i, node] / rpdenominator[i])

        # Compute capacity difference
        F = np.zeros(h)
        for l in range(h):
            capa[l] = np.sum(randprob[:, l + 1])
            F[l] = m[l] - capa[l]

        return F, randprob

    # Initial guess with random seed
    x0 = np.random.uniform(0, 10, h)

    # Solve using fsolve
    def objective(x):
        F, _ = ttl_tree_time(x)
        return F

    t, info, ier, mesg = fsolve(objective, x0, full_output=True)

    # Get final probabilities
    _, prob = ttl_tree_time(t)

    return prob


def _dtmc_solve(P: np.ndarray) -> np.ndarray:
    """
    Solve for stationary distribution of a DTMC.

    Args:
        P: Transition probability matrix

    Returns:
        Stationary distribution vector
    """
    n = P.shape[0]

    if n == 0:
        return np.array([])

    if n == 1:
        return np.array([1.0])

    # Normalize rows to make it a proper transition matrix
    row_sums = np.sum(P, axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    P = P / row_sums

    # Solve (P^T - I) * pi = 0 with sum(pi) = 1
    # Use the method: replace last equation with sum = 1
    A = P.T - np.eye(n)
    A[-1, :] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0

    try:
        pi = np.linalg.solve(A, b)
        pi = np.maximum(pi, 0)  # Ensure non-negative
        pi = pi / np.sum(pi)  # Normalize
    except np.linalg.LinAlgError:
        # Fallback: power method
        pi = np.ones(n) / n
        for _ in range(1000):
            pi_new = P.T @ pi
            pi_new = pi_new / np.sum(pi_new)
            if np.max(np.abs(pi_new - pi)) < 1e-10:
                break
            pi = pi_new
        pi = pi_new

    return pi


__all__ = [
    'cache_t_lrum',
    'cache_t_hlru',
    'cache_ttl_lrum',
    'cache_ttl_hlru',
    'cache_ttl_lrua',
]
