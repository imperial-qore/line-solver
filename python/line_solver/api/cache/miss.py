"""
Cache Miss Rate Computation Methods.

Native Python implementations of various methods for computing cache miss rates,
including exact methods, fixed-point iteration, and approximations.

Key functions:
    cache_miss: Basic miss rate computation using exact recursive method
    cache_miss_fpi: Miss rates via fixed-point iteration
    cache_miss_spm: Miss rates via saddle-point method
    cache_miss_is: Miss rates via importance sampling
    cache_mva_miss: Miss rates via Mean Value Analysis

References:
    Original MATLAB: matlab/src/api/cache/cache_miss*.m
"""

import numpy as np
from typing import Tuple, Optional
from .erec import cache_erec, cache_erec_aux
from .spm import cache_spm, cache_prob_spm


def cache_miss(gamma: np.ndarray, m: np.ndarray,
               lambd: Optional[np.ndarray] = None
               ) -> Tuple[float, Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Compute miss rates for a cache system using exact recursive method.

    Args:
        gamma: Item popularity probabilities (n x h matrix)
        m: Cache capacity vector (h,)
        lambd: Optional arrival rates per user per item (u x n x h+1)

    Returns:
        Tuple of (M, MU, MI, pi0) where:
            - M: Global miss rate
            - MU: Per-user miss rate (u,) or None
            - MI: Per-item miss rate (n,) or None
            - pi0: Per-item miss probability (n,) or None

    References:
        Original MATLAB: matlab/src/api/cache/cache_miss.m
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    n = gamma.shape[0]

    # Compute global miss rate
    ma = m.copy()
    ma[0] = ma[0] + 1
    E = cache_erec(gamma, m)
    Ea = cache_erec(gamma, ma)
    M = Ea / E if E != 0 else 1.0

    if lambd is None:
        return M, None, None, None

    lambd = np.asarray(lambd, dtype=np.float64)
    u = lambd.shape[0]  # number of users

    # Compute per-item miss probabilities
    pi0 = np.zeros(n)
    for k in range(n):
        # Remove item k from gamma
        gamma_sub = np.delete(gamma, k, axis=0)
        E_sub = cache_erec_aux(gamma_sub, m, n - 1)
        pi0[k] = E_sub / E if E != 0 else 1.0

    # Compute per-user miss rate
    MU = np.zeros(u)
    for v in range(u):
        for k in range(n):
            MU[v] += lambd[v, k, 0] * pi0[k]

    # Compute per-item miss rate
    MI = np.zeros(n)
    for k in range(n):
        MI[k] = np.sum(lambd[:, k, 0]) * pi0[k]

    return M, MU, MI, pi0


def cache_xi_fp(gamma: np.ndarray, m: np.ndarray,
                xi_init: Optional[np.ndarray] = None,
                tol: float = 1e-14, max_iter: int = 10000
                ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Fixed-point iteration for computing cache performance metrics.

    Computes cache performance metrics including Lagrange multipliers,
    miss probabilities, and hit probabilities.

    Args:
        gamma: Item popularity probabilities (n x h)
        m: Cache capacity vector (h,)
        xi_init: Optional initial guess for Lagrange multipliers
        tol: Convergence tolerance (default: 1e-14)
        max_iter: Maximum iterations (default: 10000)

    Returns:
        Tuple of (xi, pi0, pij, it) where:
            - xi: Converged Lagrange multipliers (h,)
            - pi0: Miss probability per item (n,)
            - pij: Hit probability per item per list (n x h)
            - it: Number of iterations

    References:
        Original MATLAB: matlab/src/api/cache/cache_xi_fp.m
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    n, h = gamma.shape

    # Initialize
    pi0 = np.ones(n) / (h + 1)
    pij = np.zeros((n, h))

    if xi_init is not None:
        xi = np.asarray(xi_init, dtype=np.float64).ravel()
    else:
        xi = np.zeros(h)
        for l in range(h):
            mean_gamma = np.mean(gamma[:, l])
            if mean_gamma > 0:
                xi[l] = m[l] / mean_gamma / (n + np.sum(m) - 1)
            else:
                xi[l] = 1.0

    for it in range(max_iter):
        pi0_old = pi0.copy()

        # Update xi: m / (pi0 * gamma)
        for l in range(h):
            denom = np.dot(pi0_old, gamma[:, l])
            if denom > 0:
                xi[l] = m[l] / denom
            else:
                xi[l] = tol

        # Update pij: |gamma * xi| / |1 + gamma @ xi|
        gamma_xi = gamma @ xi
        for i in range(n):
            denom = np.abs(1 + gamma_xi[i])
            if denom > tol:
                for j in range(h):
                    pij[i, j] = np.abs(gamma[i, j] * xi[j]) / denom
            else:
                pij[i, :] = 0.0

        # Update pi0
        pi0 = np.maximum(tol, 1 - np.sum(pij, axis=1))

        # Check convergence
        delta = np.linalg.norm(np.abs(1 - pi0 / pi0_old), ord=1)
        if delta < tol:
            xi[xi < 0] = tol
            return xi, pi0, pij, it + 1

    xi[xi < 0] = tol
    return xi, pi0, pij, max_iter


def cache_miss_fpi(gamma: np.ndarray, m: np.ndarray,
                   lambd: Optional[np.ndarray] = None
                   ) -> Tuple[float, Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Compute cache miss rates using fixed-point iteration method.

    Args:
        gamma: Item popularity probabilities (n x h)
        m: Cache capacity vector (h,)
        lambd: Optional arrival rates per user per item (u x n x h+1)

    Returns:
        Tuple of (M, MU, MI, pi0) where:
            - M: Global miss rate
            - MU: Per-user miss rate or None
            - MI: Per-item miss rate or None
            - pi0: Per-item miss probability or None

    References:
        Original MATLAB: matlab/src/api/cache/cache_miss_fpi.m
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    n = gamma.shape[0]

    # Compute xi using fixed-point iteration
    xi, _, _, _ = cache_xi_fp(gamma, m)

    # Compute per-item miss rates
    MI = np.zeros(n)

    if lambd is not None:
        lambd = np.asarray(lambd, dtype=np.float64)
        for i in range(n):
            gamma_xi = np.dot(gamma[i, :], xi)
            MI[i] = np.sum(lambd[:, i, 0]) / (1 + gamma_xi)
    else:
        for i in range(n):
            gamma_xi = np.dot(gamma[i, :], xi)
            MI[i] = 1.0 / (1 + gamma_xi)

    M = np.sum(MI)

    if lambd is None:
        return M, None, MI, None

    # Compute per-item miss probability
    pi0 = np.zeros(n)
    for i in range(n):
        gamma_xi = np.dot(gamma[i, :], xi)
        pi0[i] = 1.0 / (1 + gamma_xi)

    # Compute per-user miss rate
    u = lambd.shape[0]
    MU = np.zeros(u)
    for v in range(u):
        for i in range(n):
            MU[v] += lambd[v, i, 0] * pi0[i]

    return M, MU, MI, pi0


def cache_miss_spm(gamma: np.ndarray, m: np.ndarray,
                   lambd: Optional[np.ndarray] = None
                   ) -> Tuple[float, Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], float]:
    """
    Compute cache miss rates using saddle-point method.

    Args:
        gamma: Item popularity probabilities (n x h)
        m: Cache capacity vector (h,)
        lambd: Optional arrival rates per user per item (u x n x h+1)

    Returns:
        Tuple of (M, MU, MI, pi0, lE) where:
            - M: Global miss rate
            - MU: Per-user miss rate or None
            - MI: Per-item miss rate or None
            - pi0: Per-item miss probability or None
            - lE: Log of normalizing constant

    References:
        Original MATLAB: matlab/src/api/cache/cache_miss_spm.m
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    n = gamma.shape[0]

    # Compute global miss rate
    ma = m.copy()
    ma[0] = ma[0] + 1

    _, lE, xi = cache_spm(gamma, m)
    _, lEa, _ = cache_spm(gamma, ma)

    M = np.exp(lEa - lE)

    if lambd is None:
        return M, None, None, None, lE

    lambd = np.asarray(lambd, dtype=np.float64)
    u = lambd.shape[0]

    # Compute per-item miss probabilities
    pi0 = np.zeros(n)
    lE1 = np.zeros(n)

    for k in range(n):
        if np.sum(gamma[k, :]) > 0:
            # Remove item k
            gamma_sub = np.delete(gamma, k, axis=0)
            _, lE1[k], _ = cache_spm(gamma_sub, m)
            pi0[k] = np.exp(lE1[k] - lE)

            # Check validity and recompute if needed
            if pi0[k] > 1 or pi0[k] < 0:
                _, lE1[k], _ = cache_spm(gamma_sub, m)
                pi0[k] = np.exp(lE1[k] - lE)

    # Compute per-user miss rate
    MU = np.zeros(u)
    for v in range(u):
        for k in range(n):
            MU[v] += lambd[v, k, 0] * pi0[k]

    # Compute per-item miss rate
    MI = np.zeros(n)
    for k in range(n):
        if np.sum(gamma[k, :]) > 0:
            MI[k] = np.sum(lambd[:, k, 0]) * np.exp(lE1[k] - lE)

    return M, MU, MI, pi0, lE


def cache_mva_miss(p: np.ndarray, m: np.ndarray,
                   R: np.ndarray) -> Tuple[float, np.ndarray]:
    """
    Compute cache miss rates using Mean Value Analysis.

    Args:
        p: Item popularity probabilities (n,)
        m: Cache capacity vector (h,)
        R: Routing probabilities (h+1 x n)

    Returns:
        Tuple of (M, Mk) where:
            - M: Global miss rate
            - Mk: Per-item miss rate (n,)

    References:
        Original MATLAB: matlab/src/api/cache/cache_mva_miss.m
    """
    p = np.asarray(p, dtype=np.float64).ravel()
    m = np.asarray(m, dtype=np.float64).ravel()
    R = np.asarray(R, dtype=np.float64)

    n = len(p)
    h = len(m)

    # Base case
    if np.sum(m) == 0 or np.min(m) < 0:
        Mk = np.ones(n)
        M = np.dot(p, Mk)
        return M, Mk

    # Recursive computation
    w = np.zeros((n, h))

    for j in range(h):
        # Compute miss rates with reduced capacity at level j
        m_reduced = m.copy()
        m_reduced[j] = max(0, m_reduced[j] - 1)
        _, Mj = cache_mva_miss(p, m_reduced, R)

        for k in range(n):
            # Product of routing probs up to level j
            R_prod = np.prod(R[:j+1, k])
            w[k, j] = R_prod * (p[k] ** (j + 1)) * np.abs(Mj[k])

    # Compute x
    x = np.zeros(h)
    for j in range(h):
        w_sum = np.sum(np.abs(w[:, j]))
        x[j] = 1.0 / w_sum if w_sum > 0 else 0.0

    # Compute per-item miss rates
    Mk = np.ones(n)
    for k in range(n):
        for j in range(h):
            Mk[k] -= x[j] * m[j] * w[k, j]

    Mk = np.abs(Mk)
    M = np.dot(p, Mk)

    return M, Mk


__all__ = [
    'cache_miss',
    'cache_xi_fp',
    'cache_miss_fpi',
    'cache_miss_spm',
    'cache_mva_miss',
]
