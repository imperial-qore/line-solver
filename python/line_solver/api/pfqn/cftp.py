"""
Coupling From The Past (CFTP) exact/approximate stationary sampler for closed networks.

Reference: S. Kijima and T. Matsui, "Approximate/Perfect Samplers for Closed Jackson
Networks", Proc. Winter Simulation Conference, 2005.
"""

import numpy as np
from typing import Tuple


def pfqn_cftp(L: np.ndarray, N: int, S: np.ndarray = None,
              nsamples: int = 1, method: str = 'cftp') -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Exact (perfect) stationary state sampling for closed single-class multiserver
    product-form networks via monotone Coupling From The Past.

    Args:
        L: Service demands (stations,), L[i] = theta_i / mu_i
        N: Total closed population (scalar)
        S: Number of servers per station (default 1; use np.inf for infinite server).
           If None, defaults to ones.
        nsamples: Number of independent samples to draw (default 1)
        method: 'cftp' for exact/perfect sampling (default) or 'approx' for
                rapidly-mixing approximate sampler M_A

    Returns:
        Q: Empirical mean queue length per station (1 x stations)
        X: Sampled states, one per row (nsamples x stations), each row sums to N
        T: Per-sample coalescence horizon ('cftp') or mixing steps used
           ('approx'), as (nsamples,)

    Raises:
        ValueError: If less than 2 stations or non-positive demands
    """
    L = np.asarray(L, dtype=float).flatten()
    M = len(L)

    if S is None:
        S = np.ones(M)
    else:
        S = np.asarray(S, dtype=float).flatten()

    K = int(round(N))

    if M < 2:
        raise ValueError("pfqn_cftp: at least two stations required")
    if np.any(L <= 0):
        raise ValueError("pfqn_cftp: all demands L must be strictly positive")

    logfac = _precompute_logfac(L, S, K, M)
    logL = np.log(L)

    X = np.zeros((nsamples, M), dtype=int)
    T = np.zeros(nsamples, dtype=int)

    for smp in range(nsamples):
        if method.lower() == 'cftp':
            X[smp, :], T[smp] = _draw_cftp(logL, logfac, K, M)
        elif method.lower() == 'approx':
            X[smp, :], T[smp] = _draw_approx(logL, logfac, K, M)
        else:
            raise ValueError(f"pfqn_cftp: unknown method '{method}'")

    Q = np.mean(X, axis=0)
    return Q, X, T


def _precompute_logfac(L: np.ndarray, S: np.ndarray, K: int, M: int) -> np.ndarray:
    """
    Precompute cumulative log service factors.
    logfac[i, m+1] = sum_{t=1}^m log(min(t, S[i])), m = 0..K
    So log alpha_i(m) = m*log(L[i]) - logfac[i, m+1]
    """
    logfac = np.zeros((M, K + 1))
    for i in range(M):
        acc = 0.0
        for m in range(1, K + 1):
            acc += np.log(min(m, S[i]))
            logfac[i, m] = acc
    return logfac


def _draw_cftp(logL: np.ndarray, logfac: np.ndarray, K: int, M: int) -> Tuple[np.ndarray, int]:
    """One exact draw via monotone CFTP."""
    u = np.array([], dtype=float)
    Tback = 1

    while True:
        old = len(u)
        # Prepend randomness for newly exposed older steps T .. 2T-1
        new_u = np.random.rand(Tback - old)
        u = np.concatenate([new_u, u])

        xU = np.zeros(M, dtype=int)
        xU[0] = K
        xL = np.zeros(M, dtype=int)
        xL[-1] = K

        for t in range(Tback - 1, -1, -1):
            xU = _monotone_update(xU, u[Tback - 1 - t], logL, logfac, M)
            xL = _monotone_update(xL, u[Tback - 1 - t], logL, logfac, M)

        if np.array_equal(xU, xL):
            return xU, Tback

        Tback *= 2


def _draw_approx(logL: np.ndarray, logfac: np.ndarray, K: int, M: int,
                 eps: float = 1e-2) -> Tuple[np.ndarray, int]:
    """Approximate rapidly-mixing sampler M_A."""
    steps = int(np.ceil(M * (M - 1) / 2 * np.log(K / eps)))
    x = np.zeros(M, dtype=int)
    x[0] = K

    for _ in range(steps):
        p = np.random.choice(M, 2, replace=False)
        i, j = p[0], p[1]
        k = x[i] + x[j]
        l = _split_index(logL, logfac, i, j, k, np.random.rand())
        x[i] = l
        x[j] = k - l

    return x, steps


def _monotone_update(x: np.ndarray, u: float, logL: np.ndarray,
                     logfac: np.ndarray, M: int) -> np.ndarray:
    """Single uniform u in [0,1) encodes pair index and split Lambda."""
    lam = 1 + u * (M - 1)          # in [1, M)
    jj = int(np.floor(lam))        # 1-based pair index (jj, jj+1), 1 <= jj <= M-1
    if jj > M - 1:
        jj = M - 1
    Lambda = lam - jj              # fractional part in [0, 1)
    p = jj - 1                     # 0-based left station of the adjacent pair
    k = x[p] + x[p + 1]
    l = _split_index(logL, logfac, p, p + 1, k, Lambda)
    x[p] = l
    x[p + 1] = k - l
    return x


def _split_index(logL: np.ndarray, logfac: np.ndarray, i: int, j: int,
                 k: int, Lambda: float) -> int:
    """Inverse-CDF split: smallest l with Lambda <= g^k_{ij}(l)."""
    s = np.arange(k + 1)
    lw = (s * logL[i] - logfac[i, s]) + ((k - s) * logL[j] - logfac[j, k - s])
    lw = lw - np.max(lw)
    w = np.exp(lw)
    cdf = np.cumsum(w)
    cdf = cdf / cdf[-1]

    l = np.searchsorted(cdf, Lambda)
    if l > k:
        l = k

    return l
