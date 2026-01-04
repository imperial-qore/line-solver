"""
Logistic Sampling approximation for normalizing constant.

Reference:
G. Casale. Accelerating performance inference over closed systems by
asymptotic methods. ACM SIGMETRICS 2017.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np
from scipy.stats import multivariate_normal
from scipy.special import gammaln
from typing import Tuple


def _factln(n: np.ndarray) -> np.ndarray:
    """Compute log(n!) using log-gamma function."""
    return gammaln(n + 1)


def _multinomialln(n: np.ndarray) -> float:
    """Compute log of multinomial coefficient."""
    return float(_factln(np.sum(n)) - np.sum(_factln(n)))


def _allbut(y: np.ndarray, xset: int) -> np.ndarray:
    """Return all elements except at index xset."""
    mask = np.ones(len(y), dtype=bool)
    mask[xset] = False
    return y[mask]


def _pfqn_le_fpi(L: np.ndarray, N: np.ndarray) -> Tuple[np.ndarray, list]:
    """Fixed point iteration to find mode of Gaussian."""
    M, R = L.shape
    u = np.ones(M) / M
    u_1 = np.inf * np.ones(M)
    d = []

    while np.linalg.norm(u - u_1, 1) > 1e-10:
        u_1 = u.copy()
        for i in range(M):
            u[i] = 1 / (np.sum(N) + M)
            for r in range(R):
                denom = np.dot(u_1, L[:, r])
                if abs(denom) > 1e-14:
                    u[i] += (N[r] / (np.sum(N) + M)) * L[i, r] * u_1[i] / denom
        d.append(np.abs(u - u_1))

    return u, d


def _pfqn_le_fpiZ(L: np.ndarray, N: np.ndarray, Z: np.ndarray) -> Tuple[np.ndarray, float, list]:
    """Fixed point iteration with think times."""
    M, R = L.shape
    eta = np.sum(N) + M
    u = np.ones(M) / M
    v = eta + 1
    u_1 = np.inf * np.ones(M)
    v_1 = np.inf
    d = []

    while np.linalg.norm(u - u_1, 1) > 1e-10:
        u_1 = u.copy()
        v_1 = v

        for ist in range(M):
            u[ist] = 1 / eta
            for r in range(R):
                denom = Z[r] + v * np.dot(u_1, L[:, r])
                if abs(denom) > 1e-14:
                    u[ist] += (N[r] / eta) * (Z[r] + v * L[ist, r]) * u_1[ist] / denom

        xi = np.zeros(R)
        for r in range(R):
            denom = Z[r] + v * np.dot(u_1, L[:, r])
            if abs(denom) > 1e-14:
                xi[r] = N[r] / denom

        v = eta + 1
        for r in range(R):
            v -= xi[r] * Z[r]

        d.append(np.abs(u - u_1).sum() + abs(v - v_1))

    return u, v, d


def _pfqn_le_hessian(L: np.ndarray, N: np.ndarray, u0: np.ndarray) -> np.ndarray:
    """Compute Hessian of Gaussian."""
    M, R = L.shape
    Ntot = np.sum(N)
    hu = np.zeros((M - 1, M - 1))

    for i in range(M - 1):
        for j in range(M - 1):
            if i != j:
                hu[i, j] = -(Ntot + M) * u0[i] * u0[j]
                for r in range(R):
                    denom = np.dot(u0, L[:, r]) ** 2
                    if abs(denom) > 1e-14:
                        hu[i, j] += N[r] * L[i, r] * L[j, r] * u0[i] * u0[j] / denom
            else:
                hu[i, j] = (Ntot + M) * u0[i] * np.sum(_allbut(u0, i))
                for r in range(R):
                    denom = np.dot(u0, L[:, r]) ** 2
                    if abs(denom) > 1e-14:
                        L_allbut = L[np.arange(M) != i, r]
                        u_allbut = _allbut(u0, i)
                        hu[i, j] -= N[r] * L[i, r] * u0[i] * np.dot(u_allbut, L_allbut) / denom

    return hu


def _pfqn_le_hessianZ(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                       u: np.ndarray, v: float) -> np.ndarray:
    """Compute Hessian with think times."""
    K, R = L.shape
    Ntot = np.sum(N)
    A = np.zeros((K, K))

    csi = np.zeros(R)
    for r in range(R):
        denom = Z[r] + v * np.dot(u, L[:, r])
        if abs(denom) > 1e-14:
            csi[r] = N[r] / denom

    Lhat = np.zeros((K, R))
    for k in range(K):
        for r in range(R):
            Lhat[k, r] = Z[r] + v * L[k, r]

    eta = Ntot + K

    for i in range(K):
        for j in range(K):
            if i != j:
                A[i, j] = -eta * u[i] * u[j]
                for r in range(R):
                    if N[r] > 0:
                        A[i, j] += (csi[r] ** 2) * Lhat[i, r] * Lhat[j, r] * u[i] * u[j] / N[r]

    for i in range(K):
        row_sum = np.sum(A[i, :]) - A[i, i]
        A[i, i] = -row_sum

    A = A[:K-1, :K-1]

    # Extend to K dimensions
    A_ext = np.zeros((K, K))
    A_ext[:K-1, :K-1] = A
    A_ext[K-1, K-1] = 1.0

    for r in range(R):
        if N[r] > 0:
            A_ext[K-1, K-1] -= (csi[r] ** 2 / N[r]) * Z[r] * np.dot(u, L[:, r])
    A_ext[K-1, K-1] *= v

    for i in range(K - 1):
        A_ext[i, K-1] = 0
        for r in range(R):
            if N[r] > 0:
                term = (csi[r] ** 2 / N[r]) * Lhat[i, r] * np.dot(u, L[:, r]) - csi[r] * L[i, r]
                A_ext[i, K-1] += v * u[i] * term
        A_ext[K-1, i] = A_ext[i, K-1]

    return A_ext


def _simplex_fun(x: np.ndarray, L: np.ndarray, N: np.ndarray) -> float:
    """Evaluate simplex function."""
    M = len(x) + 1
    v = np.zeros(M)
    v[:M-1] = np.exp(x)
    v[M-1] = 1

    log_sum = 0
    for r in range(L.shape[1]):
        log_sum += N[r] * np.log(np.dot(v, L[:, r]))

    return np.exp(log_sum + np.sum(x) - (np.sum(N) + M) * np.log(np.sum(v)))


def pfqn_ls(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
            I: int = 100000) -> Tuple[float, float]:
    """
    Logistic sampling approximation for normalizing constant.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,) (default: zeros)
        I: Number of samples (default: 1e5)

    Returns:
        Tuple (Gn, lGn) - normalizing constant and its log

    Reference:
        G. Casale. Accelerating performance inference over closed systems by
        asymptotic methods. ACM SIGMETRICS 2017.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()

    # Filter out zero-demand stations
    Lsum = np.sum(L, axis=1)
    L = L[Lsum > 1e-4, :]

    M, R = L.shape

    if L.size == 0 or np.sum(L) < 1e-4 or len(N) == 0 or np.sum(N) == 0:
        lGn = -np.sum(_factln(N)) + np.sum(N * np.log(np.sum(Z) if Z is not None else 1))
        return np.exp(lGn), lGn

    if Z is None or len(Z) == 0 or np.sum(Z) < 1e-4:
        # No think time case
        umax, _ = _pfqn_le_fpi(L, N)
        A = _pfqn_le_hessian(L, N, umax)
        A = (A + A.T) / 2  # Symmetrize

        try:
            iA = np.linalg.inv(A)
            iA = (iA + iA.T) / 2  # Ensure symmetry
        except np.linalg.LinAlgError:
            # Singular matrix - return NaN
            return np.nan, np.nan

        x0 = np.log(umax[:M-1] / umax[M-1])

        # Generate samples
        try:
            samples = multivariate_normal.rvs(mean=x0, cov=iA, size=I)
            if I == 1:
                samples = samples.reshape(1, -1)
        except Exception:
            return np.nan, np.nan

        # Evaluate function at samples
        T = np.zeros(I)
        for i in range(I):
            T[i] = _simplex_fun(samples[i], L, N)

        # Compute density
        dpdf = multivariate_normal.pdf(samples, mean=x0, cov=iA)

        # Avoid division by zero
        valid = dpdf > 1e-300
        if not np.any(valid):
            return np.nan, np.nan

        lGn = _multinomialln(np.append(N, M-1)) + _factln(np.array([M-1]))[0] + np.log(np.mean(T[valid] / dpdf[valid]))
        Gn = np.exp(lGn)

    else:
        # With think time
        Z = np.asarray(Z, dtype=float).ravel()

        umax, vmax, _ = _pfqn_le_fpiZ(L, N, Z)
        A = _pfqn_le_hessianZ(L, N, Z, umax, vmax)
        A = (A + A.T) / 2  # Symmetrize

        try:
            iA = np.linalg.inv(A)
            iA = (iA + iA.T) / 2
        except np.linalg.LinAlgError:
            return np.nan, np.nan

        x0 = np.concatenate([np.log(umax[:M-1] / umax[M-1]), [np.log(vmax)]])

        try:
            samples = multivariate_normal.rvs(mean=x0, cov=iA, size=I)
            if I == 1:
                samples = samples.reshape(1, -1)
        except Exception:
            return np.nan, np.nan

        epsilon = 1e-10
        eN = epsilon * np.sum(N)
        eta = np.sum(N) + M * (1 + eN)
        K = M

        def h_func(x):
            """Evaluate integrand with think time."""
            xK = x[K-1]  # log(v)
            xm = x[:K-1]  # log(u_i/u_K) for i=1..K-1

            v_val = np.exp(xK)

            # Compute log of the sum term
            log_sum = 0
            for r in range(R):
                # L(K,:)*exp(x(K)) + Z
                base = L[K-1, r] * v_val + Z[r]
                # exp(x(1:K-1)) * (L(1:K-1,:)*exp(x(K)) + Z)
                inner = np.exp(xm) * (L[:K-1, r] * v_val + Z[r])
                log_sum += N[r] * np.log(base + np.sum(inner))

            result = -v_val + K * (1 + eN) * xK + log_sum + np.sum(xm)
            result -= eta * np.log(1 + np.sum(np.exp(xm)))

            return np.exp(result)

        T = np.zeros(I)
        for i in range(I):
            try:
                T[i] = h_func(samples[i])
            except (RuntimeWarning, FloatingPointError):
                T[i] = 0

        dpdf = multivariate_normal.pdf(samples, mean=x0, cov=iA)

        valid = (dpdf > 1e-300) & np.isfinite(T)
        if not np.any(valid):
            return np.nan, np.nan

        Gn = np.exp(-np.sum(gammaln(1 + N))) * np.mean(T[valid] / dpdf[valid])
        lGn = np.log(Gn) if Gn > 0 else -np.inf

    return Gn, lGn
