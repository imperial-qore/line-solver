"""
Knessl-Tier Asymptotic Expansion for Normalizing Constant.

Native Python implementation of the Knessl-Tier asymptotic expansion
method for computing normalizing constants of product-form queueing networks.

Key functions:
    pfqn_kt: Knessl-Tier asymptotic expansion

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_kt.m
    Knessl and Tier, "Asymptotic Expansions for Large Closed Queueing Networks"
"""

import numpy as np
from typing import Tuple, Optional

from .mva import pfqn_bs, pfqn_aql


# Small tolerance constant
FINE_TOL = 1e-14


def pfqn_kt(L: np.ndarray, N: np.ndarray,
            Z: Optional[np.ndarray] = None
            ) -> Tuple[float, float, np.ndarray, np.ndarray]:
    """
    Knessl-Tier asymptotic expansion for normalizing constant.

    Computes the normalizing constant using Knessl-Tier's asymptotic
    expansion, which is particularly accurate for large populations.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,), optional (default: zeros)

    Returns:
        Tuple of (G, lG, X, Q):
            G: Normalizing constant
            lG: Logarithm of normalizing constant
            X: System throughput (R,)
            Q: Mean queue lengths (M, R)

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_kt.m
    """
    if L is None or len(L) == 0 or N is None or len(N) == 0 or np.sum(N) == 0:
        return 1.0, 0.0, np.array([]), np.array([[]])

    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()

    if Z is None:
        Z = np.zeros(len(N))
    else:
        Z = np.asarray(Z, dtype=float).flatten()

    Morig, Rorig = L.shape

    # Handle self-looping customers (they would yield Uk=1)
    slcdemandfactor = 0.0

    if Rorig > 1:
        isslc = np.zeros(Rorig, dtype=bool)
        for r in range(Rorig):
            if np.count_nonzero(L[:, r]) == 1 and Z[r] == 0:
                ist = np.where(L[:, r] > 0)[0][0]
                # Replicate station for each job
                new_rows = np.tile(L[ist, :].reshape(1, -1), (int(N[r]), 1))
                L = np.vstack([L, new_rows])
                isslc[r] = True
                slcdemandfactor = N[r] * np.log(L[ist, r])

        # Remove self-looping classes
        keep_classes = ~isslc
        L = L[:, keep_classes]
        Z = Z[keep_classes]
        N = N[keep_classes]

    M, R = L.shape
    Ntot = int(np.sum(N))

    if Ntot == 0:
        return 1.0, 0.0, np.zeros(R), np.zeros((M, R))

    # Get throughput estimate
    if Ntot <= 4:
        result = pfqn_bs(L, N, Z)
        X, Q = result[0], result[1]  # XN, QN
    else:
        result = pfqn_aql(L, N, Z)
        X, Q = result[0], result[2]  # XN, QN (pfqn_aql still has old format)

    X = np.asarray(X).flatten()
    Q = np.atleast_2d(np.asarray(Q))

    # see _kb/03-api-layer.md for rationale
    u = X.astype(float).copy()
    Zc = Z.astype(float)
    Nc = N.astype(float)
    Uk = L @ u
    if np.max(Uk) >= 1:
        u = u * (1 - 1e-6) / np.max(Uk)
    converged = False
    for _ in range(200):
        Uk = L @ u
        D = 1.0 / (1.0 - Uk)
        g = u * (Zc + L.T @ D) - Nc
        if np.linalg.norm(g) <= 1e-12 * Ntot:
            converged = True
            break
        J = np.diag(Zc + L.T @ D) + u[:, None] * (L.T @ ((D ** 2)[:, None] * L))
        try:
            du = np.linalg.solve(J, -g)
        except np.linalg.LinAlgError:
            du = -np.linalg.lstsq(J, g, rcond=None)[0]
        alpha = 1.0
        while np.any(u + alpha * du <= 0) or np.max(L @ (u + alpha * du)) >= 1:
            alpha = alpha / 2
            if alpha < 1e-12:
                break
        if alpha < 1e-12:
            break
        u = u + alpha * du
    Uk = L @ u
    D = 1.0 / (1.0 - Uk)
    if converged and np.linalg.norm(u * (Zc + L.T @ D) - Nc) <= 1e-8 * Ntot:
        us = u  # exact saddle point
    else:
        us = X.astype(float).copy()  # fallback: AQL/BS throughput

    # Assemble the expansion at us
    Uk = L @ us
    D = 1.0 / np.maximum(FINE_TOL, 1.0 - Uk)
    H = np.diag(Nc / us ** 2) + L.T @ ((D ** 2)[:, None] * L)
    F = (Zc @ us - np.sum(np.log(np.maximum(FINE_TOL, 1.0 - Uk)))
         - Nc @ np.log(us))
    lG = (F - np.sum(np.log(us)) - (R / 2.0) * np.log(2 * np.pi)
          - 0.5 * np.log(np.linalg.det(H)) + slcdemandfactor)

    G = np.exp(lG)

    return G, lG, X, Q


__all__ = ['pfqn_kt']
