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

    beta = np.zeros(R)
    for r in range(R):
        beta[r] = N[r] / Ntot if Ntot > 0 else 0

    # Get throughput estimate
    if Ntot <= 4:
        X, Q = pfqn_bs(L, N, Z)
    else:
        X, Q = pfqn_aql(L, N, Z)

    X = np.asarray(X).flatten()
    Q = np.atleast_2d(np.asarray(Q))

    # Compute covariance matrix C
    delta = np.eye(R)
    C = np.zeros((R, R))

    for s in range(R):
        for r in range(R):
            SK = 0.0
            for k in range(M):
                denom = 1 - np.sum(X * L[k, :])
                denom = max(FINE_TOL, denom) ** 2
                SK += X[s] * X[r] * L[k, s] * L[k, r] / denom

            C[s, r] = delta[s, r] * beta[s] + (1.0 / Ntot) * SK

    # Compute denominator product
    Den = 1.0
    for k in range(M):
        factor = 1 - np.sum(X * L[k, :])
        Den *= max(FINE_TOL, factor)

    # Compute log normalizing constant
    det_C = np.linalg.det(C)
    det_C = max(FINE_TOL, abs(det_C))

    # Avoid log of zero for throughput
    X_safe = np.maximum(X, FINE_TOL)

    lG = (np.log((2 * np.pi) ** (-R / 2) / np.sqrt(Ntot ** R * det_C)) +
          (-Ntot * np.dot(beta, np.log(X_safe))) -
          np.log(max(FINE_TOL, abs(Den))) +
          slcdemandfactor)

    G = np.exp(lG)

    return G, lG, X, Q


__all__ = ['pfqn_kt']
