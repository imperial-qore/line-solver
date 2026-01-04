"""
RECAL (REcursive CALculation) Method for Normalizing Constant.

Native Python implementation of the RECAL algorithm for computing
the normalizing constant of closed product-form queueing networks.

Key functions:
    pfqn_recal: RECAL method for normalizing constant computation

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_recal.m
"""

import numpy as np
from typing import Tuple, Optional

from .utils import multichoose, matchrow
from .replicas import pfqn_unique


def pfqn_recal(L: np.ndarray, N: np.ndarray,
               Z: Optional[np.ndarray] = None,
               m0: Optional[np.ndarray] = None) -> Tuple[float, float]:
    """
    RECAL (REcursive CALculation) method for normalizing constant.

    Computes the normalizing constant G(N) using the RECAL recursive
    method, which is efficient for networks with moderate population sizes.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,), optional (default: zeros)
        m0: Initial multiplicity vector (M,), optional (default: ones)

    Returns:
        Tuple of (G, lG):
            G: Normalizing constant
            lG: Logarithm of normalizing constant

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_recal.m
        Conway, A.E. and Georganas, N.D. "RECAL - A New Efficient
        Algorithm for the Exact Analysis of Multiple-Chain
        Closed Queueing Networks", JACM, 1986.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=int).flatten()

    M_original, R = L.shape

    # Default multiplicity
    if m0 is None:
        m0 = np.ones(M_original)
    else:
        m0 = np.asarray(m0, dtype=float).flatten()

    # Default think times
    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float).flatten()

    # Detect and consolidate replicated stations
    result = pfqn_unique(L)
    L = result.L
    mapping = result.mapping

    M = L.shape[0]

    # Combine user-provided m0 with detected multiplicity
    m0_combined = np.zeros(M)
    for i in range(M_original):
        m0_combined[mapping[i]] += m0[i]
    m0 = m0_combined

    Ntot = int(np.sum(N))

    if Ntot == 0:
        return 1.0, 0.0

    # Compute number of combinations
    from scipy.special import comb
    num_combinations = int(comb(Ntot + (M + 1) - 1, Ntot, exact=True))

    G_1 = np.ones(num_combinations)
    G = G_1.copy()

    I_1 = multichoose(M + 1, Ntot)
    n = 0

    for r in range(R):
        for nr in range(1, N[r] + 1):
            n += 1
            I = multichoose(M + 1, (Ntot + 1) - (n + 1))

            for i in range(I.shape[0]):
                m = I[i, :].copy()
                mZ = m[:M]

                # Find matching row for think time contribution
                idx = matchrow(I_1[:, :M], mZ)
                if idx > 0:
                    G[i] = Z[r] * G_1[idx - 1] / nr
                else:
                    G[i] = 0.0

                # Station contributions
                for jst in range(M):
                    m[jst] += 1
                    idx_station = matchrow(I_1, m)
                    if idx_station > 0:
                        G[i] += (m[jst] + m0[jst] - 1) * L[jst, r] * G_1[idx_station - 1] / nr
                    m[jst] -= 1

            I_1 = I
            G_1 = G.copy()

    G_final = G[0] if len(G) > 0 else 1.0
    lG = np.log(max(G_final, 1e-300))

    return G_final, lG


__all__ = ['pfqn_recal']
