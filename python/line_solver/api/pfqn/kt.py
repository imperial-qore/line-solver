"""
Knessl-Tier Asymptotic Expansion for Normalizing Constant.

Native Python implementation of the Knessl-Tier asymptotic expansion
method for computing normalizing constants of product-form queueing networks.

Key functions:
    pfqn_kt: Knessl-Tier asymptotic expansion
    pfqn_bkt: the same, minus the Stirling remainder of each Laplaced class (BKT)

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_kt.m
    Knessl and Tier, "Asymptotic Expansions for Large Closed Queueing Networks"
"""

import numpy as np
from scipy.special import gammaln
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

    # A class with no jobs contributes a factor of 1 to G, but its saddle point
    # is u_r -> 0, where N_r*log(u_r) and N_r/u_r^2 are indeterminate and lG
    # comes back NaN. Solve the reduced model instead, as the self-looping
    # branch below already does for the classes it strips.
    if np.any(N == 0) and np.any(N > 0):
        keep = N > 0
        return pfqn_kt(L[:, keep], N[keep], Z[keep])

    Morig, Rorig = L.shape

    # Handle self-looping customers (they would yield Uk=1). Extracting
    # u_r^N_r from 1/(1-U_ist) exactly leaves L(ist,r)^N_r and raises that
    # station's factor to (1-V_ist)^-(1+N_r); classes looping at the SAME
    # station share one factor (1+sum N) and add the multinomial
    # (sum N)!/prod N_r!
    slcdemandfactor = 0.0

    if Rorig > 1:
        isslc = np.zeros(Rorig, dtype=bool)
        slcstation = np.zeros(Rorig, dtype=int)
        for r in range(Rorig):
            # detect on the original rows: appended copies must not mask a later class
            if np.count_nonzero(L[:Morig, r]) == 1 and Z[r] == 0:
                isslc[r] = True
                slcstation[r] = int(np.where(L[:Morig, r] > 0)[0][0])

        for ist in np.unique(slcstation[isslc]):
            grp = np.where(isslc & (slcstation == ist))[0]
            ntot = int(np.sum(N[grp]))
            slcdemandfactor += (float(np.dot(N[grp], np.log(L[ist, grp])))
                                + gammaln(ntot + 1.0)
                                - float(np.sum(gammaln(N[grp] + 1.0))))
            L = np.vstack([L, np.tile(L[ist, :].reshape(1, -1), (ntot, 1))])

        # Remove self-looping classes
        keep_classes = ~isslc
        L = L[:, keep_classes]
        Z = Z[keep_classes]
        N = N[keep_classes]

    M, R = L.shape
    Ntot = int(np.sum(N))

    if R == 0:
        # every class self-loops: the demand factors are the exact answer
        return float(np.exp(slcdemandfactor)), float(slcdemandfactor), np.array([]), np.array([[]])

    if Ntot == 0:
        return float(np.exp(slcdemandfactor)), float(slcdemandfactor), np.zeros(R), np.zeros((M, R))

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
    # log|H| via slogdet: det(H) of an R x R Hessian leaves double range well before
    # its logarithm does (it overflowed at R = 64, turning lG into -inf)
    _, logdetH = np.linalg.slogdet(H)
    lG = (F - np.sum(np.log(us)) - (R / 2.0) * np.log(2 * np.pi)
          - 0.5 * logdetH + slcdemandfactor)

    G = np.exp(lG)

    return G, lG, X, Q


def stirling_remainder(n):
    """s(N) = log(N!) - (N log N - N + log(2 pi N)/2), exactly, for N >= 1.

    s(1) = 1 - log(2 pi)/2, the constant pfqn_ble adds per station direction, and
    s(N) = 1/(12 N) + O(N^-2).
    """
    n = np.asarray(n, dtype=float)
    return gammaln(n + 1.0) - (n + 0.5) * np.log(n) + n - np.log(2 * np.pi) / 2


def pfqn_bkt(L: np.ndarray, N: np.ndarray,
              Z: Optional[np.ndarray] = None) -> Tuple[float, float]:
    r"""Knessl-Tier expansion corrected for the Stirling remainder (BKT).

    pfqn_kt extracts N from the generating function of G by steepest descent. On the
    demand-free integral the exact coefficient is [u^N] exp(Z u) = Z^N/N!, whereas the
    expansion returns N log Z - (N log N - N + log(2 pi N)/2), Stirling's approximation
    of log(N!) in place of log(N!). So KT lies ABOVE the exact value by the remainder
    s(N) per Laplaced class direction, and BKT subtracts sum_r s(N_r). The remainder
    is evaluated exactly from gammaln: truncating it at 1/(12 N) loses an order of
    magnitude (on the 1562 models of Cas17 sec5.3.1 the median \|error\| is 0.083 nats
    for KT, 1.9e-4 for the truncation and 1.4e-5 for the exact remainder).

    Only the classes pfqn_kt actually Laplaces are corrected: a class with no jobs is
    dropped by its recursion and a self-looping class (one nonzero demand and no think
    time) has its coefficient extracted exactly, so neither carries a remainder. The
    predicate here is pfqn_kt's own. See _kb/03-api-layer.md.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,), optional (default: zeros)

    Returns:
        (Gn, lGn), the normalizing constant and its logarithm.

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_bkt.m
        Knessl and Tier, IEEE Trans. Computers 41(4):480-488, 1992.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    Z = np.zeros(len(N)) if Z is None else np.asarray(Z, dtype=float).flatten()
    _, lG = pfqn_kt(L, N, Z)[:2]
    keep = N > 0
    Nc, Lc, Zc = N[keep], L[:, keep], Z[keep]
    if Nc.size > 1:
        isslc = (np.count_nonzero(Lc, axis=0) == 1) & (Zc == 0)
        Nc = Nc[~isslc]
    lG = float(lG - np.sum(stirling_remainder(Nc))) if Nc.size else float(lG)
    return float(np.exp(lG)), lG


__all__ = ['pfqn_kt', 'pfqn_bkt', 'stirling_remainder']
