"""
Load-Dependent MVA for Mixed Open/Closed Networks.

Native Python implementations of load-dependent MVA algorithms for
mixed queueing networks with limited load dependence.

Key functions:
    pfqn_mvaldmx: Load-dependent MVA for mixed networks
    pfqn_mvaldmx_ec: Effective capacity computation
    pfqn_mvaldms: Multi-server wrapper for mvaldmx

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_mvaldmx*.m
"""

import numpy as np
from typing import Tuple, Optional

from .schmidt import pprod, hashpop
from .utils import oner


def pfqn_mvaldmx_ec(lam: np.ndarray, D: np.ndarray, mu: np.ndarray
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute effective capacity terms for MVALDMX solver.

    Calculates the effective capacity E, E', and EC terms needed for
    load-dependent MVA with limited load dependence.

    Args:
        lam: Arrival rate vector (R,)
        D: Service demand matrix (M x R)
        mu: Load-dependent rate matrix (M x Nt)

    Returns:
        Tuple of (EC, E, Eprime, Lo):
            EC: Effective capacity matrix (M x Nt)
            E: E-function values (M x (1+Nt))
            Eprime: E-prime function values (M x (1+Nt))
            Lo: Open class load vector (M,)

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_mvaldmx_ec.m
    """
    mu = np.atleast_2d(np.asarray(mu, dtype=float))
    D = np.atleast_2d(np.asarray(D, dtype=float))
    lam = np.asarray(lam, dtype=float).flatten()

    M, Nt_orig = mu.shape

    # Compute open class load
    Lo = np.zeros(M)
    for ist in range(M):
        Lo[ist] = np.dot(lam, D[ist, :])

    # Find limited load dependence level b for each station
    b = np.zeros(M, dtype=int)
    for ist in range(M):
        # Find first index where mu equals the final value
        for k in range(Nt_orig):
            if mu[ist, k] == mu[ist, -1]:
                b[ist] = k + 1  # 1-based index
                break
        if b[ist] == 0:
            b[ist] = Nt_orig

    Nt = Nt_orig
    max_b = int(np.max(b))

    # Extend mu if needed
    mu_ext = np.zeros((M, Nt + 1 + max_b))
    mu_ext[:, :Nt] = mu
    for k in range(Nt, Nt + 1 + max_b):
        mu_ext[:, k] = mu[:, -1]

    # C = 1/mu (service time)
    C = 1.0 / np.maximum(mu_ext, 1e-12)

    # Initialize outputs
    EC = np.zeros((M, Nt))
    E = np.zeros((M, 1 + Nt))
    Eprime = np.zeros((M, 1 + Nt))

    for ist in range(M):
        b_i = int(b[ist])

        E1 = np.zeros(1 + Nt)
        E2 = np.zeros(1 + Nt)
        E3 = np.zeros(1 + Nt)
        E2prime = np.zeros(1 + Nt)

        F2 = np.zeros((1 + Nt, max(1, b_i - 1)))
        F3 = np.zeros((1 + Nt, max(1, b_i - 1)))
        F2prime = np.zeros((1 + Nt, max(1, b_i - 1)))

        for n in range(Nt + 1):
            if n >= b_i:
                # Simple formula for n >= b
                denom = 1 - Lo[ist] * C[ist, b_i - 1]
                if abs(denom) > 1e-12:
                    E[ist, n] = 1.0 / (denom ** (n + 1))
                else:
                    E[ist, n] = 1e12

                Eprime[ist, n] = C[ist, b_i - 1] * E[ist, n]

            else:  # n <= b_i - 1
                # Compute E1
                if n == 0:
                    denom = 1 - Lo[ist] * C[ist, b_i - 1]
                    if abs(denom) > 1e-12:
                        E1[n] = 1.0 / denom
                    else:
                        E1[n] = 1e12

                    for j in range(b_i - 1):
                        E1[n] *= C[ist, j] / C[ist, b_i - 1]
                else:
                    E1[n] = (1.0 / (1 - Lo[ist] * C[ist, b_i - 1])) * \
                            (C[ist, b_i - 1] / C[ist, n - 1]) * E1[n - 1]

                # Compute F2
                for n0 in range(b_i - 1):
                    if n0 == 0:
                        F2[n, n0] = 1.0
                    else:
                        idx = min(n + n0 - 1, C.shape[1] - 1)
                        F2[n, n0] = ((n + n0) / n0) * Lo[ist] * C[ist, idx] * F2[n, n0 - 1]

                # Compute E2
                E2[n] = np.sum(F2[n, :b_i - 1])

                # Compute F3
                for n0 in range(b_i - 1):
                    if n == 0 and n0 == 0:
                        F3[n, n0] = 1.0
                        for j in range(b_i - 1):
                            F3[n, n0] *= C[ist, j] / C[ist, b_i - 1]
                    elif n > 0 and n0 == 0:
                        F3[n, n0] = (C[ist, b_i - 1] / C[ist, n - 1]) * F3[n - 1, 0]
                    else:
                        F3[n, n0] = ((n + n0) / n0) * Lo[ist] * C[ist, b_i - 1] * F3[n, n0 - 1]

                # Compute E3
                E3[n] = np.sum(F3[n, :b_i - 1])

                # Compute F2prime
                for n0 in range(b_i - 1):
                    if n0 == 0:
                        idx = min(n, C.shape[1] - 1)
                        F2prime[n, n0] = C[ist, idx]
                    else:
                        idx = min(n + n0, C.shape[1] - 1)
                        F2prime[n, n0] = ((n + n0) / n0) * Lo[ist] * C[ist, idx] * F2prime[n, n0 - 1]

                # Compute E2prime
                E2prime[n] = np.sum(F2prime[n, :b_i - 1])

                # Final E and Eprime
                E[ist, n] = E1[n] + E2[n] - E3[n]

                if n < b_i - 1:
                    Eprime[ist, n] = C[ist, b_i - 1] * E1[n] + E2prime[n] - C[ist, b_i - 1] * E3[n]
                else:
                    Eprime[ist, n] = C[ist, b_i - 1] * E[ist, n]

        # Compute EC
        for n in range(1, Nt + 1):
            if E[ist, n - 1] > 0:
                EC[ist, n - 1] = C[ist, n - 1] * E[ist, n] / E[ist, n - 1]
            else:
                EC[ist, n - 1] = 0.0

    return EC, E, Eprime, Lo


def pfqn_mvaldmx(lam: np.ndarray, D: np.ndarray, N: np.ndarray,
                 Z: np.ndarray, mu: Optional[np.ndarray] = None,
                 S: Optional[np.ndarray] = None
                 ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, np.ndarray]:
    """
    Load-dependent MVA for mixed open/closed networks with limited load dependence.

    Implements the MVALDMX algorithm for analyzing mixed queueing networks
    with load-dependent service rates using limited load dependence.

    Args:
        lam: Arrival rate vector (R,) - non-zero for open classes
        D: Service demand matrix (M x R)
        N: Population vector (R,) - inf for open classes
        Z: Think time vector (R,)
        mu: Load-dependent rate matrix (M x Nt), optional
        S: Number of servers per station (M,), optional

    Returns:
        Tuple of (XN, QN, UN, CN, lGN, Pc):
            XN: System throughput (R,)
            QN: Mean queue lengths (M, R)
            UN: Utilization (M, R)
            CN: Cycle times (M, R)
            lGN: Logarithm of normalizing constant
            Pc: Marginal queue-length probabilities (M, 1+Ntot, prod(1+Nc))

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_mvaldmx.m
    """
    D = np.atleast_2d(np.asarray(D, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()
    lam = np.asarray(lam, dtype=float).flatten()

    M, R = D.shape

    # Identify open and closed classes
    openClasses = np.where(np.isinf(N))[0]
    closedClasses = np.array([i for i in range(R) if i not in openClasses], dtype=int)

    Nct = int(np.sum(N[np.isfinite(N)]))

    # Default mu and S
    if mu is None:
        mu = np.ones((M, Nct))
        S = np.ones(M)
    elif S is None:
        S = np.ones(M)

    mu = np.atleast_2d(np.asarray(mu, dtype=float))
    S = np.asarray(S, dtype=float).flatten()

    # Extend mu if needed
    if mu.shape[1] < Nct:
        mu_ext = np.zeros((M, Nct + 1))
        mu_ext[:, :mu.shape[1]] = mu
        mu_ext[:, mu.shape[1]:] = mu[:, -1].reshape(-1, 1)
        mu = mu_ext

    # Initialize outputs
    XN = np.zeros(R)
    UN = np.zeros((M, R))
    CN = np.zeros((M, R))
    QN = np.zeros((M, R))
    lGN = 0.0

    # Extend mu by one column
    mu = np.hstack([mu, mu[:, -1:]])

    # Compute effective capacity terms
    EC, E, Eprime = pfqn_mvaldmx_ec(lam, D, mu)[:3]

    C = len(closedClasses)
    if C == 0:
        # Pure open network - use standard formulas
        for r in openClasses:
            XN[r] = lam[r]
            for ist in range(M):
                rho = np.sum([lam[s] * D[ist, s] for s in openClasses])
                if rho < 1:
                    QN[ist, r] = lam[r] * D[ist, r] / (1 - rho)
                    CN[ist, r] = D[ist, r] / (1 - rho)
                    UN[ist, r] = lam[r] * D[ist, r]

        return XN, QN, UN, CN, lGN, np.array([])

    Dc = D[:, closedClasses]
    Nc = N[closedClasses].astype(int)
    Zc = Z[closedClasses]

    # Precomputed products for hashing
    prods = np.zeros(C, dtype=int)
    for r in range(C):
        prods[r] = int(np.prod(Nc[:r] + 1))

    total_states = int(np.prod(Nc + 1))

    # Initialize
    nvec = pprod(Nc)
    Pc = np.zeros((M, 1 + int(np.sum(Nc)), total_states))
    x = np.zeros((C, total_states))
    w = np.zeros((M, C, total_states))

    # Initialize Pc(0|0) = 1
    hnvec_init = hashpop(nvec, Nc, C, prods)
    for ist in range(M):
        Pc[ist, 0, hnvec_init - 1] = 1.0

    u = np.zeros((M, C))

    # Population recursion
    while np.all(nvec >= 0):
        hnvec = hashpop(nvec, Nc, C, prods)
        nc = int(np.sum(nvec))

        for ist in range(M):
            for c in range(C):
                if nvec[c] > 0:
                    hnvec_c = hashpop(oner(nvec, c), Nc, C, prods)

                    # Compute mean residence times
                    for n in range(1, nc + 1):
                        if n - 1 < EC.shape[1] and Pc[ist, n - 1, hnvec_c - 1] > 0:
                            w[ist, c, hnvec - 1] += Dc[ist, c] * n * EC[ist, n - 1] * Pc[ist, n - 1, hnvec_c - 1]

        # Compute throughput
        for c in range(C):
            denom = Zc[c] + np.sum(w[:, c, hnvec - 1])
            if denom > 0:
                x[c, hnvec - 1] = nvec[c] / denom
            else:
                x[c, hnvec - 1] = 0.0

        # Update Pc
        for ist in range(M):
            for n in range(1, nc + 1):
                for c in range(C):
                    if nvec[c] > 0:
                        hnvec_c = hashpop(oner(nvec, c), Nc, C, prods)
                        if n - 1 < EC.shape[1]:
                            Pc[ist, n, hnvec - 1] += Dc[ist, c] * EC[ist, n - 1] * x[c, hnvec - 1] * Pc[ist, n - 1, hnvec_c - 1]

            Pc[ist, 0, hnvec - 1] = max(1e-12, 1 - np.sum(Pc[ist, 1:nc + 1, hnvec - 1]))

        nvec = pprod(nvec, Nc)

    # Final performance metrics at Nc
    hnvec = hashpop(Nc, Nc, C, prods)

    # Utilization
    # MATLAB: for n=1:sum(Nc), using Eprime(ist, n) which is 1-based
    # In 0-based terms: n from 0 to sum(Nc)-1
    for c in range(C):
        hnvec_c = hashpop(oner(Nc, c), Nc, C, prods)
        for ist in range(M):
            u[ist, c] = 0.0
            for n in range(int(np.sum(Nc))):  # Fixed: was range(sum(Nc)+1), should be range(sum(Nc))
                if n < E.shape[1] and E[ist, n] > 0:
                    u[ist, c] += Dc[ist, c] * x[c, hnvec - 1] * Eprime[ist, n] / E[ist, n] * Pc[ist, n, hnvec_c - 1]

    # Throughput
    XN[closedClasses] = x[:, hnvec - 1]

    # Utilization
    UN[:, closedClasses] = u

    # Response time
    CN[:, closedClasses] = w[:, :, hnvec - 1]

    # Queue-length
    for c in range(C):
        QN[:, closedClasses[c]] = XN[closedClasses[c]] * CN[:, closedClasses[c]]

    # Open class metrics
    for r in openClasses:
        XN[r] = lam[r]
        for ist in range(M):
            QN[ist, r] = 0.0
            for n in range(int(np.sum(Nc)) + 1):
                if n < EC.shape[1]:
                    QN[ist, r] += lam[r] * D[ist, r] * (n + 1) * EC[ist, n] * Pc[ist, n, hnvec - 1]

            if lam[r] > 0:
                CN[ist, r] = QN[ist, r] / lam[r]
            else:
                CN[ist, r] = 0.0

            UN[ist, r] = 0.0
            for n in range(int(np.sum(Nc)) + 1):
                if n < E.shape[1] and E[ist, n] > 0:
                    UN[ist, r] += lam[r] * Eprime[ist, n] / E[ist, n] * Pc[ist, n, hnvec - 1]

    # Final Pc
    Pc_final = Pc[:, :, hnvec - 1]

    return XN, QN, UN, CN, lGN, Pc_final


def pfqn_mvaldms(lam: np.ndarray, D: np.ndarray, N: np.ndarray,
                 Z: np.ndarray, S: np.ndarray
                 ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Load-dependent MVA for multiserver mixed networks.

    Wrapper for pfqn_mvaldmx that adjusts utilizations to account for
    multi-server stations.

    Args:
        lam: Arrival rate vector (R,)
        D: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,)
        S: Number of servers per station (M,)

    Returns:
        Tuple of (XN, QN, UN, CN, lGN):
            XN: System throughput (R,)
            QN: Mean queue lengths (M, R)
            UN: Utilization (M, R) - adjusted for multiservers
            CN: Cycle times (M, R)
            lGN: Logarithm of normalizing constant

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_mvaldms.m
    """
    D = np.atleast_2d(np.asarray(D, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()
    lam = np.asarray(lam, dtype=float).flatten()
    S = np.asarray(S, dtype=float).flatten()

    M, R = D.shape

    # Compute closed population
    Nct = int(np.sum(N[np.isfinite(N)]))

    if len(Z) == 0 or Z is None:
        Z = np.zeros(R)

    # Handle pure open networks (Nct=0) using M/M/k formulas directly
    # pfqn_mvaldmx cannot handle empty mu matrix
    if Nct == 0:
        # All classes are open - use M/M/k formulas
        openClasses = np.where(np.isinf(N))[0]
        XN = np.zeros(R)
        QN = np.zeros((M, R))
        UN = np.zeros((M, R))
        CN = np.zeros((M, R))
        lGN = np.nan  # Normalizing constant not defined for open networks

        for r in openClasses:
            XN[r] = lam[r]
            for ist in range(M):
                if D[ist, r] > 0 and S[ist] > 0:
                    # M/M/k: rho = lambda * D / S
                    rho = lam[r] * D[ist, r] / S[ist]
                    UN[ist, r] = rho
                    if rho < 1:
                        k = int(S[ist])
                        if k == 1:
                            # M/M/1: Q = rho/(1-rho)
                            QN[ist, r] = rho / (1 - rho)
                        else:
                            # M/M/k approximation using Erlang-C
                            from ..qsys import qsys_mmk
                            mu_rate = 1.0 / D[ist, r] if D[ist, r] > 0 else float('inf')
                            result = qsys_mmk(lam[r], mu_rate, k)
                            QN[ist, r] = result.get('L', rho / (1 - rho))
                        CN[ist, r] = QN[ist, r] / lam[r] if lam[r] > 0 else 0

        return XN, QN, UN, CN, lGN

    # Build load-dependent rate matrix
    mu = np.ones((M, Nct))
    for ist in range(M):
        for n in range(Nct):
            mu[ist, n] = min(n + 1, S[ist])

    # Call mvaldmx
    XN, QN, _, CN, lGN, _ = pfqn_mvaldmx(lam, D, N, Z, mu, S)

    # Identify open and closed classes
    openClasses = np.where(np.isinf(N))[0]
    closedClasses = np.array([i for i in range(R) if i not in openClasses], dtype=int)

    # Adjust utilization for multiservers
    UN = np.zeros((M, R))
    for r in closedClasses:
        for ist in range(M):
            UN[ist, r] = XN[r] * D[ist, r] / S[ist]

    for r in openClasses:
        for ist in range(M):
            UN[ist, r] = lam[r] * D[ist, r] / S[ist]

    return XN, QN, UN, CN, lGN


__all__ = [
    'pfqn_mvaldmx',
    'pfqn_mvaldmx_ec',
    'pfqn_mvaldms',
]
