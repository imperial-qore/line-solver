"""
Multi-server and Mixed Linearizer Algorithms.

Native Python implementations of the multi-server Linearizer algorithm
(Krzesinski/Conway/De Souza-Muntz) and the mixed open/closed linearizer.

Key functions:
    pfqn_linearizerms: Multi-server Linearizer
    pfqn_linearizermx: Mixed open/closed Linearizer

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_linearizerms.m
    Original MATLAB: matlab/src/api/pfqn/pfqn_linearizermx.m
    Conway, 1989, "Fast Approximate Solution of Queueing Networks
    with Multi-Server Chain-Dependent FCFS Queues"
"""

import numpy as np
from typing import Tuple, Optional

from .linearizer import SchedStrategy, pfqn_linearizer, pfqn_gflinearizer, pfqn_egflinearizer
from .mva import pfqn_bs
from .utils import oner


def pfqn_linearizerms(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                      nservers: np.ndarray,
                      type_sched: Optional[np.ndarray] = None,
                      tol: float = 1e-8,
                      maxiter: int = 1000
                      ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Multiserver Linearizer (Krzesinski/Conway/De Souza-Muntz).

    Extends the Linearizer algorithm to handle multi-server stations
    in product-form queueing networks.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,)
        nservers: Number of servers per station (M,)
        type_sched: Scheduling strategy per station (M,), optional (default: PS)
        tol: Convergence tolerance (default: 1e-8)
        maxiter: Maximum iterations (default: 1000)

    Returns:
        Tuple of (Q, U, R, C, X, totiter):
            Q: Mean queue lengths (M, R)
            U: Utilization (M, R)
            R: Residence times (M, R)
            C: Cycle times (R,)
            X: System throughput (R,)
            totiter: Total iterations performed

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_linearizerms.m
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()
    nservers = np.asarray(nservers, dtype=float).flatten()

    M, R = L.shape

    if type_sched is None:
        type_sched = np.full(M, SchedStrategy.PS)

    if len(Z) == 0:
        Z = np.zeros(R)

    max_servers = int(np.max(nservers))

    # Initialize Q, PB, P, Delta
    Q = np.zeros((M, R, 1 + R))
    PB = np.zeros((M, 1 + R))
    P = np.zeros((M, max_servers, 1 + R))
    Delta = np.zeros((M, R, R))

    # Initial estimates using Bard-Schweitzer
    for r in range(R):
        for s in range(R + 1):
            N_1 = oner(N, s) if s < R else N.copy()
            _, q = pfqn_bs(L, N_1, Z)
            Q[:, r, s] = q.flatten()

    for ist in range(M):
        for s in range(R + 1):
            N_1 = oner(N, s) if s < R else N.copy()
            pop = np.sum(N_1)
            if nservers[ist] > 1:
                for j in range(1, int(nservers[ist])):
                    P[ist, j, s] = 2 * np.sum(Q[ist, :, s]) / (pop * (pop + 1))

                if pop > nservers[ist] - 1:
                    PB[ist, s] = 2 * np.sum(Q[ist, :, s]) / (pop + 1 - nservers[ist]) / (pop * (pop + 1))
                else:
                    PB[ist, s] = 0.0

                P[ist, 0, s] = 1 - PB[ist, s] - np.sum(P[ist, 1:int(nservers[ist]), s])

    totiter = 0

    # Main loop (2 iterations)
    for I in range(2):
        for s in range(R + 1):
            N_1 = oner(N, s) if s < R else N.copy()

            # Core iteration
            Q[:, :, s], _, _, P[:, :, s], PB[:, s], iter_count = _core_ms(
                L, M, R, N_1, Z, nservers, Q[:, :, s], P[:, :, s], PB[:, s],
                Delta, type_sched, tol, maxiter - totiter
            )
            totiter += iter_count

        # Update Delta
        for ist in range(M):
            for r in range(R):
                for s in range(R):
                    Ns = oner(N, s)
                    if Ns[r] > 0 and N[r] > 0:
                        Delta[ist, r, s] = Q[ist, r, s] / Ns[r] - Q[ist, r, R] / N[r]
                    else:
                        Delta[ist, r, s] = 0.0

    # Final Core(N)
    Q_final, W, X, _, _, iter_count = _core_ms(
        L, M, R, N, Z, nservers, Q[:, :, R], P[:, :, R], PB[:, R],
        Delta, type_sched, tol, maxiter - totiter
    )
    totiter += iter_count

    # Compute performance metrics
    U = np.zeros((M, R))
    for ist in range(M):
        for r in range(R):
            if nservers[ist] == 1:
                U[ist, r] = X[r] * L[ist, r]
            else:
                U[ist, r] = X[r] * L[ist, r] / nservers[ist]

    Q_out = Q_final
    C = N / X - Z
    R_out = W

    return Q_out, U, R_out, C, X, totiter


def _core_ms(L, M, R, N_1, Z, nservers, Q, P, PB, Delta, type_sched, tol, maxiter):
    """Core iteration for multiserver linearizer."""
    max_servers = int(np.max(nservers))
    iter_count = 0
    hasConverged = False

    while not hasConverged:
        iter_count += 1
        Qlast = Q.copy()

        # Estimate
        Q_1, P_1, PB_1 = _estimate_ms(M, R, N_1, nservers, Q, P, PB, Delta)

        # Forward MVA
        Q, W, T, P, PB = _forward_mva_ms(L, M, R, N_1, Z, nservers, type_sched, Q_1, P_1, PB_1)

        if np.linalg.norm(Q - Qlast) < tol or iter_count >= maxiter:
            hasConverged = True

    return Q, W, T, P, PB, iter_count


def _estimate_ms(M, R, N_1, nservers, Q, P, PB, Delta):
    """Estimate populations for linearizer."""
    max_servers = int(np.max(nservers))
    P_1 = np.zeros((M, max_servers, 1 + R))
    PB_1 = np.zeros((M, 1 + R))
    Q_1 = np.zeros((M, R, 1 + R))

    for ist in range(M):
        if nservers[ist] > 1:
            for j in range(int(nservers[ist])):
                for s in range(R + 1):
                    P_1[ist, j, s] = P[ist, j]

            for s in range(R + 1):
                PB_1[ist, s] = PB[ist]

        for r in range(R):
            for s in range(R):
                Ns = oner(N_1, s)
                if N_1[r] > 0:
                    Q_1[ist, r, s] = Ns[r] * (Q[ist, r] / N_1[r] + Delta[ist, r, s])
                else:
                    Q_1[ist, r, s] = 0.0

    return Q_1, P_1, PB_1


def _forward_mva_ms(L, M, R, N_1, Z, nservers, type_sched, Q_1, P_1, PB_1):
    """Forward MVA step for multiserver linearizer."""
    max_servers = int(np.max(nservers))
    W = np.zeros((M, R))
    T = np.zeros(R)
    Q = np.zeros((M, R))
    P = np.zeros((M, max_servers))
    PB = np.zeros(M)

    for ist in range(M):
        for r in range(R):
            W[ist, r] = L[ist, r] / nservers[ist]
            if L[ist, r] == 0:
                continue

            if type_sched[ist] == SchedStrategy.FCFS:
                for s in range(R):
                    W[ist, r] += (L[ist, s] / nservers[ist]) * Q_1[ist, s, r]
            else:
                for s in range(R):
                    W[ist, r] += (L[ist, r] / nservers[ist]) * Q_1[ist, s, r]

            if nservers[ist] > 1:
                for j in range(int(nservers[ist]) - 1):
                    if type_sched[ist] == SchedStrategy.FCFS:
                        for s in range(R):
                            W[ist, r] += L[ist, s] * (nservers[ist] - 1 - j) * P_1[ist, j, r]
                    else:
                        for s in range(R):
                            W[ist, r] += L[ist, r] * (nservers[ist] - 1 - j) * P_1[ist, j, r]

    for r in range(R):
        denom = Z[r] + np.sum(W[:, r])
        if denom > 0:
            T[r] = N_1[r] / denom
        else:
            T[r] = 0.0

        for ist in range(M):
            Q[ist, r] = T[r] * W[ist, r]

    for ist in range(M):
        if nservers[ist] > 1:
            P[ist, :] = 0
            for j in range(1, int(nservers[ist])):
                for s in range(R):
                    P[ist, j] += L[ist, s] * T[s] * P_1[ist, j - 1, s] / j

    for ist in range(M):
        if nservers[ist] > 1:
            PB[ist] = 0
            for s in range(R):
                ns = int(nservers[ist])
                PB[ist] += L[ist, s] * T[s] * (PB_1[ist, s] + P_1[ist, ns - 1, s]) / ns

    for ist in range(M):
        if nservers[ist] > 1:
            P[ist, 0] = 1 - PB[ist]
            for j in range(1, int(nservers[ist])):
                P[ist, 0] -= P[ist, j]

    return Q, W, T, P, PB


def pfqn_linearizermx(lam: np.ndarray, L: np.ndarray, N: np.ndarray,
                      Z: np.ndarray, nservers: np.ndarray,
                      type_sched: Optional[np.ndarray] = None,
                      tol: float = 1e-8,
                      maxiter: int = 1000,
                      method: str = 'egflin'
                      ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Linearizer for mixed open/closed queueing networks.

    Handles networks with both open and closed classes using a
    decomposition approach.

    Args:
        lam: Arrival rate vector (R,) - inf for closed classes
        L: Service demand matrix (M x R)
        N: Population vector (R,) - inf for open classes
        Z: Think time vector (R,)
        nservers: Number of servers per station (M,)
        type_sched: Scheduling strategy per station (M,), optional
        tol: Convergence tolerance (default: 1e-8)
        maxiter: Maximum iterations (default: 1000)
        method: Linearizer variant ('lin', 'gflin', 'egflin', default: 'egflin')

    Returns:
        Tuple of (QN, UN, WN, CN, XN, totiter):
            QN: Mean queue lengths (M, R)
            UN: Utilization (M, R)
            WN: Waiting times (M, R)
            CN: Cycle times (R,)
            XN: System throughput (R,)
            totiter: Total iterations

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_linearizermx.m
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()
    lam = np.asarray(lam, dtype=float).flatten()
    nservers = np.asarray(nservers, dtype=float).flatten()

    M, R = L.shape

    if type_sched is None:
        type_sched = np.full(M, SchedStrategy.PS)

    # Handle NaN
    lam = np.nan_to_num(lam, nan=0.0)
    L = np.nan_to_num(L, nan=0.0)
    Z = np.nan_to_num(Z, nan=0.0)

    openClasses = np.where(np.isinf(N))[0]
    closedClasses = np.array([i for i in range(R) if i not in openClasses], dtype=int)

    # Initialize outputs
    XN = np.zeros(R)
    UN = np.zeros((M, R))
    WN = np.zeros((M, R))
    QN = np.zeros((M, R))
    CN = np.zeros(R)

    # Open class utilization
    for r in openClasses:
        for ist in range(M):
            UN[ist, r] = lam[r] * L[ist, r]
        XN[r] = lam[r]

    UNt = np.sum(UN, axis=1)

    if len(Z) == 0:
        Z = np.zeros(R)
    else:
        Z = np.sum(Z.reshape(-1, R) if Z.ndim > 1 else Z.reshape(1, -1), axis=0)

    # Effective demands for closed classes
    if len(closedClasses) > 0:
        Dc = L[:, closedClasses] / (1 - np.tile(UNt.reshape(-1, 1), (1, len(closedClasses))))

        if np.max(nservers) == 1:
            if method == 'lin' or method == 'default':
                QNc, UNc, WNc, CNc, XNc, totiter = pfqn_linearizer(
                    Dc, N[closedClasses], Z[closedClasses], type_sched, tol, maxiter
                )
            elif method == 'gflin':
                linAlpha = 2.0
                QNc, UNc, WNc, CNc, XNc, totiter = pfqn_gflinearizer(
                    Dc, N[closedClasses], Z[closedClasses], type_sched, tol, maxiter, linAlpha
                )
            elif method == 'egflin':
                alphaM = np.zeros(R)
                for r in closedClasses:
                    alphaM[r] = 0.6 + 1.4 * np.exp(-8 * np.exp(-0.8 * N[r]))
                QNc, UNc, WNc, CNc, XNc, totiter = pfqn_egflinearizer(
                    Dc, N[closedClasses], Z[closedClasses], type_sched, tol, maxiter,
                    alphaM[closedClasses]
                )
            else:
                QNc, UNc, WNc, CNc, XNc, totiter = pfqn_linearizer(
                    Dc, N[closedClasses], Z[closedClasses], type_sched, tol, maxiter
                )
        else:
            QNc, UNc, WNc, CNc, XNc, totiter = pfqn_linearizerms(
                Dc, N[closedClasses], Z[closedClasses], nservers, type_sched, tol, maxiter
            )

        XN[closedClasses] = XNc
        QN[:, closedClasses] = QNc
        WN[:, closedClasses] = WNc
        UN[:, closedClasses] = UNc
        CN[closedClasses] = CNc

        # Update closed class utilization
        for ist in range(M):
            for r in closedClasses:
                UN[ist, r] = XN[r] * L[ist, r]
    else:
        totiter = 0
        QNc = np.array([])

    # Open class metrics
    for ist in range(M):
        for r in openClasses:
            if len(QNc) == 0:
                WN[ist, r] = L[ist, r] / (1 - UNt[ist])
            else:
                WN[ist, r] = L[ist, r] * (1 + np.sum(QNc[ist, :])) / (1 - UNt[ist])

            QN[ist, r] = WN[ist, r] * XN[r]

    CN[openClasses] = np.sum(WN[:, openClasses], axis=0)

    return QN, UN, WN, CN, XN, totiter


def sprod(R: int, n: int) -> Tuple[int, np.ndarray, np.ndarray, np.ndarray]:
    """
    Initialize state product space iterator.

    Used for enumerating states in multi-class queueing networks.

    Args:
        R: Number of classes
        n: Total population constraint

    Returns:
        Tuple of (s, nvec, SD, D):
            s: Current state index
            nvec: Current state vector
            SD: Upper bounds
            D: Direction vector
    """
    nvec = np.zeros(R, dtype=int)
    nvec[0] = n
    SD = np.full(R, n, dtype=int)
    D = np.arange(R, dtype=int)
    s = 0
    return s, nvec, SD, D


def sprod_next(s: int, SD: np.ndarray, D: np.ndarray) -> Tuple[int, np.ndarray]:
    """
    Get next state in product space.

    Args:
        s: Current state index
        SD: Upper bounds
        D: Direction vector

    Returns:
        Tuple of (s, nvec) for next state, or (-1, nvec) if exhausted
    """
    R = len(SD)
    n = SD[0]

    if s < 0:
        return s, np.zeros(R, dtype=int)

    # Count total combinations
    from scipy.special import comb
    total = int(comb(n + R - 1, R - 1))

    s += 1
    if s >= total:
        return -1, np.zeros(R, dtype=int)

    # Convert s to nvec using stars and bars
    nvec = np.zeros(R, dtype=int)
    remaining = n
    for i in range(R - 1):
        for k in range(remaining + 1):
            c = int(comb(remaining - k + R - i - 2, R - i - 2))
            if s < c:
                nvec[i] = k
                remaining -= k
                break
            s -= c
    nvec[R - 1] = remaining

    return 0, nvec


def multinomialln(n: np.ndarray) -> float:
    """Compute log of multinomial coefficient."""
    from scipy.special import gammaln
    return float(gammaln(np.sum(n) + 1) - np.sum(gammaln(np.asarray(n) + 1)))


def pfqn_conwayms(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                   nservers: np.ndarray,
                   type_sched: Optional[np.ndarray] = None,
                   tol: float = 1e-8,
                   maxiter: int = 1000
                   ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Conway (1989) multiserver Linearizer approximation for FCFS queues.

    Implements the algorithm from Conway (1989), "Fast Approximate Solution
    of Queueing Networks with Multi-Server Chain-Dependent FCFS Queues".

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,)
        nservers: Number of servers per station (M,)
        type_sched: Scheduling strategy per station (M,), optional (default: FCFS)
        tol: Convergence tolerance (default: 1e-8)
        maxiter: Maximum iterations (default: 1000)

    Returns:
        Tuple of (Q, U, R, C, X, totiter):
            Q: Mean queue lengths (M, R)
            U: Utilization (M, R)
            R: Residence times (M, R)
            C: Cycle times (R,)
            X: System throughput (R,)
            totiter: Total iterations performed

    References:
        Conway, A. E., "Fast Approximate Solution of Queueing Networks with
        Multi-Server Chain-Dependent FCFS Queues", Performance Evaluation,
        Vol. 8, 1989, pp. 141-159.
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()
    nservers = np.asarray(nservers, dtype=float).flatten()

    M, R = L.shape

    if type_sched is None:
        type_sched = np.full(M, SchedStrategy.FCFS)

    if len(Z) == 0:
        Z = np.zeros(R)
    Z = np.sum(Z.reshape(-1, R) if Z.ndim > 1 else Z.reshape(1, -1), axis=0)

    max_servers = int(np.max(nservers))

    # Initialize Q, PB, P, Delta indexed by population reduction
    Q = np.zeros((M, R, 1 + R))
    PB = np.zeros((M, 1 + R))
    P = np.zeros((M, max_servers, 1 + R))
    Delta = np.zeros((M, R, R))

    # Initialize queue lengths
    for ist in range(M):
        for r in range(R):
            for s in range(R + 1):
                N_1 = oner(N, s) if s < R else N.copy()
                Q[ist, r, s] = N_1[r] / M

    # Initialize probabilities
    for ist in range(M):
        for s in range(R + 1):
            N_1 = oner(N, s) if s < R else N.copy()
            pop = np.sum(N_1)
            if nservers[ist] > 1:
                for j in range(1, int(nservers[ist])):
                    if pop * (pop + 1) > 0:
                        P[ist, j, s] = 2 * np.sum(Q[ist, :, s]) / (pop * (pop + 1))

                denom = (pop + 1 - nservers[ist]) * pop * (pop + 1)
                if denom > 0:
                    PB[ist, s] = 2 * np.sum(Q[ist, :, s]) / denom
                else:
                    PB[ist, s] = 0.0

                P[ist, 0, s] = max(0, 1 - PB[ist, s] - np.sum(P[ist, 1:int(nservers[ist]), s]))

    totiter = 0

    # Main loop (2 iterations for linearizer)
    for I in range(2):
        for s in range(R + 1):
            N_1 = oner(N, s) if s < R else N.copy()

            # Core iteration
            Q[:, :, s], W_temp, T_temp, P[:, :, s], PB[:, s], iter_count = _conway_core(
                L, M, R, N_1, Z, nservers, Q[:, :, s], P[:, :, s], PB[:, s],
                Delta, W_temp if 'W_temp' in dir() else L.copy(), type_sched, tol, maxiter - totiter
            )
            totiter += iter_count

        # Update Delta
        for ist in range(M):
            for r in range(R):
                for s in range(R):
                    Ns = oner(N, s)
                    if N[s] > 2 and N[r] > 0 and Ns[r] > 0:
                        Delta[ist, r, s] = Q[ist, r, s] / Ns[r] - Q[ist, r, R] / N[r]
                    else:
                        Delta[ist, r, s] = 0.0

    # Final Core(N)
    Q_final, W, X, P_final, PB_final, iter_count = _conway_core(
        L, M, R, N, Z, nservers, Q[:, :, R], P[:, :, R], PB[:, R],
        Delta, L.copy(), type_sched, tol, maxiter - totiter
    )
    totiter += iter_count

    # Compute performance metrics
    U = np.zeros((M, R))
    for ist in range(M):
        for r in range(R):
            if nservers[ist] == 1:
                U[ist, r] = X[r] * L[ist, r]
            else:
                U[ist, r] = X[r] * L[ist, r] / nservers[ist]

    Q_out = Q_final
    C = N / np.maximum(X, 1e-12) - Z
    R_out = W

    return Q_out, U, R_out, C, X, totiter


def _conway_core(L, M, R, N_1, Z, nservers, Q, P, PB, Delta, W, type_sched, tol, maxiter):
    """Core iteration for Conway multiserver linearizer."""
    max_servers = int(np.max(nservers))
    hasConverged = False
    iter_count = 0

    while not hasConverged:
        Qlast = Q.copy()

        # Estimate
        Q_1, P_1, PB_1, T_1 = _conway_estimate(M, R, N_1, nservers, Q, P, PB, Delta, W)

        # Forward MVA
        Q, W, T, P, PB = _conway_forward_mva(L, M, R, N_1, Z, nservers, type_sched, Q_1, P_1, PB_1, T_1)

        if np.linalg.norm(Q - Qlast) < tol or iter_count >= maxiter:
            hasConverged = True

        iter_count += 1

    return Q, W, T, P, PB, iter_count


def _conway_estimate(M, R, N_1, nservers, Q, P, PB, Delta, W):
    """Estimate populations for Conway linearizer."""
    max_servers = int(np.max(nservers))
    P_1 = np.zeros((M, max_servers, 1 + R))
    PB_1 = np.zeros((M, 1 + R))
    Q_1 = np.zeros((M, R, 1 + R))
    T_1 = np.zeros((R, 1 + R))

    for ist in range(M):
        if nservers[ist] > 1:
            for j in range(int(nservers[ist])):
                for s in range(R + 1):
                    P_1[ist, j, s] = P[ist, j]
            for s in range(R + 1):
                PB_1[ist, s] = PB[ist]

        for r in range(R):
            for s in range(R):
                Ns = oner(N_1, s)
                if N_1[r] > 0:
                    Q_1[ist, r, s] = Ns[r] * (Q[ist, r] / N_1[r] + Delta[ist, r, s])
                else:
                    Q_1[ist, r, s] = 0.0

    # Compute T_1
    for r in range(R):
        for s in range(R):
            Nr = oner(N_1, r)
            for ist in range(M):
                if W[ist, s] > 0 and N_1[s] > 0:
                    T_1[s, r] = Nr[s] * (Q[ist, s] / N_1[s] + Delta[ist, r, s]) / W[ist, s]
                    break

    return Q_1, P_1, PB_1, T_1


def _conway_forward_mva(L, M, R, N_1, Z, nservers, type_sched, Q_1, P_1, PB_1, T_1):
    """Forward MVA step for Conway multiserver linearizer."""
    max_servers = int(np.max(nservers))
    W = np.zeros((M, R))
    T = np.zeros(R)
    Q = np.zeros((M, R))
    P = np.zeros((M, max_servers))
    PB = np.zeros(M)
    XR = np.zeros((M, R))
    XE = np.zeros((M, R, R))

    mu = 1.0 / np.maximum(L, 1e-12)

    # Compute F matrix for each class
    F = []
    for r in range(R):
        F_r = np.zeros((M, R))
        for ist in range(M):
            den = np.dot(L[ist, :], T_1[:, r])
            if den > 0:
                for c in range(R):
                    F_r[ist, c] = T_1[c, r] * L[ist, c] / den
        F.append(F_r)

    # Compute XR (multi-server specific)
    for ist in range(M):
        for r in range(R):
            if nservers[ist] > 1:
                ns = int(nservers[ist])
                XR[ist, r] = 0.0
                C_val = 0.0

                s, nvec, SD, D = sprod(R, ns)
                while s >= 0:
                    Nr = oner(N_1, r)
                    if np.all(nvec <= Nr):
                        Ai = np.exp(multinomialln(nvec) + np.dot(nvec, np.log(np.maximum(F[r][ist, :], 1e-300))))
                        C_val += Ai
                        mu_sum = np.dot(mu[ist, :], nvec)
                        if mu_sum > 0:
                            XR[ist, r] += Ai / mu_sum

                    s, nvec = sprod_next(s, SD, D)

                if C_val > 0:
                    XR[ist, r] /= C_val

    # Compute XE
    for ist in range(M):
        for r in range(R):
            if nservers[ist] > 1:
                ns = int(nservers[ist])
                for c in range(R):
                    XE[ist, r, c] = 0.0
                    Cx = 0.0

                    s, nvec, SD, D = sprod(R, ns)
                    while s >= 0:
                        Nr = oner(N_1, r)
                        if np.all(nvec <= Nr) and nvec[c] >= 1:
                            Aix = np.exp(multinomialln(nvec) + np.dot(nvec, np.log(np.maximum(F[r][ist, :], 1e-300))))
                            Cx += Aix
                            mu_sum = np.dot(mu[ist, :], nvec)
                            if mu_sum > 0:
                                XE[ist, r, c] += Aix / mu_sum

                        s, nvec = sprod_next(s, SD, D)

                    if Cx > 0:
                        XE[ist, r, c] /= Cx

    # Compute residence time
    for ist in range(M):
        for r in range(R):
            if nservers[ist] == 1:
                if type_sched[ist] == SchedStrategy.FCFS:
                    W[ist, r] = L[ist, r]
                    for c in range(R):
                        W[ist, r] += L[ist, c] * Q_1[ist, c, r]
                else:
                    W[ist, r] = L[ist, r]
                    for c in range(R):
                        W[ist, r] += L[ist, r] * Q_1[ist, c, r]
            else:
                W[ist, r] = L[ist, r] + PB_1[ist, r] * XR[ist, r]
                for c in range(R):
                    W[ist, r] += XE[ist, r, c] * (Q_1[ist, c, r] - L[ist, c] * T_1[c, r])

    # Compute throughputs and queue lengths
    for r in range(R):
        denom = Z[r] + np.sum(W[:, r])
        if denom > 0:
            T[r] = N_1[r] / denom
        else:
            T[r] = 0.0

        for ist in range(M):
            Q[ist, r] = T[r] * W[ist, r]

    # Compute marginal probabilities
    for ist in range(M):
        if nservers[ist] > 1:
            ns = int(nservers[ist])
            P[ist, :] = 0
            for j in range(1, ns):
                for c in range(R):
                    if j > 0:
                        P[ist, j] += L[ist, c] * T[c] * P_1[ist, j - 1, c] / j

    for ist in range(M):
        if nservers[ist] > 1:
            ns = int(nservers[ist])
            PB[ist] = 0
            for c in range(R):
                PB[ist] += L[ist, c] * T[c] * (PB_1[ist, c] + P_1[ist, ns - 1, c]) / ns

    for ist in range(M):
        if nservers[ist] > 1:
            P[ist, 0] = max(0, 1 - PB[ist])
            for j in range(1, int(nservers[ist])):
                P[ist, 0] = max(0, P[ist, 0] - P[ist, j])

    return Q, W, T, P, PB


__all__ = [
    'pfqn_linearizerms',
    'pfqn_linearizermx',
    'pfqn_conwayms',
]
