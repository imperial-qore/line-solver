"""
Schmidt's Exact MVA for General Scheduling Disciplines.

Native Python implementation of Schmidt's exact MVA algorithm for
product-form queueing networks with general scheduling disciplines
(PS, FCFS, INF) including multi-server stations.

Key functions:
    pfqn_schmidt: Schmidt's exact MVA algorithm

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_schmidt.m
"""

import numpy as np
from typing import Tuple, Optional, Union, List
from dataclasses import dataclass
from enum import IntEnum, Enum

from .utils import oner


class SchedStrategy(IntEnum):
    """Scheduling strategies for queueing stations (integer values for MVA)."""
    INF = 0      # Infinite server (delay)
    FCFS = 1     # First Come First Serve
    PS = 2       # Processor Sharing
    LCFS = 3     # Last Come First Serve


def _normalize_sched(sched: np.ndarray) -> np.ndarray:
    """
    Normalize scheduling strategy array to integer values.

    Handles both string-based and integer-based SchedStrategy values.

    Args:
        sched: Array of scheduling strategies (string or int)

    Returns:
        Array of integer scheduling strategy values
    """
    result = np.zeros(len(sched), dtype=int)

    # Mapping from string to int
    str_to_int = {
        'INF': SchedStrategy.INF,
        'FCFS': SchedStrategy.FCFS,
        'PS': SchedStrategy.PS,
        'LCFS': SchedStrategy.LCFS,
    }

    for i, s in enumerate(sched):
        if isinstance(s, (int, np.integer)):
            result[i] = int(s)
        elif isinstance(s, str):
            result[i] = str_to_int.get(s, SchedStrategy.FCFS)
        elif hasattr(s, 'value'):
            # Handle Enum types
            val = s.value
            if isinstance(val, (int, np.integer)):
                result[i] = int(val)
            elif isinstance(val, str):
                result[i] = str_to_int.get(val, SchedStrategy.FCFS)
            else:
                result[i] = SchedStrategy.FCFS
        else:
            result[i] = SchedStrategy.FCFS

    return result


@dataclass
class SchmidtResult:
    """Result from Schmidt's exact MVA."""
    XN: np.ndarray  # System throughput (R,)
    QN: np.ndarray  # Mean queue lengths (M, R)
    UN: np.ndarray  # Utilization (M, R)
    CN: np.ndarray  # Cycle times (M, R)


def pprod(N: np.ndarray, Nmax: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Generate population vectors for MVA recursion.

    When called with one argument, initializes to zero vector.
    When called with two arguments, increments to next population vector.
    Returns -1 when iteration is complete.

    Args:
        N: Current population vector or max population
        Nmax: Maximum population per class (if incrementing)

    Returns:
        Next population vector, or array of -1 if done
    """
    if Nmax is None:
        # Initialize to zero vector
        return np.zeros(len(N), dtype=int)

    # Increment to next population
    N = np.asarray(N, dtype=int).copy()
    Nmax = np.asarray(Nmax, dtype=int)

    R = len(N)
    for r in range(R):
        if N[r] < Nmax[r]:
            N[r] += 1
            return N
        else:
            N[r] = 0

    # All combinations exhausted
    return -np.ones(R, dtype=int)


def hashpop(nvec: np.ndarray, Nc: np.ndarray, C: int, prods: np.ndarray) -> int:
    """
    Hash population vector to linear index.

    Args:
        nvec: Population vector
        Nc: Maximum population per class
        C: Number of classes
        prods: Precomputed products for hashing

    Returns:
        Linear index (1-based for MATLAB compatibility)
    """
    nvec = np.asarray(nvec, dtype=int)
    idx = 1  # 1-based indexing
    for r in range(C):
        idx += int(nvec[r] * prods[r])
    return idx


def pfqn_schmidt(D: np.ndarray, N: np.ndarray, S: np.ndarray,
                 sched: np.ndarray, v: Optional[np.ndarray] = None
                 ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Schmidt's exact MVA for networks with general scheduling disciplines.

    Implements Schmidt's exact Mean Value Analysis algorithm for product-form
    queueing networks with PS, FCFS, or INF scheduling disciplines, including
    support for multi-server stations.

    Args:
        D: Service demand matrix (M x R)
        N: Population vector (R,)
        S: Number of servers per station (M,) or (M x R)
        sched: Scheduling discipline per station (M,) - SchedStrategy values
        v: Visit ratio matrix (M x R), optional (default: ones)

    Returns:
        Tuple of (XN, QN, UN, CN):
            XN: System throughput (M x R) - same per station for closed networks
            QN: Mean queue lengths (M, R)
            UN: Utilization (M, R) - returns zeros (not computed)
            CN: Cycle times / response times (M, R)

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_schmidt.m
    """
    D = np.atleast_2d(np.asarray(D, dtype=float))
    N = np.asarray(N, dtype=int).flatten()
    S = np.asarray(S, dtype=int)
    sched = _normalize_sched(np.atleast_1d(sched))

    M, R = D.shape
    closedClasses = np.arange(R)

    # Initialize outputs
    XN = np.zeros(R)
    UN = np.zeros((M, R))
    CN = np.zeros((M, R))
    QN = np.zeros((M, R))

    # Default visit ratios
    if v is None:
        v = np.ones((M, R))
    else:
        v = np.atleast_2d(np.asarray(v, dtype=float))

    # Compute pure service times from demands
    S_pure = D / np.maximum(v, 1e-12)

    C = len(closedClasses)
    Dc = D[:, closedClasses]
    Nc = N[closedClasses]

    # Precomputed products for hashing
    prods = np.zeros(C, dtype=int)
    for r in range(C):
        prods[r] = int(np.prod(Nc[:r] + 1))

    # Initialize population recursion
    kvec = pprod(Nc)

    # Initialize L (mean queue-length) and Pc (state probabilities)
    total_states = int(np.prod(Nc + 1))
    L = [np.zeros((R, total_states)) for _ in range(M)]
    Pc = [None] * M

    # Handle S dimension
    if S.ndim == 1:
        S_mat = np.tile(S.reshape(-1, 1), (1, R))
    else:
        S_mat = S

    for ist in range(M):
        if sched[ist] == SchedStrategy.INF:
            pass  # No Pc needed
        elif sched[ist] == SchedStrategy.PS:
            if not np.all(S_mat[ist, :] == 1):
                Pc[ist] = np.zeros((1 + int(np.sum(Nc)), total_states))
        elif sched[ist] == SchedStrategy.FCFS:
            if np.all(D[ist, :] == D[ist, 0]):  # Class-independent
                if not np.all(S_mat[ist, :] == 1):
                    Pc[ist] = np.zeros((1 + int(np.sum(Nc)), total_states))
            else:  # Class-dependent
                Pc[ist] = np.zeros((total_states, total_states))

    # Per-station throughput and waiting time arrays
    x = np.zeros((M, C, total_states))
    w = np.zeros((M, C, total_states))

    # Initialize Pc(0|0) = 1
    hkvec_init = hashpop(kvec, Nc, C, prods)
    for ist in range(M):
        if Pc[ist] is not None:
            Pc[ist][0, hkvec_init - 1] = 1.0

    u = np.zeros((M, C))

    # Population recursion
    while np.all(kvec >= 0) and np.all(kvec <= Nc):
        nc = int(np.sum(kvec))
        hkvec = hashpop(kvec, Nc, C, prods)

        for ist in range(M):
            ns = int(S_mat[ist, 0]) if S_mat.ndim == 2 else int(S[ist])

            for c in range(C):
                if kvec[c] > 0:
                    hkvec_c = hashpop(oner(kvec, c), Nc, C, prods)

                    if sched[ist] == SchedStrategy.INF:
                        w[ist, c, hkvec - 1] = D[ist, c]

                    elif sched[ist] == SchedStrategy.PS:
                        if ns == 1:
                            totalQL = np.sum(L[ist][:, hkvec_c - 1])
                            w[ist, c, hkvec - 1] = Dc[ist, c] * (1 + totalQL)
                        else:
                            totalQL = np.sum(L[ist][:, hkvec_c - 1])
                            w[ist, c, hkvec - 1] = (Dc[ist, c] / ns) * (1 + totalQL)
                            if Pc[ist] is not None:
                                for j in range(1, ns):
                                    w[ist, c, hkvec - 1] += (ns - 1 - (j - 1)) * Pc[ist][j, hkvec_c - 1] * (Dc[ist, c] / ns)

                    elif sched[ist] == SchedStrategy.FCFS:
                        if np.all(D[ist, :] == D[ist, 0]):  # Product-form case
                            if ns == 1:
                                totalQL = np.sum(L[ist][:, hkvec_c - 1])
                                w[ist, c, hkvec - 1] = Dc[ist, c] * (1 + totalQL)
                            else:
                                totalQL = np.sum(L[ist][:, hkvec_c - 1])
                                w[ist, c, hkvec - 1] = (Dc[ist, c] / ns) * (1 + totalQL)
                                if Pc[ist] is not None:
                                    for j in range(1, ns):
                                        w[ist, c, hkvec - 1] += (ns - 1 - (j - 1)) * Pc[ist][j, hkvec_c - 1] * (Dc[ist, c] / ns)
                        else:  # Class-dependent
                            if ns == 1:
                                totalQL = np.sum(L[ist][:, hkvec_c - 1])
                                w[ist, c, hkvec - 1] = Dc[ist, c] * (1 + totalQL)
                            else:
                                # Multi-server FCFS class-dependent case
                                nvec_inner = pprod(kvec)
                                while np.all(nvec_inner >= 0):
                                    if nvec_inner[c] > 0:
                                        hnvec_c = hashpop(oner(nvec_inner, c), Nc, C, prods)
                                        n_sum = int(np.sum(nvec_inner))
                                        if n_sum <= ns:
                                            Bcn = S_pure[ist, c]
                                        else:
                                            sumVal = np.dot(nvec_inner, S_pure[ist, :])
                                            Bcn = S_pure[ist, c] + max(0, n_sum - ns) / max(ns * (n_sum - 1), 1e-12) * (sumVal - S_pure[ist, c])

                                        if Pc[ist] is not None:
                                            w[ist, c, hkvec - 1] += Bcn * Pc[ist][hnvec_c - 1, hkvec_c - 1]

                                    nvec_inner = pprod(nvec_inner, kvec)

        # Compute throughputs
        for c in range(C):
            denom = 0.0
            for ist in range(M):
                denom += v[ist, c] * w[ist, c, hkvec - 1]

            for ist in range(M):
                if denom > 0:
                    x[ist, c, hkvec - 1] = v[ist, c] * kvec[c] / denom
                else:
                    x[ist, c, hkvec - 1] = 0.0

        # Update queue lengths
        for ist in range(M):
            for c in range(C):
                L[ist][c, hkvec - 1] = x[ist, c, hkvec - 1] * w[ist, c, hkvec - 1]

            ns = int(S_mat[ist, 0]) if S_mat.ndim == 2 else int(S[ist])

            if sched[ist] == SchedStrategy.PS and ns > 1:
                if Pc[ist] is not None:
                    for n in range(1, min(ns, int(np.sum(kvec))) + 1):
                        for c in range(C):
                            if kvec[c] > 0:
                                hkvec_c = hashpop(oner(kvec, c), Nc, C, prods)
                                Pc[ist][n, hkvec - 1] += Dc[ist, c] * (1.0 / n) * x[ist, c, hkvec - 1] * Pc[ist][n - 1, hkvec_c - 1]

                    Pc[ist][0, hkvec - 1] = max(1e-12, 1 - np.sum(Pc[ist][1:min(ns, int(np.sum(kvec))) + 1, hkvec - 1]))

            elif sched[ist] == SchedStrategy.FCFS:
                if np.all(D[ist, :] == D[ist, 0]) and ns > 1:
                    if Pc[ist] is not None:
                        for n in range(1, min(ns, int(np.sum(kvec)))):
                            for c in range(C):
                                if kvec[c] > 0:
                                    hkvec_c = hashpop(oner(kvec, c), Nc, C, prods)
                                    Pc[ist][n, hkvec - 1] += Dc[ist, c] * (1.0 / n) * x[ist, c, hkvec - 1] * Pc[ist][n - 1, hkvec_c - 1]

                        Pc[ist][0, hkvec - 1] = max(1e-12, 1 - np.sum(Pc[ist][1:min(ns, int(np.sum(kvec))) + 1, hkvec - 1]))

        kvec = pprod(kvec, Nc)

    # Final throughput computation
    hkvec = hashpop(Nc, Nc, C, prods)
    for c in range(C):
        totalRT = 0.0
        for ist in range(M):
            totalRT += w[ist, c, hkvec - 1]

        if totalRT > 0:
            XN[c] = Nc[c] / totalRT
        else:
            XN[c] = 0.0

    # Replicate throughput for all stations
    if M > 1:
        XN = np.tile(XN.reshape(1, -1), (M, 1))
    else:
        XN = XN.reshape(1, -1)

    # Response times
    CN[:, closedClasses] = w[:, :, hkvec - 1]

    # Queue lengths
    for ist in range(M):
        QN[ist, closedClasses] = L[ist][closedClasses, hkvec - 1]

    return XN, QN, UN, CN


def pfqn_schmidt_ext(D: np.ndarray, N: np.ndarray, S: np.ndarray,
                      sched: np.ndarray, v: Optional[np.ndarray] = None
                      ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extended Schmidt MVA algorithm with queue-aware alpha corrections.

    A queue-aware version of the Schmidt algorithm that precomputes alpha values
    for improved accuracy in networks with class-dependent FCFS scheduling.

    Reference:
        R. Schmidt, "An approximate MVA algorithm for exponential,
        class-dependent multiple server stations," Performance Evaluation,
        vol. 29, no. 4, pp. 245-254, 1997.

    Args:
        D: Service demand matrix (M x R)
        N: Population vector (R,)
        S: Number of servers per station (M,) or (M x R)
        sched: Scheduling discipline per station (M,) - SchedStrategy values
        v: Visit ratio matrix (M x R), optional (default: ones)

    Returns:
        Tuple of (XN, QN, UN, CN):
            XN: System throughput (M x R)
            QN: Mean queue lengths (M, R)
            UN: Utilization (M, R)
            CN: Cycle times / response times (M, R)
    """
    D = np.atleast_2d(np.asarray(D, dtype=float))
    N = np.asarray(N, dtype=int).flatten()
    S = np.asarray(S, dtype=int)
    sched = _normalize_sched(np.atleast_1d(sched))

    M, R = D.shape
    closedClasses = np.arange(R)

    # Initialize outputs
    XN = np.zeros(R)
    UN = np.zeros((M, R))
    CN = np.zeros((M, R))
    QN = np.zeros((M, R))

    # Default visit ratios
    if v is None:
        v = np.ones((M, R))
    else:
        v = np.atleast_2d(np.asarray(v, dtype=float))

    C = len(closedClasses)
    Dc = D[:, closedClasses]
    Nc = N[closedClasses]

    # Handle S dimension
    if S.ndim == 1:
        S_mat = np.tile(S.reshape(-1, 1), (1, R))
    else:
        S_mat = S

    # Precompute alphas for class-dependent FCFS stations
    alphas = [[None for _ in range(R)] for _ in range(M)]
    for ist in range(M):
        if sched[ist] == SchedStrategy.FCFS:
            # Check if class-dependent
            class_independent = np.all(D[ist, :] == D[ist, 0])
            if not class_independent:
                for r in range(R):
                    # Create modified problem with tagged customer
                    D_mod = np.zeros((M, R + 1))
                    N_mod = np.zeros(R + 1, dtype=int)

                    for k in range(R):
                        N_mod[k] = N[k] - 1 if k == r else N[k]
                        for j in range(M):
                            D_mod[j, k] = D[j, k]
                            if ist == j:
                                D_mod[j, R] = D[j, r]
                            else:
                                D_mod[j, R] = 0

                    N_mod[R] = 1
                    sched_mod = np.append(sched, SchedStrategy.FCFS)

                    # Run basic Schmidt to get alpha values
                    _, _, U_alpha, _, = pfqn_schmidt(D_mod, N_mod, S, sched_mod)
                    alphas[ist][r] = U_alpha

    # Precomputed products for hashing
    prods = np.zeros(C, dtype=int)
    for r in range(C):
        prods[r] = int(np.prod(Nc[:r] + 1))

    # Initialize population recursion
    kvec = pprod(Nc)

    # Initialize L (mean queue-length) and Pc (state probabilities)
    total_states = int(np.prod(Nc + 1))
    L = [np.zeros((R, total_states)) for _ in range(M)]
    Pc = [None] * M

    for ist in range(M):
        if sched[ist] == SchedStrategy.INF:
            pass  # No Pc needed
        elif sched[ist] == SchedStrategy.PS:
            if not np.all(S_mat[ist, :] == 1):
                Pc[ist] = np.zeros((1 + int(np.sum(Nc)), total_states))
        elif sched[ist] == SchedStrategy.FCFS:
            if np.all(D[ist, :] == D[ist, 0]):  # Class-independent
                if not np.all(S_mat[ist, :] == 1):
                    Pc[ist] = np.zeros((1 + int(np.sum(Nc)), total_states))
            else:  # Class-dependent
                Pc[ist] = np.zeros((total_states, total_states))

    # Per-station throughput and waiting time arrays
    x = np.zeros((C, total_states))
    w = np.zeros((M, C, total_states))

    # Initialize Pc(0|0) = 1
    hkvec_init = hashpop(kvec, Nc, C, prods)
    for ist in range(M):
        if Pc[ist] is not None:
            Pc[ist][0, hkvec_init - 1] = 1.0

    # Population recursion
    while np.all(kvec >= 0) and np.all(kvec <= Nc):
        hkvec = hashpop(kvec, Nc, C, prods)

        for ist in range(M):
            ns = int(S_mat[ist, 0]) if S_mat.ndim == 2 else int(S[ist])

            for c in range(C):
                if kvec[c] > 0:
                    hkvec_c = hashpop(oner(kvec, c), Nc, C, prods)

                    if sched[ist] == SchedStrategy.INF:
                        w[ist, c, hkvec - 1] = D[ist, c]

                    elif sched[ist] == SchedStrategy.PS:
                        if ns == 1:
                            totalQL = np.sum(L[ist][:, hkvec_c - 1])
                            w[ist, c, hkvec - 1] = Dc[ist, c] * (1 + totalQL)
                        else:
                            totalQL = np.sum(L[ist][:, hkvec_c - 1])
                            w[ist, c, hkvec - 1] = (Dc[ist, c] / ns) * (1 + totalQL)
                            if Pc[ist] is not None:
                                for j in range(1, ns):
                                    w[ist, c, hkvec - 1] += (ns - 1 - (j - 1)) * Pc[ist][j, hkvec_c - 1] * (Dc[ist, c] / ns)

                    elif sched[ist] == SchedStrategy.FCFS:
                        class_independent = np.all(D[ist, :] == D[ist, 0])

                        if not class_independent:
                            if ns == 1:
                                totalQL = np.sum(L[ist][:, hkvec_c - 1])
                                w[ist, c, hkvec - 1] = Dc[ist, c] * (1 + totalQL)
                            else:
                                nvec_inner = pprod(kvec)
                                while np.all(nvec_inner >= 0):
                                    if nvec_inner[c] > 0:
                                        hnvec_c = hashpop(oner(nvec_inner, c), Nc, C, prods)

                                        # Use extended Bcn with alpha
                                        if Nc[c] > 1 and alphas[ist][c] is not None:
                                            Bcn = _get_bcn_ext(alphas[ist][c], D, ist, c, nvec_inner, C, ns, R)
                                        else:
                                            Bcn = _get_bcn(D, ist, c, nvec_inner, C, ns)

                                        if Pc[ist] is not None:
                                            w[ist, c, hkvec - 1] += Bcn * Pc[ist][hnvec_c - 1, hkvec_c - 1]

                                    nvec_inner = pprod(nvec_inner, kvec)
                        else:
                            # Class-independent: treat like PS
                            if ns == 1:
                                totalQL = np.sum(L[ist][:, hkvec_c - 1])
                                w[ist, c, hkvec - 1] = Dc[ist, c] * (1 + totalQL)
                            else:
                                totalQL = np.sum(L[ist][:, hkvec_c - 1])
                                w[ist, c, hkvec - 1] = (Dc[ist, c] / ns) * (1 + totalQL)
                                if Pc[ist] is not None:
                                    for j in range(1, ns):
                                        w[ist, c, hkvec - 1] += (ns - 1 - (j - 1)) * Pc[ist][j, hkvec_c - 1] * (Dc[ist, c] / ns)

        # Compute throughputs
        for c in range(C):
            sumW = np.sum(w[:M, c, hkvec - 1])
            if sumW > 0:
                x[c, hkvec - 1] = kvec[c] / sumW
            else:
                x[c, hkvec - 1] = 0.0

        # Update queue lengths
        for ist in range(M):
            for c in range(C):
                L[ist][c, hkvec - 1] = x[c, hkvec - 1] * w[ist, c, hkvec - 1]

            ns = int(S_mat[ist, 0]) if S_mat.ndim == 2 else int(S[ist])

            if sched[ist] == SchedStrategy.PS and ns > 1:
                if Pc[ist] is not None:
                    for n in range(1, min(ns, int(np.sum(kvec))) + 1):
                        for c in range(C):
                            if kvec[c] > 0:
                                hkvec_c = hashpop(oner(kvec, c), Nc, C, prods)
                                Pc[ist][n, hkvec - 1] += Dc[ist, c] * (1.0 / n) * x[c, hkvec - 1] * Pc[ist][n - 1, hkvec_c - 1]

                    Pc[ist][0, hkvec - 1] = max(1e-12, 1 - np.sum(Pc[ist][1:min(ns, int(np.sum(kvec))) + 1, hkvec - 1]))

            elif sched[ist] == SchedStrategy.FCFS:
                class_independent = np.all(D[ist, :] == D[ist, 0])

                if not class_independent:
                    nvec_inner = pprod(kvec)
                    nvec_inner = pprod(nvec_inner, kvec)
                    sumOfAllProbs = 0.0
                    while np.all(nvec_inner >= 0):
                        hnvec = hashpop(nvec_inner, Nc, C, prods)

                        prob = 0.0
                        for r in range(C):
                            if nvec_inner[r] > 0:
                                hnvec_c = hashpop(oner(nvec_inner, r), Nc, C, prods)
                                hkvec_c = hashpop(oner(kvec, r), Nc, C, prods)

                                if Nc[r] > 1 and alphas[ist][r] is not None:
                                    Bcn = _get_bcn_ext(alphas[ist][r], D, ist, r, nvec_inner, C, ns, R)
                                else:
                                    Bcn = _get_bcn(D, ist, r, nvec_inner, C, ns)

                                capacity_inv = 1.0 / np.sum(nvec_inner)
                                x_ir = x[r, hkvec - 1]
                                prob_c = Pc[ist][hnvec_c - 1, hkvec_c - 1] if Pc[ist] is not None else 0.0
                                classProb = Bcn * capacity_inv * x_ir * prob_c
                                prob += classProb

                        if Pc[ist] is not None:
                            Pc[ist][hnvec - 1, hkvec - 1] = prob
                        sumOfAllProbs += prob
                        nvec_inner = pprod(nvec_inner, kvec)

                    if Pc[ist] is not None:
                        Pc[ist][0, hkvec - 1] = max(1e-12, 1 - sumOfAllProbs)
                else:
                    if ns > 1 and Pc[ist] is not None:
                        for n in range(1, min(ns, int(np.sum(kvec)))):
                            for c in range(C):
                                if kvec[c] > 0:
                                    hkvec_c = hashpop(oner(kvec, c), Nc, C, prods)
                                    Pc[ist][n, hkvec - 1] += Dc[ist, c] * (1.0 / n) * x[c, hkvec - 1] * Pc[ist][n - 1, hkvec_c - 1]

                        Pc[ist][0, hkvec - 1] = max(1e-12, 1 - np.sum(Pc[ist][1:min(ns, int(np.sum(kvec))) + 1, hkvec - 1]))

        kvec = pprod(kvec, Nc)

    # Extract final results
    hkvec_final = hashpop(Nc, Nc, C, prods)

    # Throughput
    XN[closedClasses] = x[:C, hkvec_final - 1]
    if M > 1:
        XN = np.tile(XN.reshape(1, -1), (M, 1))
    else:
        XN = XN.reshape(1, -1)

    # Utilization
    for m in range(M):
        ns = int(S_mat[m, 0]) if S_mat.ndim == 2 else int(S[m])
        for c in range(C):
            UN[m, c] = (D[m, c] * XN[0, c]) / ns

    # Response time
    CN[:, closedClasses] = w[:M, :C, hkvec_final - 1]

    # Queue length
    for ist in range(M):
        QN[ist, closedClasses] = L[ist][closedClasses, hkvec_final - 1]

    return XN, QN, UN, CN


def _get_bcn(D: np.ndarray, i: int, c: int, nvec: np.ndarray, C: int, ns: int) -> float:
    """Standard Bcn calculation for Schmidt algorithm."""
    Bcn = D[i, c]
    n_sum = int(np.sum(nvec))
    if n_sum > 1:
        eps_val = 1e-12
        sumVal = np.dot(nvec, D[i, :C])
        Bcn = Bcn + max(0, n_sum - ns) / max(ns * (n_sum - 1), eps_val) * (sumVal - D[i, c])
    return Bcn


def _get_bcn_ext(u: np.ndarray, D: np.ndarray, i: int, c: int, nvec: np.ndarray,
                  C: int, ns: int, R: int) -> float:
    """Extended Bcn with alpha corrections for Schmidt algorithm."""
    # u is the utilization matrix from the alpha computation (M x R+1)
    # The last column (R+1) is the tagged customer utilization
    if u is None:
        return _get_bcn(D, i, c, nvec, C, ns)

    weightedProb = 0.0
    totalNonPinnedTime = np.sum(u[i, :C]) - u[i, C] if u.shape[1] > C else 0.0

    if totalNonPinnedTime > 0:
        for s in range(C):
            if D[i, s] > 0:
                prob = u[i, s] / totalNonPinnedTime
                weightedProb += prob / D[i, s]

    if weightedProb > 0:
        meanInterdepartureTime = 1.0 / (ns * weightedProb)
    else:
        meanInterdepartureTime = 0.0

    Bcn = D[i, c]
    n_sum = int(np.sum(nvec))

    if n_sum > 1:
        Bcn = Bcn + max(0, n_sum - ns) * meanInterdepartureTime

    if np.isnan(Bcn) or np.isinf(Bcn):
        Bcn = 0.0

    return Bcn


__all__ = [
    'pfqn_schmidt',
    'pfqn_schmidt_ext',
    'SchmidtResult',
    'pprod',
    'hashpop',
]
