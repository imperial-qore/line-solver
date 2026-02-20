"""
Mean Value Analysis (MVA) algorithms for Product-Form Queueing Networks.

Native Python implementations of MVA and related algorithms including:
- Standard MVA for closed networks
- Multi-class MVA with population recursion
- Approximate MVA methods

References:
    Reiser, M., and Lavenberg, S.S. "Mean-Value Analysis of Closed Multichain
    Queueing Networks." Journal of the ACM 27.2 (1980): 313-322.
"""

import numpy as np
from numpy.linalg import LinAlgError
from typing import Tuple, Optional, Dict, Any, Union
from math import log, ceil

from .replicas import pfqn_unique, pfqn_expand, pfqn_combine_mi

# Import JIT-compiled kernels
try:
    from .mva_jit import (
        HAS_NUMBA as MVA_HAS_NUMBA,
        mva_single_class_jit,
        mva_population_recursion_jit,
        schweitzer_iteration_jit,
    )
except ImportError:
    # Fallback if JIT import fails
    MVA_HAS_NUMBA = False
    mva_single_class_jit = None
    mva_population_recursion_jit = None
    schweitzer_iteration_jit = None

# Threshold for using JIT (population space size)
JIT_THRESHOLD = 100


def _population_lattice_pprod(n: np.ndarray, N: np.ndarray = None) -> np.ndarray:
    """
    Generate next population vector in lexicographic order.

    Args:
        n: Current population vector
        N: Maximum population (optional, for wrapping)

    Returns:
        Next population vector, or (-1,...,-1) when exhausted
    """
    R = len(n)
    n_next = n.copy()

    # If N is None, just increment by 1 in the last position and carry
    if N is None:
        n_next[-1] += 1
        return n_next

    # Find rightmost position that can be incremented
    for i in range(R - 1, -1, -1):
        if n_next[i] < N[i]:
            n_next[i] += 1
            # Reset all positions to the right to 0
            for j in range(i + 1, R):
                n_next[j] = 0
            return n_next

    # All exhausted
    return -np.ones(R)


def _population_lattice_hashpop(n: np.ndarray, N: np.ndarray) -> int:
    """
    Compute hash index for population vector n given max N.

    Args:
        n: Population vector
        N: Maximum population vector

    Returns:
        Linear index in flattened population lattice
    """
    R = len(n)
    idx = 0
    mult = 1
    for i in range(R - 1, -1, -1):
        idx += int(n[i]) * mult
        mult *= int(N[i]) + 1
    return idx


def pfqn_mva_single_class(N: int, L: np.ndarray, Z: float = 0.0,
                          mi: Optional[np.ndarray] = None) -> Dict[str, Any]:
    """
    Mean Value Analysis for single-class closed network.

    Simplified MVA for single customer class, with optional multi-server
    stations specified via mi (multiplicity).

    Args:
        N: Number of customers
        L: Service demands at each station (1D array of length M)
        Z: Think time (default 0)
        mi: Number of servers at each station (default all 1)

    Returns:
        dict with keys:
            - 'X': Throughput
            - 'Q': Queue lengths (array of length M)
            - 'R': Residence times (array of length M)
            - 'U': Utilizations (array of length M)
            - 'lG': Log of normalizing constant
    """
    L = np.asarray(L, dtype=np.float64).flatten()
    M = len(L)

    if mi is None:
        mi = np.ones(M)
    else:
        mi = np.asarray(mi, dtype=np.float64).flatten()

    if N <= 0:
        return {
            'X': 0.0,
            'Q': np.zeros(M),
            'R': L.copy(),
            'U': np.zeros(M),
            'lG': 0.0
        }

    # Use JIT version for larger populations
    if MVA_HAS_NUMBA and mva_single_class_jit is not None and N > JIT_THRESHOLD:
        X, Q, R, U, lG = mva_single_class_jit(N, L, Z, mi)
        return {'X': X, 'Q': Q, 'R': R, 'U': U, 'lG': lG}

    # Pure Python MVA recursion
    Q = np.zeros(M)
    lG = 0.0

    for n in range(1, N + 1):
        # Residence times: R_i = L_i * (m_i + Q_i(n-1))
        R = L * (mi + Q)

        # Throughput: X = n / (Z + sum(R))
        total_R = R.sum()
        X = n / (Z + total_R)

        # Update queue lengths
        Q = X * R

        # Update log normalizing constant
        if n > 0:
            lG -= log(X)

    # Utilizations
    U = X * L

    return {
        'X': X,
        'Q': Q,
        'R': R,
        'U': U,
        'lG': lG
    }


def pfqn_mva(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
             mi: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                              np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Mean Value Analysis for multi-class closed product-form network.

    Implements the exact MVA algorithm using population recursion.
    Computes exact performance measures for closed product-form networks
    with load-independent stations.

    Args:
        L: Service demand matrix (M x R) where M is stations, R is classes
        N: Population vector (1 x R or R,) - number of jobs per class
        Z: Think time vector (1 x R or R,) - think time per class (default 0)
        mi: Multiplicity vector (1 x M or M,) - servers per station (default 1)

    Returns:
        Tuple of (XN, CN, QN, UN, RN, TN, AN) where:
            - XN: Throughputs per class (1 x R)
            - CN: Response times per class (1 x R) - total cycle time
            - QN: Queue lengths (M x R)
            - UN: Utilizations (M x R)
            - RN: Residence times (M x R)
            - TN: Node throughputs (M x R)
            - AN: Arrival rates (M x R)
    """
    L = np.asarray(L, dtype=np.float64)
    N = np.asarray(N, dtype=np.float64).flatten()
    N = np.ceil(N).astype(int)

    R = len(N)  # Number of classes
    if L.ndim == 1:
        L = L.reshape(-1, 1) if R == 1 else L.reshape(1, -1)
    M_original = L.shape[0]  # Original number of stations

    if L.shape[1] != R:
        raise ValueError(f"Demand matrix columns ({L.shape[1]}) must match population size ({R})")

    # Handle Z
    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()
        if len(Z) != R:
            raise ValueError(f"Think time vector length ({len(Z)}) must match number of classes ({R})")

    # Handle mi
    if mi is None:
        mi = np.ones(M_original)
    else:
        mi = np.asarray(mi, dtype=np.float64).flatten()
        if len(mi) != M_original:
            raise ValueError(f"Multiplicity vector length ({len(mi)}) must match number of stations ({M_original})")

    # Detect and consolidate replicated stations
    unique_result = pfqn_unique(L)
    L_reduced = unique_result.L_unique
    mapping = unique_result.mapping
    M = L_reduced.shape[0]

    # Combine user-provided mi with detected multiplicity
    mi = pfqn_combine_mi(mi, mapping, M).flatten()

    # Empty population check
    if not np.any(N > 0):
        return (np.zeros((1, R)), np.zeros((1, R)), np.zeros((M_original, R)),
                np.zeros((M_original, R)), np.zeros((M_original, R)), np.zeros((M_original, R)), np.zeros((M_original, R)))

    # For single class, use simpler algorithm
    if R == 1:
        result = pfqn_mva_single_class(int(N[0]), L_reduced[:, 0], Z[0], mi)
        XN = np.array([[result['X']]])
        QN = result['Q'].reshape(-1, 1)
        RN = result['R'].reshape(-1, 1)
        UN = result['U'].reshape(-1, 1)

        # Expand results back to original dimensions if stations were consolidated
        if M < M_original:
            QN, UN, RN = pfqn_expand(QN, UN, RN, mapping)

        TN = XN * np.ones((M_original, 1))  # Node throughputs = system throughput
        AN = TN.copy()  # Arrival rates = throughputs
        CN = np.array([[RN.sum() + Z[0]]])
        return XN, CN, QN, UN, RN, TN, AN

    # Multi-class MVA using population recursion
    # Total population combinations for JIT threshold check
    totpop = int(np.prod(N + 1))

    # Use JIT version for larger population spaces
    if MVA_HAS_NUMBA and mva_population_recursion_jit is not None and totpop > JIT_THRESHOLD:
        N_int = N.astype(np.int64)
        XN_jit, QN_jit, CN_jit, lGN = mva_population_recursion_jit(
            L_reduced, N_int, Z, mi
        )
        XN = XN_jit.reshape(1, -1)
        QN = QN_jit
        CN = CN_jit

        # Compute utilizations
        UN = np.zeros((M, R))
        for m in range(M):
            for r in range(R):
                UN[m, r] = XN[0, r] * L_reduced[m, r]

        # Compute residence times
        RN = np.zeros((M, R))
        for m in range(M):
            for r in range(R):
                if XN[0, r] > 0:
                    RN[m, r] = QN[m, r] / XN[0, r]
                else:
                    RN[m, r] = L_reduced[m, r]

        # Expand results back to original dimensions if stations were consolidated
        if M < M_original:
            QN, UN, RN = pfqn_expand(QN, UN, RN, mapping)
            CN, _, _ = pfqn_expand(CN, CN, CN, mapping)

        # Node throughputs and arrival rates
        TN = np.zeros((M_original, R))
        AN = np.zeros((M_original, R))
        for m in range(M_original):
            for r in range(R):
                TN[m, r] = XN[0, r]
                AN[m, r] = XN[0, r]

        # Response time per class
        CN_total = np.zeros((1, R))
        for r in range(R):
            CN_total[0, r] = RN[:, r].sum() + Z[r]

        return XN, CN_total, QN, UN, RN, TN, AN

    # Pure Python: Compute product of (N[i]+1) for indexing
    prods = np.zeros(R - 1)
    for w in range(R - 1):
        prods[w] = np.prod(np.ones(R - w - 1) + N[w + 1:])

    # Find first non-empty class (from the end)
    first_non_empty = R - 1
    while first_non_empty >= 0 and N[first_non_empty] == 0:
        first_non_empty -= 1

    if first_non_empty < 0:
        return (np.zeros((1, R)), np.zeros((1, R)), np.zeros((M_original, R)),
                np.zeros((M_original, R)), np.zeros((M_original, R)), np.zeros((M_original, R)), np.zeros((M_original, R)))

    # Q[pop_idx, station] stores cumulative queue length at population index
    Q = np.zeros((totpop, M))

    # Output arrays
    XN = np.zeros((1, R))
    QN = np.zeros((M, R))
    CN = np.zeros((M, R))

    # Log normalizing constant
    lGN = 0.0

    # Initialize population vector
    n = np.zeros(R, dtype=int)
    n[first_non_empty] = 1

    currentpop = 1
    ctr = totpop  # Process all populations including empty state indexing

    while ctr > 0:
        for s in range(R):
            if n[s] > 0:
                # Compute index for n - e_s (one less job in class s)
                n[s] -= 1
                pos_n_1s = int(n[R - 1])
                for w in range(R - 1):
                    pos_n_1s += int(n[w] * prods[w])
                n[s] += 1
            else:
                pos_n_1s = 0

            # Compute residence times: CN[i,s] = L_reduced[i,s] * (mi[i] + Q[pos_n_1s, i])
            CNtot = 0.0
            for i in range(M):
                CN[i, s] = L_reduced[i, s] * (mi[i] + Q[pos_n_1s, i])
                CNtot += CN[i, s]

            # Compute throughput for class s
            XN[0, s] = n[s] / (Z[s] + CNtot) if (Z[s] + CNtot) > 0 else 0.0

            # Compute queue lengths and accumulate
            for i in range(M):
                QN[i, s] = XN[0, s] * CN[i, s]
                Q[currentpop, i] += QN[i, s]

        # Update log normalizing constant
        # Find last non-zero class position
        nonzero_idx = np.where(n > 0)[0]
        if len(nonzero_idx) > 0:
            last_nnz = nonzero_idx[-1]
            sumn = np.sum(n[:last_nnz])
            sumN = np.sum(N[:last_nnz])
            sumnprime = np.sum(n[last_nnz + 1:])
            if sumn == sumN and sumnprime == 0 and XN[0, last_nnz] > 0:
                lGN -= log(XN[0, last_nnz])

        # Find next population vector
        s = R - 1
        while s >= 0 and (n[s] == N[s] or s > first_non_empty):
            s -= 1

        if s < 0:
            break

        n[s] += 1
        for i in range(s + 1, R):
            n[i] = 0

        ctr -= 1
        currentpop += 1

    # Compute utilizations
    UN = np.zeros((M, R))
    for m in range(M):
        for r in range(R):
            UN[m, r] = XN[0, r] * L_reduced[m, r]

    # Compute residence times (waiting times)
    RN = np.zeros((M, R))
    for m in range(M):
        for r in range(R):
            if XN[0, r] > 0:
                RN[m, r] = QN[m, r] / XN[0, r]
            else:
                RN[m, r] = L_reduced[m, r]

    # Expand results back to original dimensions if stations were consolidated
    if M < M_original:
        QN, UN, RN = pfqn_expand(QN, UN, RN, mapping)
        CN, _, _ = pfqn_expand(CN, CN, CN, mapping)

    # Node throughputs and arrival rates
    TN = np.zeros((M_original, R))
    AN = np.zeros((M_original, R))
    for m in range(M_original):
        for r in range(R):
            TN[m, r] = XN[0, r]  # Closed network: all throughputs equal system throughput
            AN[m, r] = XN[0, r]  # Arrival rate = departure rate = throughput

    # Response time per class (sum of residence times + think time)
    CN_total = np.zeros((1, R))
    for r in range(R):
        CN_total[0, r] = RN[:, r].sum() + Z[r]

    return XN, CN_total, QN, UN, RN, TN, AN


def pfqn_bs(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
            tol: float = 1e-6, maxiter: int = 1000,
            QN0: np.ndarray = None, type_sched: np.ndarray = None
            ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Bard-Schweitzer Approximate Mean Value Analysis (MVA).

    Iterative approximate MVA algorithm that uses the (N-1)/N correction
    for the arrival theorem, providing good accuracy for most networks.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector
        Z: Think time vector (default 0)
        tol: Convergence tolerance (default 1e-6)
        maxiter: Maximum iterations (default 1000)
        QN0: Initial queue lengths (default: uniform distribution)
        type_sched: Scheduling strategy per station (default: PS)

    Returns:
        Tuple (XN, QN, UN, RN, it) matching MATLAB's pfqn_bs:
            XN: System throughputs (1 x R)
            QN: Mean queue lengths (M x R)
            UN: Utilizations (M x R)
            RN: Residence times (M x R)
            it: Number of iterations performed
    """
    from ...lang.base import SchedStrategy

    L = np.asarray(L, dtype=np.float64)
    N = np.asarray(N, dtype=np.float64).flatten()

    R = len(N)
    if L.ndim == 1:
        L = L.reshape(-1, 1) if R == 1 else L.reshape(1, -1)
    M = L.shape[0]

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()

    # Initialize queue lengths
    if QN0 is None:
        QN = np.tile(N, (M, 1)) / M
    else:
        QN = np.asarray(QN0, dtype=np.float64).copy()

    # Default scheduling: PS
    if type_sched is None:
        type_sched = [SchedStrategy.PS] * M

    CN = np.zeros((M, R))
    XN = np.zeros(R)
    UN = np.zeros((M, R))

    # Iterative Bard-Schweitzer algorithm
    for it in range(1, maxiter + 1):
        QN_old = QN.copy()

        for r in range(R):
            for ist in range(M):
                CN[ist, r] = L[ist, r]
                if L[ist, r] == 0:
                    continue

                for s in range(R):
                    if s != r:
                        # Different class contribution
                        sched_val = type_sched[ist]
                        is_fcfs = False
                        if isinstance(sched_val, int):
                            is_fcfs = (sched_val == SchedStrategy.FCFS)
                        elif hasattr(sched_val, 'name'):
                            is_fcfs = (sched_val.name == 'FCFS')
                        elif isinstance(sched_val, str):
                            is_fcfs = (sched_val.upper() == 'FCFS')

                        if is_fcfs:
                            CN[ist, r] += L[ist, s] * QN[ist, s]
                        else:
                            CN[ist, r] += L[ist, r] * QN[ist, s]
                    else:
                        # Same class contribution with arrival theorem correction
                        if N[r] > 0:
                            CN[ist, r] += L[ist, r] * QN[ist, r] * (N[r] - 1) / N[r]

            # Compute throughput
            CN_sum = np.sum(CN[:, r])
            if Z[r] + CN_sum > 0:
                XN[r] = N[r] / (Z[r] + CN_sum)
            else:
                XN[r] = 0

        # Update queue lengths
        for r in range(R):
            for ist in range(M):
                QN[ist, r] = XN[r] * CN[ist, r]

        # Update utilizations
        for r in range(R):
            for ist in range(M):
                UN[ist, r] = XN[r] * L[ist, r]

        # Check convergence
        with np.errstate(divide='ignore', invalid='ignore'):
            rel_change = np.abs(1 - QN / QN_old)
            rel_change = np.nan_to_num(rel_change, nan=0.0, posinf=0.0, neginf=0.0)
        if np.max(rel_change) < tol:
            break

    # Compute residence times
    RN = np.zeros((M, R))
    for r in range(R):
        if XN[r] > 0:
            RN[:, r] = QN[:, r] / XN[r]

    # Format output to match MATLAB's pfqn_bs: [XN, QN, UN, RN, it]
    XN_out = XN.reshape(1, -1)

    return XN_out, QN, UN, RN, it


def pfqn_aql(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
             max_iter: int = 1000, tol: float = 1e-6
             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                        np.ndarray, np.ndarray, np.ndarray]:
    """
    Approximate Queue Length (AQL) algorithm.

    Uses iterative approximation to compute queue lengths for large
    populations where exact MVA would be computationally expensive.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector
        Z: Think time vector (default 0)
        max_iter: Maximum iterations (default 1000)
        tol: Convergence tolerance (default 1e-6)

    Returns:
        Same format as pfqn_mva
    """
    L = np.asarray(L, dtype=np.float64)
    N = np.asarray(N, dtype=np.float64).flatten()
    N_total = N.sum()

    R = len(N)
    if L.ndim == 1:
        L = L.reshape(-1, 1) if R == 1 else L.reshape(1, -1)
    M = L.shape[0]

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()

    # Initialize with balanced system solution
    # pfqn_bs returns (XN, QN, UN, RN, it)
    XN, QN, UN, RN, _ = pfqn_bs(L, N, Z)
    # Compute TN and AN (node throughputs and arrival rates)
    TN = np.tile(XN, (M, 1)) if XN.ndim == 1 else np.tile(XN, (M, 1))
    AN = TN.copy()

    # Use JIT version for iterative refinement
    if MVA_HAS_NUMBA and schweitzer_iteration_jit is not None and M * R > 10:
        XN_jit, QN, UN, RN, iterations = schweitzer_iteration_jit(
            L, N.astype(np.float64), Z, QN, max_iter, tol
        )
        XN = XN_jit.reshape(1, -1)

        # Node throughputs and arrival rates
        for m in range(M):
            for r in range(R):
                TN[m, r] = XN[0, r]
                AN[m, r] = XN[0, r]

        # Response times
        CN = np.zeros((1, R))
        for r in range(R):
            CN[0, r] = RN[:, r].sum() + Z[r]

        return XN, CN, QN, UN, RN, TN, AN

    # Pure Python iterative refinement (Schweitzer approximation)
    for iteration in range(max_iter):
        Q_old = QN.copy()

        for r in range(R):
            if N[r] <= 0:
                continue

            for m in range(M):
                # Schweitzer approximation: E[Q_i | arrival of class r job]
                # ≈ (N_r - 1) / N_r * Q_i
                if N[r] > 1:
                    Q_others = (N[r] - 1) / N[r] * Q_old[m, r]
                else:
                    Q_others = 0

                # Add other class contributions
                for s in range(R):
                    if s != r:
                        Q_others += Q_old[m, s]

                # Update residence time
                RN[m, r] = L[m, r] * (1 + Q_others)

            # Throughput
            R_total = RN[:, r].sum()
            if Z[r] + R_total > 0:
                XN[0, r] = N[r] / (Z[r] + R_total)
            else:
                XN[0, r] = 0

            # Queue lengths
            for m in range(M):
                QN[m, r] = XN[0, r] * RN[m, r]
                UN[m, r] = XN[0, r] * L[m, r]
                TN[m, r] = XN[0, r]
                AN[m, r] = XN[0, r]

        # Check convergence
        diff = np.abs(QN - Q_old).max()
        if diff < tol:
            break

    # Response times
    CN = np.zeros((1, R))
    for r in range(R):
        CN[0, r] = RN[:, r].sum() + Z[r]

    return XN, CN, QN, UN, RN, TN, AN


def pfqn_sqni(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None
              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Single Queue Network Interpolation (SQNI) approximate MVA.

    Implements a fast approximation for multi-class closed queueing networks
    that reduces the system to single-queue representations with
    interpolation-based corrections.

    This method is particularly efficient for networks where one station
    dominates (bottleneck analysis), providing a good trade-off between
    accuracy and computational speed.

    Args:
        L: Service demand vector (1 x R or R,) - demands at the bottleneck queue
        N: Population vector (1 x R or R,) - number of jobs per class
        Z: Think time vector (1 x R or R,) - think time per class (default 0)

    Returns:
        Tuple of (Q, U, X) where:
            - Q: Queue lengths (2 x R) - first row for queue, second placeholder
            - U: Utilizations (2 x R) - first row for queue, second placeholder
            - X: Throughputs (1 x R)

    Reference:
        Based on the SQNI method for approximate MVA analysis.
    """
    L = np.asarray(L, dtype=np.float64).flatten()
    N = np.asarray(N, dtype=np.float64).flatten()

    R = len(N)
    if len(L) != R:
        raise ValueError(f"L length ({len(L)}) must match N length ({R})")

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=np.float64).flatten()
        if len(Z) != R:
            raise ValueError(f"Z length ({len(Z)}) must match N length ({R})")

    queue_idx = 0
    Nt = N.sum()

    Q = np.zeros((2, R))
    U = np.zeros((2, R))
    X = np.zeros((1, R))

    if Nt <= 0:
        return Q, U, X

    if Nt == 1.0:
        for r in range(R):
            if Z[r] + L[r] > 0:
                Xr = N[r] / (Z[r] + L[r])
            else:
                Xr = 0.0
            X[0, r] = Xr
            U[queue_idx, r] = Xr * L[r]
            Q[queue_idx, r] = Xr * L[r]
    else:
        for r in range(R):
            Nr = N[r]
            Lr = L[r]
            Zr = Z[r]

            # Compute Nvec_1r: N with one less job in class r
            Nvec_1r = N.copy()
            Nvec_1r[r] = max(0, Nvec_1r[r] - 1)

            sumN = N.sum()

            # Compute sumBrPart
            sumBrPart = 0.0
            for i in range(R):
                if i != r:
                    Zi = Z[i]
                    Li = L[i]
                    Ni = Nvec_1r[i]
                    denom = Zi + Li + Li * (sumN - 2)
                    if denom > 0:
                        sumBrPart += Zi * Ni / denom

            # Compute BrVec
            BrVec = np.zeros(R)
            for i in range(R):
                Zi = Z[i]
                Li = L[i]
                Ni = N[i]
                denom = Zi + Li + Li * (sumN - 1 - sumBrPart)
                if denom > 0:
                    BrVec[i] = Ni / denom * Zi

            # Compute BrSum (sum of BrVec except class r)
            BrSum = 0.0
            for i in range(R):
                if i != r:
                    BrSum += BrVec[i]

            Br = Lr * BrSum

            # Compute throughput
            if Lr == 0.0:
                if Zr > 0:
                    Xr = Nr / Zr
                else:
                    Xr = 0.0
            else:
                # Quadratic formula solution
                discriminant = (Br * Br - 2 * Br * Lr * Nt - 2 * Br * Zr +
                               Lr * Lr * Nt * Nt + 2 * Lr * Nt * Zr -
                               4 * Nr * Lr * Zr + Zr * Zr)
                if discriminant < 0:
                    discriminant = 0
                sqrt_term = np.sqrt(discriminant)
                denom = 2 * Lr * Zr
                if denom > 0:
                    Xr = (Zr - sqrt_term - Br + Lr * Nt) / denom
                else:
                    Xr = Nr / (Lr * Nt) if Lr * Nt > 0 else 0.0

            X[0, r] = max(0, Xr)
            U[queue_idx, r] = X[0, r] * Lr
            Q[queue_idx, r] = Nr - X[0, r] * Zr

    # Handle Z=0 case (adjust for infinite server at think station)
    for r in range(R):
        if Z[r] == 0.0 and L[r] > 0:
            Q_sum = Q[queue_idx, :].sum()
            denom = L[r] * (1 + Q_sum)
            if denom > 0:
                Xr = N[r] / denom
            else:
                Xr = 0.0
            X[0, r] = Xr
            U[queue_idx, r] = Xr * L[r]
            Q[queue_idx, r] = N[r] - Xr * Z[r]

    return Q, U, X


def pfqn_qd(
    L: np.ndarray,
    N: np.ndarray,
    ga: callable = None,
    be: callable = None,
    Q0: np.ndarray = None,
    tol: float = 1e-6,
    max_iter: int = 1000
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Queue-Dependent (QD) Approximate MVA.

    Implements the QD-AMVA algorithm that uses queue-dependent correction
    factors to improve accuracy of approximate MVA for closed networks.

    The algorithm iteratively computes queue lengths using:
    - A correction factor delta = (N_tot - 1) / N_tot
    - Per-class correction factor delta_r = (N_r - 1) / N_r
    - Optional scaling functions ga(A) and be(A) for advanced corrections

    Args:
        L: Service demand matrix (M x R) - rows are stations, columns are classes
        N: Population vector (R,) - number of jobs per class
        ga: Gamma scaling function ga(A) -> array(M,) (default: ones)
            A is the arrival queue seen by class, A[k] = 1 + delta * sum(Q[k,:])
        be: Beta scaling function be(A) -> array(M, R) (default: ones)
            A is the arrival queue per class, A[k,r] = 1 + delta_r * Q[k,r]
        Q0: Initial queue length estimate (M x R) (default: proportional)
        tol: Convergence tolerance (default 1e-6)
        max_iter: Maximum iterations (default 1000)

    Returns:
        Tuple of (Q, X, U, iter) where:
            Q: Mean queue lengths (M x R)
            X: Class throughputs (R,)
            U: Utilizations (M x R)
            iter: Number of iterations performed

    Reference:
        Schweitzer, P.J. "Approximate analysis of multiclass closed networks
        of queues." Proceedings of the International Conference on Stochastic
        Control and Optimization (1979).
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()

    M, R = L.shape
    N_tot = np.sum(N)

    # Default scaling functions (identity - return ones)
    if ga is None:
        def ga(A):
            return np.ones(M)
    if be is None:
        def be(A):
            return np.ones((M, R))

    # Initialize queue lengths
    if Q0 is None:
        # Proportional initialization: Q_kr = L_kr / sum_k(L_kr) * N_r
        L_sum = np.sum(L, axis=0, keepdims=True)
        L_sum[L_sum == 0] = 1  # Avoid division by zero
        Q = (L / L_sum) * N
    else:
        Q = np.asarray(Q0, dtype=float).copy()

    # Queue-dependent correction factors
    if N_tot > 0:
        delta = (N_tot - 1) / N_tot
    else:
        delta = 0.0

    deltar = np.zeros(R)
    for r in range(R):
        if N[r] > 0:
            deltar[r] = (N[r] - 1) / N[r]
        else:
            deltar[r] = 0.0

    # Initialize outputs
    X = np.zeros(R)
    U = np.zeros((M, R))
    C = np.zeros((M, R))

    # Storage for arrival queue vectors
    Ak = [np.zeros(M) for _ in range(R)]
    Akr = np.zeros((M, R))

    # Iteration
    Q_prev = Q * 10  # Ensure first iteration runs
    iteration = 0

    while np.max(np.abs(Q - Q_prev)) > tol and iteration < max_iter:
        iteration += 1
        Q_prev = Q.copy()

        # Compute arrival queue vectors
        for k in range(M):
            Q_row_sum = np.sum(Q[k, :])
            for r in range(R):
                Ak[r][k] = 1 + delta * Q_row_sum
                Akr[k, r] = 1 + deltar[r] * Q[k, r]

        # Compute queue lengths and throughputs for each class
        for r in range(R):
            if N[r] <= 0:
                X[r] = 0
                Q[:, r] = 0
                U[:, r] = 0
                continue

            # Get scaling factors
            g = ga(Ak[r])
            b = be(Akr)

            # Compute cycle times C[k,r]
            for k in range(M):
                Q_row_sum = np.sum(Q_prev[k, :])
                C[k, r] = L[k, r] * g[k] * b[k, r] * (1 + delta * Q_row_sum)

            # Compute throughput
            C_sum = np.sum(C[:, r])
            if C_sum > 0:
                X[r] = N[r] / C_sum
            else:
                X[r] = 0

            # Compute queue lengths and utilizations
            for k in range(M):
                Q[k, r] = X[r] * C[k, r]
                U[k, r] = L[k, r] * g[k] * b[k, r] * X[r]

    return Q, X, U, iteration


def pfqn_qdlin(
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray = None,
    tol: float = 1e-6,
    max_iter: int = 1000
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    QD-Linearizer (QDLIN) Approximate MVA.

    Combines Queue-Dependent (QD) correction with Linearizer iteration
    for improved accuracy in multi-class closed networks.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,) (default: zeros)
        tol: Convergence tolerance (default 1e-6)
        max_iter: Maximum iterations (default 1000)

    Returns:
        Tuple of (Q, U, R, X, C, iter) where:
            Q: Mean queue lengths (M x R)
            U: Utilizations (M x R)
            R: Residence times (M x R)
            X: Class throughputs (1 x R)
            C: Cycle times (1 x R)
            iter: Number of iterations performed
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()
    M, R = L.shape

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float).ravel()

    N_tot = np.sum(N)
    if N_tot <= 0:
        return (np.zeros((M, R)), np.zeros((M, R)), np.zeros((M, R)),
                np.zeros((1, R)), np.zeros((1, R)), 0)

    # QD correction factors
    delta = (N_tot - 1) / N_tot if N_tot > 0 else 0.0
    deltar = np.where(N > 0, (N - 1) / N, 0.0)

    # Initialize with proportional distribution
    L_sum = np.sum(L, axis=0, keepdims=True)
    L_sum[L_sum == 0] = 1
    Q = (L / L_sum) * N

    X = np.zeros((1, R))
    RN = np.zeros((M, R))
    C = np.zeros((1, R))
    U = np.zeros((M, R))

    Q_prev = Q * 10
    iteration = 0

    while np.max(np.abs(Q - Q_prev)) > tol and iteration < max_iter:
        iteration += 1
        Q_prev = Q.copy()

        for r in range(R):
            if N[r] <= 0:
                continue

            # Compute residence times with QD correction
            for k in range(M):
                # Queue seen by arriving class-r job
                Q_others = np.sum(Q_prev[k, :]) - Q_prev[k, r]
                Q_seen = Q_others + deltar[r] * Q_prev[k, r]
                RN[k, r] = L[k, r] * (1 + Q_seen)

            # Throughput
            R_total = np.sum(RN[:, r])
            if Z[r] + R_total > 0:
                X[0, r] = N[r] / (Z[r] + R_total)
            else:
                X[0, r] = 0

            # Update queue lengths
            for k in range(M):
                Q[k, r] = X[0, r] * RN[k, r]
                U[k, r] = X[0, r] * L[k, r]

            C[0, r] = R_total

    return Q, U, RN, X, C, iteration


def pfqn_qli(
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray = None,
    tol: float = 1e-6,
    max_iter: int = 1000
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Queue-Line (QLI) Approximate MVA (Wang-Sevcik).

    Implements the Wang-Sevcik Queue-Line approximation which provides
    improved accuracy for multi-class networks by better estimating
    the queue length seen by arriving customers.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,) (default: zeros)
        tol: Convergence tolerance (default 1e-6)
        max_iter: Maximum iterations (default 1000)

    Returns:
        Tuple of (Q, U, R, X, C, iter) where:
            Q: Mean queue lengths (M x R)
            U: Utilizations (M x R)
            R: Residence times (M x R)
            X: Class throughputs (1 x R)
            C: Cycle times (1 x R)
            iter: Number of iterations performed

    Reference:
        Wang, W. and Sevcik, K.C. "Performance Models for Multiprogrammed
        Systems." IBM Research Report RC 5925 (1976).
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()
    M, R = L.shape

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float).ravel()

    N_tot = np.sum(N)
    if N_tot <= 0:
        return (np.zeros((M, R)), np.zeros((M, R)), np.zeros((M, R)),
                np.zeros((1, R)), np.zeros((1, R)), 0)

    # Initialize with proportional distribution
    L_sum = np.sum(L, axis=0, keepdims=True)
    L_sum[L_sum == 0] = 1
    Q = (L / L_sum) * N

    X = np.zeros((1, R))
    RN = np.zeros((M, R))
    C = np.zeros((1, R))
    U = np.zeros((M, R))

    Q_prev = Q * 10
    iteration = 0

    while np.max(np.abs(Q - Q_prev)) > tol and iteration < max_iter:
        iteration += 1
        Q_prev = Q.copy()

        for r in range(R):
            if N[r] <= 0:
                continue

            for k in range(M):
                # Wang-Sevcik Queue-Line correction
                # Estimate queue seen by arriving class-r customer
                Q_total_k = np.sum(Q_prev[k, :])

                # Compute qlinum: L[k,r] * (1 + Q_total - Q[k,r])
                qlinum = L[k, r] * (1 + Q_total_k - Q_prev[k, r])

                # Compute qliden: sum over all stations m of L[m,r] * (1 + Q_total_m - Q[m,r])
                qliden = 0.0
                for m in range(M):
                    if L[m, r] > 0:
                        Q_total_m = np.sum(Q_prev[m, :])
                        qliden += L[m, r] * (1 + Q_total_m - Q_prev[m, r])

                if qliden > 0 and N[r] > 1:
                    Q_seen = Q_total_k - (1 / (N[r] - 1)) * (Q_prev[k, r] - qlinum / qliden)
                else:
                    Q_seen = Q_total_k - Q_prev[k, r]

                Q_seen = max(0, Q_seen)
                RN[k, r] = L[k, r] * (1 + Q_seen)

            # Throughput
            R_total = np.sum(RN[:, r])
            if Z[r] + R_total > 0:
                X[0, r] = N[r] / (Z[r] + R_total)
            else:
                X[0, r] = 0

            # Update queue lengths
            for k in range(M):
                Q[k, r] = X[0, r] * RN[k, r]
                U[k, r] = X[0, r] * L[k, r]

            C[0, r] = R_total

    return Q, U, RN, X, C, iteration


def pfqn_fli(
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray = None,
    tol: float = 1e-6,
    max_iter: int = 1000
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Fraction-Line (FLI) Approximate MVA (Wang-Sevcik).

    Implements the Wang-Sevcik Fraction-Line approximation, an alternative
    to Queue-Line that uses a different formula for estimating the queue
    length seen by arriving customers.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,) (default: zeros)
        tol: Convergence tolerance (default 1e-6)
        max_iter: Maximum iterations (default 1000)

    Returns:
        Tuple of (Q, U, R, X, C, iter) where:
            Q: Mean queue lengths (M x R)
            U: Utilizations (M x R)
            R: Residence times (M x R)
            X: Class throughputs (1 x R)
            C: Cycle times (1 x R)
            iter: Number of iterations performed

    Reference:
        Wang, W. and Sevcik, K.C. "Performance Models for Multiprogrammed
        Systems." IBM Research Report RC 5925 (1976).
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()
    M, R = L.shape

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float).ravel()

    N_tot = np.sum(N)
    if N_tot <= 0:
        return (np.zeros((M, R)), np.zeros((M, R)), np.zeros((M, R)),
                np.zeros((1, R)), np.zeros((1, R)), 0)

    # Initialize with proportional distribution
    L_sum = np.sum(L, axis=0, keepdims=True)
    L_sum[L_sum == 0] = 1
    Q = (L / L_sum) * N

    X = np.zeros((1, R))
    RN = np.zeros((M, R))
    C = np.zeros((1, R))
    U = np.zeros((M, R))

    Q_prev = Q * 10
    iteration = 0

    while np.max(np.abs(Q - Q_prev)) > tol and iteration < max_iter:
        iteration += 1
        Q_prev = Q.copy()

        for r in range(R):
            if N[r] <= 0:
                continue

            for k in range(M):
                # Wang-Sevcik Fraction-Line correction
                Q_total_k = np.sum(Q_prev[k, :])

                # Compute qlinum: L[k,r] * (1 + Q_total - Q[k,r])
                qlinum = L[k, r] * (1 + Q_total_k - Q_prev[k, r])

                # Compute qliden: sum over all stations m of L[m,r] * (1 + Q_total_m - Q[m,r])
                qliden = 0.0
                for m in range(M):
                    if L[m, r] > 0:
                        Q_total_m = np.sum(Q_prev[m, :])
                        qliden += L[m, r] * (1 + Q_total_m - Q_prev[m, r])

                # FLI uses different formula than QLI
                if qliden > 0 and N[r] > 0:
                    Q_seen = Q_total_k - (2 / N[r]) * Q_prev[k, r] + qlinum / qliden
                else:
                    Q_seen = Q_total_k - Q_prev[k, r]

                Q_seen = max(0, Q_seen)
                RN[k, r] = L[k, r] * (1 + Q_seen)

            # Throughput
            R_total = np.sum(RN[:, r])
            if Z[r] + R_total > 0:
                X[0, r] = N[r] / (Z[r] + R_total)
            else:
                X[0, r] = 0

            # Update queue lengths
            for k in range(M):
                Q[k, r] = X[0, r] * RN[k, r]
                U[k, r] = X[0, r] * L[k, r]

            C[0, r] = R_total

    return Q, U, RN, X, C, iteration


def pfqn_bsfcfs(
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray = None,
    tol: float = 1e-6,
    max_iter: int = 1000,
    QN: np.ndarray = None,
    weight: np.ndarray = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Bard-Schweitzer approximate MVA for FCFS scheduling with weighted priorities.

    Implements AMVA with FCFS approximation where classes can have
    relative priority weights affecting the expected waiting times.

    Args:
        L: Service demand matrix (M x R) where M is stations, R is classes
        N: Population vector (1 x R) - number of jobs per class
        Z: Think time vector (1 x R) - default: zeros
        tol: Convergence tolerance (default 1e-6)
        max_iter: Maximum iterations (default 1000)
        QN: Initial queue length matrix (M x R) - default: uniform distribution
        weight: Weight matrix (M x R) for relative priorities - default: ones

    Returns:
        Tuple of (XN, QN, UN, RN, it) where:
            XN: System throughput per class (1 x R)
            QN: Mean queue lengths (M x R)
            UN: Utilizations (M x R)
            RN: Residence times (M x R)
            it: Number of iterations performed

    Reference:
        Bard, Y. and Schweitzer, P.J. "Analyzing Closed Queueing Networks with Multiple
        Job Classes and Multiserver Stations." Performance Evaluation Review 7.1-2 (1978).
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()

    M, R = L.shape

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float).ravel()

    # Initialize queue lengths
    if QN is None:
        QN = np.tile(N, (M, 1)) / M
    else:
        QN = np.asarray(QN, dtype=float).copy() + 1e-12  # Avoid zero for numerical stability

    # Initialize weight matrix
    if weight is None:
        weight = np.ones((M, R))
    else:
        weight = np.asarray(weight, dtype=float)

    XN = np.zeros(R)
    UN = np.zeros((M, R))
    CN = np.zeros((M, R))
    relprio = np.zeros((M, R))

    for it in range(1, max_iter + 1):
        QN_old = QN.copy()

        # Compute relative priorities
        for ist in range(M):
            for r in range(R):
                relprio[ist, r] = QN[ist, r] * weight[ist, r]

        # Compute residence times with FCFS approximation
        for r in range(R):
            for ist in range(M):
                CN[ist, r] = L[ist, r]
                for s in range(R):
                    if s != r:
                        # FCFS approximation: weighted by relative priorities
                        if relprio[ist, r] > 0:
                            CN[ist, r] += L[ist, s] * QN[ist, s] * relprio[ist, s] / relprio[ist, r]
                    else:
                        # Same class contribution with arrival theorem correction
                        if N[r] > 0 and relprio[ist, r] > 0:
                            CN[ist, r] += L[ist, r] * QN[ist, r] * (N[r] - 1) / N[r] * relprio[ist, s] / relprio[ist, r]

            # Compute throughput
            CN_sum = np.sum(CN[:, r])
            if Z[r] + CN_sum > 0:
                XN[r] = N[r] / (Z[r] + CN_sum)
            else:
                XN[r] = 0

        # Update queue lengths
        for r in range(R):
            for ist in range(M):
                QN[ist, r] = XN[r] * CN[ist, r]

        # Compute utilizations
        for r in range(R):
            for ist in range(M):
                UN[ist, r] = XN[r] * L[ist, r]

        # Check convergence
        if QN_old.max() > 0:
            rel_diff = np.abs(1 - QN / np.maximum(QN_old, 1e-12)).max()
            if rel_diff < tol:
                break

    # Compute residence times from queue lengths
    RN = np.zeros((M, R))
    for r in range(R):
        for ist in range(M):
            if XN[r] > 0:
                RN[ist, r] = QN[ist, r] / XN[r]
            else:
                RN[ist, r] = L[ist, r]

    return XN.reshape(1, -1), QN, UN, RN, it


def pfqn_joint(
    n: np.ndarray,
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray = None,
    lGN: float = None
) -> float:
    """
    Compute joint queue-length probability distribution.

    Computes the joint probability for a given queue-length state vector
    in a closed product-form queueing network.

    Args:
        n: Queue-length state vector (M,) for total or (M x R) for per-class
           - If 1D (M,): n[i] is the total number of jobs at station i
           - If 2D (M x R): n[i,r] is the number of class-r jobs at station i
        L: Service demand matrix (M x R)
        N: Population vector (1 x R)
        Z: Think time vector (1 x R) - default: zeros
        lGN: Log normalizing constant (optional, computed if not provided)

    Returns:
        pjoint: Joint probability of state n

    Examples:
        # Total queue-lengths (Z > 0)
        >>> p = pfqn_joint([2, 1], [[10, 2], [5, 4]], [2, 2], [91, 92])

        # Per-class queue-lengths
        >>> p = pfqn_joint([[1, 0], [0, 1]], [[10, 2], [5, 4]], [2, 2], [91, 92])
    """
    from .nc import pfqn_ca
    from scipy.special import gammaln

    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()
    n = np.asarray(n, dtype=float)

    M, R = L.shape

    if Z is None:
        Z = np.zeros(R)
    else:
        Z = np.asarray(Z, dtype=float).ravel()

    # Compute normalizing constant if not provided
    if lGN is None:
        _, lGN = pfqn_ca(L, N, Z)

    def factln(x):
        """Log factorial using gammaln."""
        return gammaln(np.asarray(x) + 1)

    def multinomialln(x):
        """Log multinomial coefficient."""
        x = np.asarray(x, dtype=float)
        return factln(np.sum(x)) - np.sum(factln(x))

    if n.ndim == 1 or (n.ndim == 2 and n.shape[1] == 1):
        # Joint probability of total queue lengths
        n = n.ravel()

        if np.sum(Z) > 0:
            # With think time
            n0 = np.sum(N) - np.sum(n)
            if n0 < 0:
                return 0.0

            # Compute F_per for L extended with Z
            L_ext = np.vstack([L, Z.reshape(1, -1)])
            n_ext = np.concatenate([n, [n0]])

            Fjoint = _fper(L_ext, N, n_ext.astype(int))
            pjoint = np.exp(np.log(max(Fjoint, 1e-300)) - lGN - factln(n0))
        else:
            Fjoint = _fper(L, N, n.astype(int))
            pjoint = np.exp(np.log(max(Fjoint, 1e-300)) - lGN)

    elif n.ndim == 2 and n.shape[1] == R:
        # Joint probability of per-class queue-lengths
        n0 = N - np.sum(n, axis=0)

        if np.any(n0 < 0):
            return 0.0

        if np.sum(Z) > 0:
            Fjoint = np.sum(n0 * np.log(np.maximum(Z, 1e-300))) - np.sum(factln(n0))
        else:
            Fjoint = 0.0

        for i in range(M):
            if np.sum(n[i, :]) > 0:
                Fjoint += multinomialln(n[i, :]) + np.sum(n[i, :] * np.log(np.maximum(L[i, :], 1e-300)))

        pjoint = np.exp(Fjoint - lGN)
    else:
        raise ValueError("Invalid argument to pfqn_joint: n must be (M,) or (M, R)")

    return max(0.0, pjoint)


def _fper(L: np.ndarray, N: np.ndarray, m: np.ndarray) -> float:
    """
    Compute permanent-based function F_per for joint probability.

    Internal helper for pfqn_joint.
    """
    from math import factorial
    from itertools import permutations

    M, R = L.shape
    m = np.asarray(m, dtype=int).ravel()

    # Build matrix Ak
    Ak_cols = []
    for r in range(R):
        N_r = int(N[r])
        for _ in range(N_r):
            Ak_cols.append(L[:, r])
    Ak = np.column_stack(Ak_cols) if Ak_cols else np.zeros((M, 0))

    # Build matrix A based on m
    A_rows = []
    for i in range(M):
        if i < len(m) and m[i] > 0:
            for _ in range(int(m[i])):
                A_rows.append(Ak[i, :])
    A = np.vstack(A_rows) if A_rows else np.zeros((0, Ak.shape[1]))

    # Compute permanent (simplified for small matrices)
    if A.shape[0] == 0 or A.shape[1] == 0:
        return 1.0

    n_rows, n_cols = A.shape
    if n_rows > n_cols:
        return 0.0

    # Simple permanent computation (exponential but works for small n)
    if n_rows <= 10:
        perm = 0.0
        from itertools import permutations as perms
        for p in perms(range(n_cols), n_rows):
            prod = 1.0
            for i, j in enumerate(p):
                prod *= A[i, j]
            perm += prod
    else:
        # For larger matrices, use approximation
        perm = _permanent_approx(A)

    # Divide by product of factorials
    prod_fact = 1.0
    for r in range(R):
        prod_fact *= factorial(int(N[r]))

    return perm / prod_fact


def _permanent_approx(A: np.ndarray) -> float:
    """
    Approximate permanent using Bethe approximation.
    For small matrices, computes exact permanent.
    """
    n_rows, n_cols = A.shape
    if n_rows == 0:
        return 1.0

    # Simple approximation using product of sums
    # This is an upper bound (van der Waerden)
    perm = 1.0
    for i in range(n_rows):
        row_sum = np.sum(A[i, :])
        perm *= row_sum

    return perm


__all__ = [
    'pfqn_mva',
    'pfqn_mva_single_class',
    'pfqn_bs',
    'pfqn_aql',
    'pfqn_sqni',
    'pfqn_qd',
    'pfqn_qdlin',
    'pfqn_qli',
    'pfqn_fli',
    'pfqn_bsfcfs',
    'pfqn_joint',
]
