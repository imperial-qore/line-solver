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
    with load-independent stations. Standard arrival theorem; for the interlocked-flow
    correction of Franks (1999), Ch. 4, Eq. (4.7) call pfqn_mva_ilock instead.

    Args:
        L: Service demand matrix (M x R) where M is stations, R is classes
        N: Population vector (1 x R or R,) - number of jobs per class
        Z: Think time vector (1 x R or R,) - think time per class (default 0)
        mi: Additive term of the residence-time recursion
            C(i,s)=L(i,s)*(mi(i)+Qarv), 1 for a queueing station (default 1).
            THIS IS NOT A SERVER COUNT: mi(i)=c inflates the residence time by c
            rather than adding c servers. For multiserver stations call
            pfqn_mvams(lambda, L, N, Z, mi, S), which passes S to the
            load-dependent recursion with mu(i,n)=min(n,S(i)).

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

    # see _kb/03-api-layer.md for rationale
    L_reduced = L
    mapping = np.arange(M_original)
    M = M_original

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
    totpop = int(np.prod(N + 1))

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
                qarv = Q[pos_n_1s, i]
                CN[i, s] = L_reduced[i, s] * (mi[i] + qarv)
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
                # see _kb/03-api-layer.md for rationale
                RN[m, r] = CN[m, r]

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



def pfqn_mva_ilock(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
                   mi: np.ndarray = None, IL: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                              np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Exact MVA recursion carrying the interlocked-flow correction.

    The correction of Franks (1999), Ch. 4, Eq. (4.7) replaces the arrival theorem
    term Q(n-1_s,i) by a per-class weighted sum, so the recursion has to carry
    per-class queue lengths that pfqn_mva does not need. Closed single-server
    models only.

    The discounted arrival-instant queue is floored at the in-service component,
    as in lqns MVA::queueOnly_adjusted, so the correction damps itself out as a
    station saturates. That is a self-limiting guard, NOT a hard capacity test:
    sum_s XN[s]*L[i,s] <= mi[i] is still asserted nowhere.
    See git show 449847e7b:_kb/log.md.

    Args:
        L: Service demand matrix (M x R) where M is stations, R is classes
        N: Population vector (1 x R or R,) - number of jobs per class
        Z: Think time vector (1 x R or R,) - think time per class (default 0)
        mi: Additive term of the residence-time recursion
            C(i,s)=L(i,s)*(mi(i)+Qarv), 1 for a queueing station (default 1).
            THIS IS NOT A SERVER COUNT: mi(i)=c inflates the residence time by c
            rather than adding c servers. For multiserver stations call
            pfqn_mvams(lambda, L, N, Z, mi, S), which passes S to the
            load-dependent recursion with mu(i,n)=min(n,S(i)).
        IL: Interlock matrix (R x R), IL[r,s] is the share of the class-s queue that a
            class-r arrival cannot see, because that work was itself caused by the class-r
            request (Franks 1999, Eq. 4.7). Required; pass None to pfqn_mva instead.

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

    # see _kb/03-api-layer.md for rationale
    L_reduced = L
    mapping = np.arange(M_original)
    M = M_original

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

    if IL is None or np.size(IL) == 0:
        raise ValueError("an interlock matrix is required; use pfqn_mva for the standard arrival theorem")
    IL = np.asarray(IL, dtype=np.float64)
    if IL.shape != (R, R):
        raise ValueError(f"the interlock matrix must be {R}x{R}, got {IL.shape}")
    ILw = np.clip(1.0 - IL, 0.0, 1.0)
    np.fill_diagonal(ILw, 1.0)  # a request always sees its own class in full

    # Multi-class MVA using population recursion
    totpop = int(np.prod(N + 1))

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
    # per-class queue lengths, needed by the interlock
    Qc = np.zeros((totpop, M, R))
    # per-class in-service component, the interlock's floor
    Uc = np.zeros((totpop, M, R))

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
                # In-service protection, as in lqns MVA::queueOnly_adjusted: the
                # discount bites on the WAITING part only, never on the job already
                # in service, so it damps itself out as the station saturates.
                qarv = float(np.sum(np.maximum(ILw[s, :] * Qc[pos_n_1s, i, :],
                                               Uc[pos_n_1s, i, :])))
                CN[i, s] = L_reduced[i, s] * (mi[i] + qarv)
                CNtot += CN[i, s]

            # Compute throughput for class s
            XN[0, s] = n[s] / (Z[s] + CNtot) if (Z[s] + CNtot) > 0 else 0.0

            # Compute queue lengths and accumulate
            for i in range(M):
                QN[i, s] = XN[0, s] * CN[i, s]
                Q[currentpop, i] += QN[i, s]
                Qc[currentpop, i, s] = QN[i, s]
                Uc[currentpop, i, s] = XN[0, s] * L_reduced[i, s]

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

    lGN = np.nan  # the interlock leaves the model outside product form

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
                # see _kb/03-api-layer.md for rationale
                RN[m, r] = CN[m, r]

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
        tol: Convergence tolerance (default 1e-6); 'cn' or NaN selects the
            published Linearizer termination test of Chandy and Neuse,
            Commun. ACM 25(2), 1982, i.e. the cutoff pfqn_cntol(N) applied to
            max_{i,r}|dQ(i,r)|/N_r instead of the relative-change metric used
            by default. This is the test LQNS runs, since it sets it in
            SchweitzerCommon.
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
    from .cntol import is_cntol, pfqn_cntol

    L = np.asarray(L, dtype=np.float64)
    N = np.asarray(N, dtype=np.float64).flatten()

    cntest = is_cntol(tol)
    if cntest:
        tol = pfqn_cntol(N)

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
        if cntest:
            # Chandy and Neuse (1982), p.129: absolute queue-length change
            # scaled by the class population, over the non-empty classes only.
            nz = N > 0
            if not np.any(nz):
                break
            change = np.max(np.abs(QN[:, nz] - QN_old[:, nz]) / N[nz])
        else:
            with np.errstate(divide='ignore', invalid='ignore'):
                rel_change = np.abs(1 - QN / QN_old)
                rel_change = np.nan_to_num(rel_change, nan=0.0, posinf=0.0, neginf=0.0)
            change = np.max(rel_change)
        if change < tol:
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
             tol: float = 1e-7, max_iter: int = 1000, QN0: np.ndarray = None
             ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                        np.ndarray, np.ndarray, np.ndarray]:
    """
    Aggregate Queue Length (AQL) approximate MVA for closed product-form networks.

    Port of matlab/src/api/pfqn/pfqn_aql.m, cross-checked against
    jar/src/main/java/jline/api/pfqn/mva/Pfqn_aql.java. The fixed point carries
    K+1 population points (the full population and each N - e_s) and a
    correction gamma(k,s) = Q_0(k)/sum(N) - Q_s(k)/(sum(N)-1) that removes the
    Schweitzer proportionality error, in the manner of Linearizer.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (default 0)
        tol: Relative tolerance on the full-population queue lengths (default 1e-7)
        max_iter: Maximum iterations (default 1000)
        QN0: Warm start for the queue lengths (M x R), optional

    Returns:
        Tuple of (XN, CN, QN, UN, RN, TN, AN); AN holds the arrival-instant
        queue lengths Q_s(k), as in the MATLAB reference.
    """
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

    if QN0 is None or np.size(QN0) == 0:
        Q0 = np.tile(N, (M, 1)) / M
    else:
        Q0 = np.asarray(QN0, dtype=np.float64).reshape(M, R) + np.finfo(float).eps

    # Q[t], R[t], X[t] hold the solution at population N (t = 0) and at
    # N - e_{t-1} (t = 1..R). Q is aggregate (per station), as in the reference.
    # Q{t+1}(k,1)=QN0(k) in the reference: MATLAB linear indexing on an (M,R)
    # array reads the FIRST column, so every population point starts there.
    Qt = [Q0[:, 0].copy() for _ in range(R + 1)]
    Rt = [np.zeros((M, R)) for _ in range(R + 1)]
    Xt = [np.zeros(R) for _ in range(R + 1)]
    gamma = np.zeros((M, R))

    it = 0
    while True:
        Q_olditer = Qt[0].copy()
        it += 1
        for t in range(R + 1):
            n = N.copy()
            if t > 0:
                n[t - 1] = max(n[t - 1] - 1.0, 0.0)
            ntot = n.sum()
            for k in range(M):
                for s in range(R):
                    Rt[t][k, s] = L[k, s] * (1.0 + (ntot - 1.0) * (
                        (Qt[t][k] / ntot if ntot > 0 else 0.0) - gamma[k, s]))
            for s in range(R):
                den = Z[s] + Rt[t][:, s].sum()
                Xt[t][s] = n[s] / den if den > 0 else 0.0
            for k in range(M):
                Qt[t][k] = float(Xt[t] @ Rt[t][k, :])
        Ntot = N.sum()
        for k in range(M):
            for s in range(R):
                gamma[k, s] = (Qt[0][k] / Ntot if Ntot > 0 else 0.0) - (
                    Qt[s + 1][k] / (Ntot - 1.0) if Ntot > 1 else 0.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            rel = np.abs((Q_olditer - Qt[0]) / np.where(Qt[0] != 0, Qt[0], np.inf))
        if np.nanmax(rel) < tol or it == max_iter:
            break

    XN = Xt[0].reshape(1, -1)
    RN = Rt[0]
    UN = np.zeros((M, R))
    QN = np.zeros((M, R))
    AN = np.zeros((M, R))
    for k in range(M):
        for s in range(R):
            UN[k, s] = XN[0, s] * L[k, s]
            QN[k, s] = UN[k, s] * (1.0 + Qt[s + 1][k])
            AN[k, s] = Qt[s + 1][k]
    TN = np.tile(XN, (M, 1))
    CN = (RN.sum(axis=0) + Z).reshape(1, -1)

    return XN, CN, QN, UN, RN, TN, AN


def pfqn_sqni(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None
              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Square-root Non-iterative (SQNI) approximate MVA.

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
        # A Z=0 class (self-looping) has no delay to interpolate through: its
        # queue length is its whole population and it is solved after the loop.
        for r in range(R):
            if Z[r] == 0.0:
                Q[queue_idx, r] = N[r]

        for r in range(R):
            if Z[r] == 0.0:
                continue
            Nr = N[r]
            Lr = L[r]
            Zr = Z[r]

            # Compute Nvec_1r: N with one less job in class r
            Nvec_1r = N.copy()
            Nvec_1r[r] = max(0, Nvec_1r[r] - 1)

            sumN = N.sum()

            # sumBrPart runs over EVERY class, class r included, as in MATLAB
            # pfqn_sqni; skipping r shifted X by 1.6% on a 2-class model.
            sumBrPart = 0.0
            for i in range(R):
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
        # Joint probability of total queue lengths. The permanent identity owns
        # this branch (pfqn_jointmarg); the aggregated think time is one extra
        # infinite-server row, which is exact by the multinomial theorem.
        n = n.ravel()

        if np.sum(Z) > 0:
            n0 = np.sum(N) - np.sum(n)
            if n0 < 0:
                return 0.0
            L_ext = np.vstack([L, Z.reshape(1, -1)])
            n_ext = np.concatenate([n, [n0]])
            pjoint, _ = pfqn_jointmarg(n_ext, L_ext, N, [M], lGN)
        else:
            pjoint, _ = pfqn_jointmarg(n, L, N, [], lGN)

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


def pfqn_jointmarg(
    n: np.ndarray,
    L: np.ndarray,
    N: np.ndarray,
    infset=None,
    lGN: float = None,
    engine: str = 'exact'
):
    """
    Joint probability of the per-station TOTAL queue lengths.

    Joint probability that station i holds n[i] jobs IN TOTAL, all classes
    summed out, in a closed multiclass product-form network::

        P(n_1,...,n_M) = perm(A) / ( prod_r N_r! * prod_{j in infset} n_j! * G(N) )

    with A the demand matrix whose column r is repeated N[r] times and whose
    row i is repeated n[i] times, so A is square of order sum(N). Unlike
    pfqn_joint, which takes the delay as a single aggregated row, every
    infinite-server station keeps its own row here and contributes its own
    1/n_j!: the queueing stations contribute the n_i! that the permanent
    identity supplies, the infinite servers do not.

    Args:
        n: (M,) per-station total queue lengths, infinite servers included;
           sum(n) must equal sum(N)
        L: (M, R) demand matrix, infinite-server rows included
        N: (R,) per-class populations
        infset: row indices of L that are infinite-server stations, empty by
            default (every station is a queue)
        lGN: log normalizing constant; computed with pfqn_ca when omitted,
            aggregating the infinite-server rows into the think time (which is
            exact: the delay stations aggregate by the multinomial theorem, so
            G does not depend on how they are split)
        engine: 'exact' (default), 'spm', 'bethe', 'heur', 'huberlaw' or
            'adapart'. 'spm' is the only engine that does not expand the
            matrix to order sum(N): it takes the row-replicated matrix with
            the class populations as column multiplicities, which is the
            regime its saddle-point expansion is asymptotically exact in, so
            its cost does not grow with the population and its relative error
            is O((R-1)/min(N)). Measured on a 3-station 2-class model, 12.8%
            at N = (1,1), 4.2% at (3,3), 2.1% at (6,6); it degrades the other
            way round, when the class count grows at fixed population (2.7% at
            R = 2, 21% at R = 7, both at N_r = 3), because R-1 is the
            dimension being expanded in. The bias is nearly constant across
            the lattice, so a caller that renormalizes a full sweep keeps far
            less of it: total variation distance 5.0e-3 at N = (1,1), 8.4e-4
            at (3,3), 4.3e-4 at (5,5), better than 'bethe' and 'heur' at every
            population measured

    Returns:
        (pjoint, lpjoint): the joint probability and its logarithm, which
        survives populations the probability itself underflows at

    The identity holds for load-independent single-server queues plus infinite
    servers. Multiserver and load-dependent stations break the n_i! factor and
    are the caller's responsibility to exclude.

    ZERO ELEMENTS are safe under the exact engine and only under it: a station
    holding no jobs contributes no row, a class with no jobs contributes no
    column, a zero demand is an ordinary zero entry of A, and the permanent of
    the empty matrix is 1. The approximate engines are REFUSED on a matrix with
    a structural zero rather than having it floored at eps: Sinkhorn scaling
    needs full support, and the Bethe gap is a state-dependent lower bound that
    does not cancel when the estimates are normalized against each other.

    References:
        H. J. Ryser, "Combinatorial Mathematics", Carus Mathematical
        Monographs 14, Mathematical Association of America, 1963.
    """
    from .nc import pfqn_ca, pfqn_perm
    from scipy.special import gammaln

    L = np.atleast_2d(np.asarray(L, dtype=float))
    n = np.asarray(n, dtype=float).ravel()
    N = np.asarray(N, dtype=float).ravel()
    M, R = L.shape

    if infset is None:
        infset = []
    infset = np.asarray(infset, dtype=int).ravel()
    if engine is None or engine == '':
        engine = 'exact'
    engine = str(engine).lower()

    if n.size != M:
        raise ValueError("pfqn_jointmarg: the occupancy vector has %d entries but L has %d rows."
                         % (n.size, M))
    if N.size != R:
        raise ValueError("pfqn_jointmarg: the population vector has %d entries but L has %d columns."
                         % (N.size, R))
    if np.any(n < 0):
        raise ValueError("pfqn_jointmarg: the occupancy vector has a negative entry.")
    if infset.size > 0 and (np.any(infset < 0) or np.any(infset >= M)):
        raise ValueError("pfqn_jointmarg: infset indexes a station outside 0..%d." % (M - 1))

    # Infeasible occupancies are not an error: the caller sweeps a lattice.
    if int(round(np.sum(n))) != int(round(np.sum(N))):
        return 0.0, -np.inf

    if lGN is None or not np.isfinite(lGN):
        isinfrow = np.zeros(M, dtype=bool)
        isinfrow[infset] = True
        Lq = L[~isinfrow, :]
        Z = np.sum(L[isinfrow, :], axis=0) if np.any(isinfrow) else np.zeros(R)
        _, lGN = pfqn_ca(Lq, N.reshape(1, -1), Z.reshape(1, -1))

    if np.sum(N) == 0:
        lpjoint = -lGN
        return float(np.exp(lpjoint)), float(lpjoint)

    # The expanded matrix is square of order sum(N). 'spm' works on the
    # unexpanded form, and building this would throw away the very property
    # that makes it independent of the population.
    A = None if engine == 'spm' else _replicate_demands(L, N, n)

    if engine != 'exact':
        zero = _first_zero_demand(L, N, n)
        if zero is not None:
            raise ValueError(
                "pfqn_jointmarg: the '%s' permanent engine cannot be applied: the demand of class %d "
                "at station %d is zero, so the replicated matrix has no full support. Use engine 'exact'."
                % (engine, zero[1] + 1, zero[0] + 1))

    if engine == 'exact':
        F = pfqn_perm(A)
    elif engine == 'spm':
        # Never the expanded A: the saddle point is asymptotic in the column
        # multiplicities, which are the class populations themselves.
        from ..perm import perm_spm
        Ar, mr = _replicate_rows(L, N, n)
        F = perm_spm(Ar, mr)
    elif engine == 'bethe':
        from ..perm import perm_bethe
        F = perm_bethe(A)
    elif engine == 'heur':
        from ..perm import perm_heur
        F = perm_heur(A)
    elif engine == 'huberlaw':
        from ..perm import HuberLawSampler
        F = HuberLawSampler(A, solve=True).value
    elif engine == 'adapart':
        from ..perm import AdaPartSampler
        F = AdaPartSampler(A, solve=True).value
    else:
        raise ValueError("pfqn_jointmarg: unrecognized permanent engine '%s'. "
                         "Use exact, spm, bethe, heur, huberlaw or adapart." % engine)

    if F <= 0:
        return 0.0, -np.inf

    lpjoint = float(np.log(F) - np.sum(gammaln(N + 1)) - np.sum(gammaln(n[infset] + 1)) - lGN)
    return float(np.exp(lpjoint)), lpjoint


def _replicate_rows(L: np.ndarray, N: np.ndarray, n: np.ndarray):
    """
    Row i of L repeated n[i] times, with the class populations as multiplicities.

    The same matrix _replicate_demands expands, one step earlier: perm(Ar, m)
    equals perm(A), and perm_spm wants the unexpanded form because its
    expansion is asymptotic in m. A class with no jobs is dropped rather than
    passed with multiplicity zero, so a zero demand in such a column cannot
    trip the full-support check.
    """
    keepc = np.flatnonzero(np.asarray(N).ravel() > 0)
    m = np.asarray(N).ravel()[keepc]
    rows = []
    for i in range(L.shape[0]):
        rows.extend([L[i, keepc]] * int(round(n[i])))
    Ar = np.vstack(rows) if rows else np.zeros((0, keepc.size))
    return Ar, m


def _replicate_demands(L: np.ndarray, N: np.ndarray, n: np.ndarray) -> np.ndarray:
    """
    Column r of L repeated N[r] times, then row i of that repeated n[i] times.

    A station holding no jobs and a class holding no jobs each drop out here,
    which is what makes a zero entry of the occupancy vector free of any
    special case: the result stays square of order sum(N).
    """
    M = L.shape[0]
    cols = []
    for r in range(L.shape[1]):
        cols.extend([L[:, r]] * int(round(N[r])))
    Ak = np.column_stack(cols) if cols else np.zeros((M, 0))

    rows = []
    for i in range(M):
        rows.extend([Ak[i, :]] * int(round(n[i])))
    return np.vstack(rows) if rows else np.zeros((0, Ak.shape[1]))


def _first_zero_demand(L: np.ndarray, N: np.ndarray, n: np.ndarray):
    """
    First (station, class) whose zero demand actually reaches the replicated
    matrix. A class with no jobs or a station with no jobs contributes nothing,
    so its zeros are irrelevant.
    """
    for i in range(L.shape[0]):
        if n[i] == 0:
            continue
        for r in range(L.shape[1]):
            if N[r] > 0 and L[i, r] <= 0:
                return i, r
    return None


__all__ = [
    'pfqn_mva',
    'pfqn_mva_single_class',
    'pfqn_bs',
    'pfqn_aql',
    'pfqn_sqni',
    'pfqn_qli',
    'pfqn_fli',
    'pfqn_joint',
    'pfqn_jointmarg',
]
