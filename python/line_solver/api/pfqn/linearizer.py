"""
Linearizer Approximate MVA for Product-Form Networks.

Implements the linearizer family of approximate MVA methods for closed
queueing networks where exact MVA becomes computationally prohibitive.
Provides near-exact accuracy with reduced computational complexity.
"""

import numpy as np
from typing import Tuple, List, Optional, Union
from enum import Enum

# Scheduling strategy enum
class SchedStrategy(Enum):
    """Scheduling strategies for queueing stations."""
    INF = "INF"      # Infinite server (delay)
    FCFS = "FCFS"    # First Come First Serve
    PS = "PS"        # Processor Sharing
    LCFS = "LCFS"    # Last Come First Serve
    SIRO = "SIRO"    # Service In Random Order


def _oner(N: np.ndarray, indices: List[int]) -> np.ndarray:
    """
    Create population vector with one less job in specified classes.

    Args:
        N: Population vector
        indices: List of class indices to decrement (-1 means no decrement)

    Returns:
        Modified population vector
    """
    N_1 = N.copy().astype(float).ravel()
    for idx in indices:
        if idx >= 0 and idx < len(N_1):
            N_1[idx] = max(0, N_1[idx] - 1)
    return N_1


def pfqn_linearizer(
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray,
    sched_type: List[str],
    tol: float = 1e-8,
    maxiter: int = 1000,
    QN0: np.ndarray = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Linearizer approximate MVA algorithm.

    Args:
        L: Demand matrix (M x R) - rows are stations, columns are classes
        N: Population vector (R,)
        Z: Think time vector (R,)
        sched_type: List of scheduling strategies per station ('FCFS', 'PS', etc.)
        tol: Convergence tolerance; 'cn' or NaN selects the published
            Linearizer termination test of Chandy and Neuse,
            Commun. ACM 25(2), 1982, p.129: each Core call stops when
            max_{i,r}|dQ(i,r)|/N_r falls below pfqn_cntol evaluated at
            the population Core is running at, rather than on the
            Frobenius norm of dQ. See pfqn_cntol.
        maxiter: Maximum iterations

    Returns:
        Tuple of (Q, U, W, T, C, X, iterations):
            Q: Queue lengths (M x R)
            U: Utilizations (M x R)
            W: Waiting times (M x R)
            T: Station throughputs (M x R)
            C: Response times (1 x R)
            X: Class throughputs (1 x R)
            iterations: Number of iterations
    """
    return pfqn_gflinearizer(L, N, Z, sched_type, tol, maxiter, alpha=1.0, QN0=QN0)


def pfqn_gflinearizer(
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray,
    sched_type: List[str],
    tol: float = 1e-8,
    maxiter: int = 1000,
    alpha: float = 1.0,
    QN0: np.ndarray = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    General-form linearizer approximate MVA.

    Args:
        L: Demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,)
        sched_type: List of scheduling strategies per station
        tol: Convergence tolerance; 'cn' or NaN selects the published
            Linearizer termination test of Chandy and Neuse,
            Commun. ACM 25(2), 1982, p.129: each Core call stops when
            max_{i,r}|dQ(i,r)|/N_r falls below pfqn_cntol evaluated at
            the population Core is running at, rather than on the
            Frobenius norm of dQ. See pfqn_cntol.
        maxiter: Maximum iterations
        alpha: Linearization parameter (scalar)

    Returns:
        Same as pfqn_linearizer
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()
    Z = np.asarray(Z, dtype=float).ravel()

    R = len(N)
    alpha_vec = np.full(R, alpha)

    return pfqn_egflinearizer(L, N, Z, sched_type, tol, maxiter, alpha_vec, QN0=QN0)


def pfqn_egflinearizer(
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray,
    sched_type: List[str],
    tol: float = 1e-8,
    maxiter: int = 1000,
    alpha: np.ndarray = None,
    QN0: np.ndarray = None,
    npasses: int = 3
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Extended general-form linearizer with class-specific parameters.

    Args:
        L: Demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,)
        sched_type: List of scheduling strategies per station
        tol: Convergence tolerance; 'cn' or NaN selects the published
            Linearizer termination test of Chandy and Neuse,
            Commun. ACM 25(2), 1982, p.129: each Core call stops when
            max_{i,r}|dQ(i,r)|/N_r falls below pfqn_cntol evaluated at
            the population Core is running at, rather than on the
            Frobenius norm of dQ. See pfqn_cntol.
        maxiter: Maximum iterations
        alpha: Class-specific linearization parameters (R,)
        npasses: Number of Delta refresh passes (3 is the Chandy-Neuse fixed
            rule; pfqn_scat passes 1)

    Returns:
        Tuple of (Q, U, W, T, C, X, iterations), where T is the throughput PER
        REFERENCE VISIT (X broadcast over the stations the class visits), not the
        per-station throughput: visits are folded into L here and cannot be
        recovered from it. MATLAB returns no T at all and has the caller build
        `T = V .* X`; the slot exists because three Python callers unpack seven
        values. See the comment at the return statement.
    """
    from .mva import pfqn_bs
    from .cntol import is_cntol

    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()
    Z = np.asarray(Z, dtype=float).ravel()

    # Carried as NaN so that the pfqn_bs warm-start below inherits the same
    # test; the cutoff itself is population-dependent and so is recomputed
    # inside each Core call rather than once here.
    cntest = is_cntol(tol)
    if cntest:
        tol = float('nan')

    M, R = L.shape

    if alpha is None:
        alpha = np.ones(R)
    else:
        alpha = np.asarray(alpha, dtype=float).ravel()

    # Handle empty or zero demand case
    if L.size == 0 or np.max(L) == 0:
        X = N / Z
        Q = np.zeros((M, R))
        U = np.zeros((M, R))
        W = np.zeros((M, R))
        T = np.zeros((M, R))
        C = np.zeros(R)
        return Q, U, W, T, C, X.reshape(1, -1), 0

    # Initialize Q arrays for each station
    # Q[i] is (R x (1+R)) matrix: Q[i][r, s+1] = queue length of class r when one job of class s is removed
    Q = [np.zeros((R, 1 + R)) for _ in range(M)]
    Delta = [np.zeros((R, R)) for _ in range(M)]

    # Initial estimates using balanced system
    # pfqn_bs returns (XN, QN, UN, RN, it) - QN is at index 1
    for s in range(-1, R):
        N_1 = _oner(N, [s])
        if QN0 is None:
            init_result = pfqn_bs(L, N_1, Z)
        else:
            # warm-start the Bard-Schweitzer initialization from the supplied Q
            init_result = pfqn_bs(L, N_1, Z, tol, maxiter, QN0)
        Q_init = init_result[1]  # QN from pfqn_bs
        for i in range(M):
            for r in range(R):
                Q[i][r, 1 + s] = Q_init[i, r]

    totiter = 0

    # Main outer loop (3 refresh passes for linearization, 1 for SCAT)
    for I in range(npasses):
        for s in range(-1, R):
            N_1 = _oner(N, [s])
            Q1 = np.zeros((M, R))
            for i in range(M):
                for j in range(R):
                    Q1[i, j] = Q[i][j, 1 + s]

            Q_new, W_new, T_new, iter_count = _egflinearizer_core(
                L, M, R, N_1, Z, Q1, Delta, sched_type, tol, maxiter - totiter, alpha,
                cntest
            )

            for i in range(M):
                for j in range(R):
                    Q[i][j, 1 + s] = Q_new[i, j]

            totiter += iter_count

        # Update delta
        for i in range(M):
            for r in range(R):
                if N[r] == 1:
                    # see _kb/03-api-layer.md for rationale
                    Q[i][r, 1 + r] = 0
                # see _kb/03-api-layer.md for rationale
                for s in range(R):
                    Ns = _oner(N, [s])
                    if Ns[r] > 0:
                        Delta[i][r, s] = (
                            Q[i][r, 1 + s] / np.power(Ns[r], alpha[r]) -
                            Q[i][r, 0] / np.power(N[r], alpha[r])
                        )
                    else:
                        # see _kb/03-api-layer.md for rationale
                        Delta[i][r, s] = -Q[i][r, 0] / np.power(N[r], alpha[r])

    # Final core iteration with full population
    Q1 = np.zeros((M, R))
    for i in range(M):
        for j in range(R):
            Q1[i, j] = Q[i][j, 0]

    Q_final, W, T, iter_count = _egflinearizer_core(
        L, M, R, N, Z, Q1, Delta, sched_type, tol, maxiter - totiter, alpha, cntest
    )
    totiter += iter_count

    # Compute performance metrics
    X = T.copy()
    U = np.zeros((M, R))
    for i in range(M):
        for r in range(R):
            U[i, r] = X[r] * L[i, r]

    C = np.zeros(R)
    for r in range(R):
        if X[r] > 0:
            C[r] = N[r] / X[r] - Z[r]

    # THE FOURTH SLOT IS PER-REFERENCE-VISIT THROUGHPUT, NOT UTILIZATION. It used
    # to return `U` a second time, so every caller that wrote it into a TN matrix
    # reported the utilization in the throughput column -- and, downstream,
    # `sn_get_arvr_from_tput` propagated that through the routing matrix, so ArvR
    # came out as the utilization of the PREDECESSOR station. Invisible on a
    # product-form model (which takes `sn_deaggregate_chain_results` instead) and
    # on any model whose service is class-independent; reached by a HETEROGENEOUS
    # FCFS station, which is exactly where the linearizer family is needed.
    #
    # MATLAB has no fourth output at all -- `pfqn_egflinearizer.m` returns
    # [Q,U,W,C,X,totiter] and `solver_amva.m:291` builds `T = V .* repmat(X,M,1)`
    # in the CALLER, from visits this layer does not have. The value returned here
    # is therefore X broadcast to the stations the class visits, i.e. the
    # throughput PER REFERENCE VISIT; a caller reporting per-station throughput
    # must still scale it by V, as the callers in solver_mva.py now do.
    Tref = np.zeros((M, R))
    for i in range(M):
        for r in range(R):
            if L[i, r] > 0:
                Tref[i, r] = X[r]

    return Q_final, U, W, Tref, C.reshape(1, -1), X.reshape(1, -1), totiter


def _egflinearizer_core(
    L: np.ndarray,
    M: int,
    R: int,
    N_1: np.ndarray,
    Z: np.ndarray,
    Q: np.ndarray,
    Delta: List[np.ndarray],
    sched_type: List[str],
    tol: float,
    maxiter: int,
    alpha: np.ndarray,
    cntest: bool = False
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Core iteration for extended general-form linearizer.

    Returns:
        Tuple of (Q, W, T, iterations)
    """
    Q = Q.copy()
    W = L.copy()
    T = None
    has_converged = False
    iter_count = 0

    if cntest:
        from .cntol import pfqn_cntol
        # Chandy and Neuse (1982), p.129 and appendix: the cutoff is a function
        # of the population Core is running at, so it is recomputed here rather
        # than once for the whole Linearizer.
        tol = pfqn_cntol(N_1)
        nz = np.asarray(N_1, dtype=float).ravel() > 0

    while not has_converged:
        Q_last = Q.copy()

        # Estimate population
        Q_1 = _egflinearizer_estimate(L, M, R, N_1, Z, Q, Delta, W, alpha)

        # Forward MVA
        Q, W, T = _egflinearizer_forward_mva(L, M, R, sched_type, N_1, Z, Q_1)

        # Check convergence
        if cntest:
            # max_{i,r} |dQ(i,r)| / N_r over the non-empty classes; an empty
            # class would divide by zero and it carries no jobs to converge.
            if not np.any(nz):
                diff_norm = 0.0
            else:
                diff_norm = float(np.max(np.abs(Q[:, nz] - Q_last[:, nz]) / np.asarray(N_1, dtype=float).ravel()[nz]))
        else:
            diff_norm = np.linalg.norm(Q - Q_last)
        if diff_norm < tol or iter_count > maxiter:
            has_converged = True

        iter_count += 1

    return Q, W, T, iter_count


def _egflinearizer_estimate(
    L: np.ndarray,
    M: int,
    R: int,
    N_1: np.ndarray,
    Z: np.ndarray,
    Q: np.ndarray,
    Delta: List[np.ndarray],
    W: np.ndarray,
    alpha: np.ndarray
) -> List[np.ndarray]:
    """
    Estimate intermediate queue lengths for linearizer.

    Returns:
        List of Q_1 matrices for each station
    """
    Q_1 = [np.zeros((R, 1 + R)) for _ in range(M)]

    for i in range(M):
        for r in range(R):
            for s in range(R):
                Ns = _oner(N_1, [s])
                if N_1[r] > 0 and Ns[r] > 0:
                    Q_1[i][r, 1 + s] = (
                        np.power(Ns[r], alpha[r]) *
                        (Q[i, r] / np.power(N_1[r], alpha[r]) + Delta[i][r, s])
                    )

    return Q_1


def _egflinearizer_forward_mva(
    L: np.ndarray,
    M: int,
    R: int,
    sched_type: List[str],
    N_1: np.ndarray,
    Z: np.ndarray,
    Q_1: List[np.ndarray]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Forward MVA step for linearizer.

    Returns:
        Tuple of (Q, W, T)
    """
    W = np.zeros((M, R))
    T = np.zeros(R)
    Q = np.zeros((M, R))

    # see _kb/03-api-layer.md for rationale
    for i in range(M):
        for r in range(R):
            sum_Q = np.sum(Q_1[i][:, 1 + r])
            W[i, r] = L[i, r] * (1 + sum_Q)

    # Compute throughputs and queue lengths
    for r in range(R):
        W_col_sum = np.sum(W[:, r])
        if Z[r] + W_col_sum > 0:
            T[r] = N_1[r] / (Z[r] + W_col_sum)
        for i in range(M):
            Q[i, r] = T[r] * W[i, r]

    return Q, W, T


__all__ = [
    'pfqn_linearizer',
    'pfqn_gflinearizer',
    'pfqn_egflinearizer',
    'SchedStrategy',
]
