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
    maxiter: int = 1000
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Linearizer approximate MVA algorithm.

    Args:
        L: Demand matrix (M x R) - rows are stations, columns are classes
        N: Population vector (R,)
        Z: Think time vector (R,)
        sched_type: List of scheduling strategies per station ('FCFS', 'PS', etc.)
        tol: Convergence tolerance
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
    return pfqn_gflinearizer(L, N, Z, sched_type, tol, maxiter, alpha=1.0)


def pfqn_gflinearizer(
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray,
    sched_type: List[str],
    tol: float = 1e-8,
    maxiter: int = 1000,
    alpha: float = 1.0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    General-form linearizer approximate MVA.

    Args:
        L: Demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,)
        sched_type: List of scheduling strategies per station
        tol: Convergence tolerance
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

    return pfqn_egflinearizer(L, N, Z, sched_type, tol, maxiter, alpha_vec)


def pfqn_egflinearizer(
    L: np.ndarray,
    N: np.ndarray,
    Z: np.ndarray,
    sched_type: List[str],
    tol: float = 1e-8,
    maxiter: int = 1000,
    alpha: np.ndarray = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Extended general-form linearizer with class-specific parameters.

    Args:
        L: Demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,)
        sched_type: List of scheduling strategies per station
        tol: Convergence tolerance
        maxiter: Maximum iterations
        alpha: Class-specific linearization parameters (R,)

    Returns:
        Tuple of (Q, U, W, T, C, X, iterations)
    """
    from .mva import pfqn_bs

    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).ravel()
    Z = np.asarray(Z, dtype=float).ravel()

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
    for s in range(-1, R):
        N_1 = _oner(N, [s])
        init_result = pfqn_bs(L, N_1, Z)
        Q_init = init_result[2]  # QN from pfqn_bs
        for i in range(M):
            for r in range(R):
                Q[i][r, 1 + s] = Q_init[i, r]

    totiter = 0

    # Main outer loop (3 iterations for linearization)
    for I in range(3):
        for s in range(-1, R):
            N_1 = _oner(N, [s])
            Q1 = np.zeros((M, R))
            for i in range(M):
                for j in range(R):
                    Q1[i, j] = Q[i][j, 1 + s]

            Q_new, W_new, T_new, iter_count = _egflinearizer_core(
                L, M, R, N_1, Z, Q1, Delta, sched_type, tol, maxiter - totiter, alpha
            )

            for i in range(M):
                for j in range(R):
                    Q[i][j, 1 + s] = Q_new[i, j]

            totiter += iter_count

        # Update delta
        for i in range(M):
            for r in range(R):
                if N[r] == 1:
                    for j in range(R):
                        Q[i][j, 1 + r] = 0
                for s in range(R):
                    if N[s] > 1:
                        Ns = _oner(N, [s])
                        if N[r] > 0 and Ns[r] > 0:
                            Delta[i][r, s] = (
                                Q[i][r, 1 + s] / np.power(Ns[r], alpha[r]) -
                                Q[i][r, 0] / np.power(N[r], alpha[r])
                            )

    # Final core iteration with full population
    Q1 = np.zeros((M, R))
    for i in range(M):
        for j in range(R):
            Q1[i, j] = Q[i][j, 0]

    Q_final, W, T, iter_count = _egflinearizer_core(
        L, M, R, N, Z, Q1, Delta, sched_type, tol, maxiter - totiter, alpha
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

    return Q_final, U, W, U, C.reshape(1, -1), X.reshape(1, -1), totiter


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
    alpha: np.ndarray
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

    while not has_converged:
        Q_last = Q.copy()

        # Estimate population
        Q_1 = _egflinearizer_estimate(L, M, R, N_1, Z, Q, Delta, W, alpha)

        # Forward MVA
        Q, W, T = _egflinearizer_forward_mva(L, M, R, sched_type, N_1, Z, Q_1)

        # Check convergence
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

    # Compute residence times
    for i in range(M):
        for r in range(R):
            sched = sched_type[i].upper() if isinstance(sched_type[i], str) else sched_type[i]

            if sched == 'FCFS':
                W[i, r] = L[i, r]
                if L[i, r] != 0:
                    for s in range(R):
                        W[i, r] += L[i, s] * Q_1[i][s, 1 + r]
            else:
                # PS or other strategies
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
