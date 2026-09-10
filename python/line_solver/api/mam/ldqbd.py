"""
Level-Dependent Quasi-Birth-Death (LDQBD) process solver.

Native Python implementation for analyzing LDQBD processes using matrix
continued fractions. Implements Algorithm 1 from "A Simple Algorithm for
the Rate Matrices of Level-Dependent QBD Processes" by Phung-Duc,
Masuyama, Kasahara, and Takahashi (2010), QTNA Conference.

An LDQBD has a block-tridiagonal infinitesimal generator:
    Q^(0)_1  Q^(0)_0   O        O      ...   O
    Q^(1)_2  Q^(1)_1  Q^(1)_0   O      ...   O
    O       Q^(2)_2  Q^(2)_1  Q^(2)_0 ...   O
    ...

where Q_0^(n) are upward transitions, Q_1^(n) are local transitions,
and Q_2^(n) are downward transitions.
"""

import numpy as np
from scipy import linalg
from dataclasses import dataclass
from typing import List, Tuple, Optional


@dataclass
class LdqbdResult:
    """Result of LDQBD solver.

    Attributes:
        R: List of rate matrices R^(1), R^(2), ..., R^(N)
        pi: Stationary distribution vector [pi_0, pi_1, ..., pi_N]
        pi_cells: Per-level stationary vectors, still resolved by phase; pi is
            their row sums. A caller that needs a phase-conditional quantity --
            the mean of something that depends on the phase and not only on the
            level -- cannot recover it from the aggregated vector and must read
            these.
    """
    R: List[np.ndarray]
    pi: np.ndarray
    pi_cells: Optional[List[np.ndarray]] = None


@dataclass
class LdqbdOptions:
    """Options for LDQBD solver.

    Attributes:
        epsilon: Convergence tolerance (default: 1e-10)
        max_iter: Maximum number of iterations (default: 1000)
        verbose: Print debug information (default: False)
    """
    epsilon: float = 1e-10
    max_iter: int = 1000
    verbose: bool = False


def ldqbd(Q0: List[np.ndarray],
         Q1: List[np.ndarray],
         Q2: List[np.ndarray],
         options: Optional[LdqbdOptions] = None) -> LdqbdResult:
    """
    Solve a level-dependent QBD process.

    Args:
        Q0: List of upward transition matrices [Q0^(0), Q0^(1), ..., Q0^(N-1)]
            Q0^(n) has shape (states_n, states_{n+1})
        Q1: List of local transition matrices [Q1^(0), Q1^(1), ..., Q1^(N)]
            Q1^(n) has shape (states_n, states_n)
        Q2: List of downward transition matrices [Q2^(1), Q2^(2), ..., Q2^(N)]
            Q2^(n) has shape (states_n, states_{n-1})
        options: LdqbdOptions instance (uses defaults if None)

    Returns:
        LdqbdResult containing rate matrices R and stationary distribution pi

    Algorithm:
        1. Compute rate matrices backward using continued fraction recursion
        2. Compute stationary distribution forward using R matrices
    """
    if options is None:
        options = LdqbdOptions()

    # see _kb/03-api-layer.md for rationale
    N = len(Q1) - 1  # Maximum level

    if options.verbose:
        print(f"LD-QBD Solver: N={N} levels, epsilon={options.epsilon}")

    # Validate input dimensions
    if len(Q0) != N:
        raise ValueError(f"Q0 must have {N} matrices (levels 0 to {N-1})")
    if len(Q2) != N:
        raise ValueError(f"Q2 must have {N} matrices (levels 1 to {N})")

    # Compute all R matrices using backward recursion
    R = _compute_all_rate_matrices(N, Q0, Q1, Q2, options)

    if options.verbose:
        for n in range(N):
            print(f"  R^({n+1}) computed ({R[n].shape})")

    # Compute stationary distribution
    pi_cells = []
    pi = _compute_stationary_dist(R, Q0, Q1, Q2, N, options, pi_cells)

    return LdqbdResult(R, pi, pi_cells)


def _compute_all_rate_matrices(N: int,
                              Q0: List[np.ndarray],
                              Q1: List[np.ndarray],
                              Q2: List[np.ndarray],
                              options: LdqbdOptions) -> List[np.ndarray]:
    """
    Compute all R matrices using backward recursion.

    For heterogeneous LDQBD:
        R^(N) = Q0^(N-1) * (-Q1^(N))^{-1}
        R^(n) = Q0^(n-1) * (-Q1^(n) - R^(n+1) * Q2^(n+1))^{-1}  for n = N-1,...,1

    Args:
        N: Maximum level
        Q0, Q1, Q2: List of transition matrices
        options: Solver options

    Returns:
        List of R matrices [R^(1), R^(2), ..., R^(N)]
    """
    R = [None] * N

    # Start at level N (boundary): R^(N) = Q0^(N-1) * (-Q1^(N))^{-1}
    Q0_Nminus1 = np.asarray(Q0[N - 1], dtype=np.float64)
    Q1_N = np.asarray(Q1[N], dtype=np.float64)

    U_N = -Q1_N
    R[N - 1] = _safe_matrix_divide(Q0_Nminus1, U_N)

    # Backward recursion for n = N-1 down to 1
    for n in range(N - 1, 0, -1):
        Q0_nminus1 = np.asarray(Q0[n - 1], dtype=np.float64)
        Q1_n = np.asarray(Q1[n], dtype=np.float64)
        Q2_nplus1 = np.asarray(Q2[n], dtype=np.float64)  # This is Q2^(n+1)
        R_nplus1 = R[n]

        # Compute R^(n+1) * Q2^(n+1)
        RQ2 = R_nplus1 @ Q2_nplus1

        # Compute U = -Q1^(n) - R^(n+1) * Q2^(n+1)
        U = -Q1_n - RQ2

        # Compute R^(n) = Q0^(n-1) * U^{-1}
        R[n - 1] = _safe_matrix_divide(Q0_nminus1, U)

    return R


def _safe_matrix_divide(A: np.ndarray, U: np.ndarray) -> np.ndarray:
    """
    Safely compute A * U^{-1}.

    Uses pseudo-inverse if U is singular.

    Args:
        A: Numerator matrix
        U: Denominator matrix (to be inverted)

    Returns:
        Product A * U^{-1}
    """
    A = np.asarray(A, dtype=np.float64)
    U = np.asarray(U, dtype=np.float64)

    if U.shape == (1, 1):
        # Scalar case
        u_val = U[0, 0]
        if abs(u_val) > 1e-14:
            return A / u_val
        else:
            return np.zeros_like(A)
    else:
        # Matrix case
        try:
            # Try standard inversion first
            det = np.linalg.det(U)
            if abs(det) > 1e-14:
                U_inv = np.linalg.inv(U)
                return A @ U_inv
            else:
                # Use pseudo-inverse for singular matrices
                U_pinv = np.linalg.pinv(U)
                return A @ U_pinv
        except np.linalg.LinAlgError:
            # If inversion fails, use pseudo-inverse
            U_pinv = np.linalg.pinv(U)
            return A @ U_pinv


def _compute_stationary_dist(R: List[np.ndarray],
                            Q0: List[np.ndarray],
                            Q1: List[np.ndarray],
                            Q2: List[np.ndarray],
                            N: int,
                            options: LdqbdOptions,
                            pi_cells_out: Optional[List[np.ndarray]] = None) -> np.ndarray:
    """
    Compute stationary distribution from rate matrices.

    Algorithm:
        1. Boundary equation: pi_0 * (Q1^(0) + R^(1)*Q2^(1)) = 0
        2. Forward recursion: pi_n = pi_{n-1} * R^(n)
        3. Normalize: sum(pi) = 1

    Args:
        R: List of rate matrices
        Q0, Q1, Q2: Transition matrices
        N: Maximum level
        options: Solver options

    Returns:
        Stationary distribution vector [pi_0, pi_1, ..., pi_N]
    """
    # Check if all levels have same dimension (homogeneous) and scalar
    dims = [Q1_n.shape[0] for Q1_n in Q1]
    is_homogeneous = len(set(dims)) == 1
    is_scalar = all(d == 1 for d in dims)

    if is_scalar and is_homogeneous:
        # Scalar case: direct computation
        pi = np.zeros(N + 1)

        # Start with pi_0 = 1 (will normalize later)
        pi[0] = 1.0

        # Forward recursion: pi_n = pi_{n-1} * R^(n)
        for n in range(1, N + 1):
            if R[n - 1].shape == (1, 1):
                pi[n] = pi[n - 1] * R[n - 1][0, 0]
            else:
                pi[n] = 0.0

        # Normalize
        total = np.sum(pi)
        if total > 0:
            pi = pi / total

        if pi_cells_out is not None:
            for n in range(N + 1):
                pi_cells_out.append(np.array([[pi[n]]]))
        return pi
    else:
        # Heterogeneous or matrix case
        pi_cells = [np.ones((1, 1)) for _ in range(N + 1)]

        # Initialize pi_0 based on its dimension
        if Q1[0].shape == (1, 1):
            # Scalar level 0: start with pi_0 = 1
            pi_cells[0] = np.array([[1.0]])
        else:
            # Matrix level 0: solve boundary equation
            Q1_0 = np.asarray(Q1[0], dtype=np.float64)
            Q2_1 = np.asarray(Q2[0], dtype=np.float64)
            A = Q1_0 + (R[0] @ Q2_1)

            # Find left null space (solve pi * A = 0)
            pi0 = _solve_left_null_space(A)
            pi_cells[0] = pi0

        # Forward recursion: pi_n = pi_{n-1} * R^(n)
        for n in range(1, N + 1):
            pi_cells[n] = pi_cells[n - 1] @ R[n - 1]

        # Normalize and convert to scalar probabilities
        total = sum(np.sum(pi_cells[n]) for n in range(N + 1))

        pi = np.zeros(N + 1)
        for n in range(N + 1):
            if total > 0:
                pi_cells[n] = pi_cells[n] / total
                pi[n] = np.sum(pi_cells[n])
            else:
                pi[n] = 1.0 / (N + 1)

        if pi_cells_out is not None:
            pi_cells_out.extend(pi_cells)
        return pi


def _solve_left_null_space(A: np.ndarray) -> np.ndarray:
    """
    Boundary vector of the level-0 block: the solution of
    pi_0 * (Q1^(0) + R^(1) Q2^(1)) = 0.

    That matrix is the generator of the process censored on level 0, so its
    stationary distribution IS the boundary vector and ctmc_solve is the right
    instrument. Taking the eigenvector of the smallest-magnitude eigenvalue of
    A.T is not: eigenvalues near zero are not separated from the true null
    direction on a stiff block, and a complex conjugate pair returns a vector
    with no probabilistic meaning at all. No caller reached this branch before
    the bgchain method, whose level 0 carries the arrival and environment phases.

    Args:
        A: the censored generator

    Returns:
        Left null space vector (row vector)
    """
    from ..mc.ctmc import ctmc_solve, ctmc_makeinfgen
    A = np.asarray(A, dtype=np.float64)
    pi = np.asarray(ctmc_solve(ctmc_makeinfgen(A)), dtype=np.float64).reshape(1, -1)
    return pi
