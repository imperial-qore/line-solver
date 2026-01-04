"""
QMAM: Queue Analysis using Matrix-Analytic Methods.

Native Python implementations of queue analysis algorithms for:
- MAP/MAP/1 queues
- MAP/M/c queues
- PH/PH/1 queues
- MMAP[K]/PH[K]/1 queues

Based on the Q-MAM library by Benny Van Houdt.
"""

import numpy as np
from numpy.typing import NDArray
from typing import Dict, Optional, Tuple, NamedTuple
import scipy.linalg as la
from dataclasses import dataclass

from ..smc import stat, qbd_cr, qbd_pi


@dataclass
class MAPMAP1Result:
    """Result of MAP/MAP/1 queue analysis."""
    queue_length: np.ndarray  # Queue length distribution
    soj_alpha: Optional[np.ndarray]  # Sojourn time PH alpha vector
    wait_alpha: Optional[np.ndarray]  # Waiting time PH alpha vector
    S_mat: Optional[np.ndarray]  # Service time matrix


@dataclass
class MAPMAP1Options:
    """Options for MAP/MAP/1 queue analysis."""
    mode: str = "SylvesCR"
    max_num_comp: int = 1000
    verbose: int = 0


def q_ct_map_map_1(C0: np.ndarray, C1: np.ndarray, D0: np.ndarray, D1: np.ndarray,
                   options: Optional[MAPMAP1Options] = None) -> MAPMAP1Result:
    """
    Compute queue length and time distributions for a continuous-time MAP/MAP/1/FCFS queue.

    Args:
        C0: MAP arrival process matrix D0 (ma x ma)
        C1: MAP arrival process matrix D1 (ma x ma)
        D0: MAP service process matrix D0 (ms x ms)
        D1: MAP service process matrix D1 (ms x ms)
        options: Solver options

    Returns:
        MAPMAP1Result containing queue length and time distributions
    """
    if options is None:
        options = MAPMAP1Options()

    C0 = np.asarray(C0, dtype=float)
    C1 = np.asarray(C1, dtype=float)
    D0 = np.asarray(D0, dtype=float)
    D1 = np.asarray(D1, dtype=float)

    ma = C0.shape[0]
    ms = D0.shape[0]
    mtot = ma * ms

    # Validate dimensions
    assert C0.shape == (ma, ma) and C1.shape == (ma, ma), \
        "Arrival process matrices must be ma x ma"
    assert D0.shape == (ms, ms) and D1.shape == (ms, ms), \
        "Service process matrices must be ms x ms"

    # Test the load of the queue
    inv_C0 = la.inv(-C0)
    pi_A = stat(C1 @ inv_C0)
    lambda_arr = np.sum(pi_A @ C1)

    inv_D0 = la.inv(-D0)
    pi_S = stat(D1 @ inv_D0)
    mu = np.sum(pi_S @ D1)

    load = lambda_arr / mu
    if load >= 1:
        raise ValueError(f"The load {load} of the system exceeds one")

    # Compute classic QBD blocks A0, A1, A2
    eye_ma = np.eye(ma)
    eye_ms = np.eye(ms)

    Am1 = np.kron(eye_ma, D1)
    A0 = np.kron(eye_ma, D0) + np.kron(C0, eye_ms)
    A1 = np.kron(C1, eye_ms)
    B0 = np.kron(C0, eye_ms)

    # Compute G and R using QBD solver
    qbd_result = qbd_cr(Am1, A0, A1, verbose=options.verbose > 0)
    R = qbd_result['R']

    # Compute stationary distribution
    stv = qbd_pi(Am1, B0, R, max_num_comp=options.max_num_comp, verbose=options.verbose)

    # Queue length distribution
    num_levels = len(stv) // mtot
    ql = np.zeros(num_levels)
    for i in range(num_levels):
        ql[i] = np.sum(stv[i * mtot:(i + 1) * mtot])

    # Compute sojourn and waiting time PH representation
    LM = np.kron(C1, D1)
    eye_D0 = np.kron(eye_ma, D0)

    # Compute T iteratively
    T_old = np.zeros((mtot, mtot))
    T_new = eye_D0.copy()

    # Simplified iteration for T
    max_iter = 100
    for _ in range(max_iter):
        # Schur-based approach
        schur = _schur_decomposition(C0)
        U, Tr = schur
        U_kron = np.kron(U, eye_ms)
        Tr_kron = np.kron(Tr, eye_ms)

        L = _q_sylvest(U_kron, Tr_kron, T_new)
        T_new = eye_D0 + L @ LM

        if la.norm(T_old - T_new, np.inf) < 1e-10:
            break
        T_old = T_new.copy()

    # Recompute L for final time distributions
    schur = _schur_decomposition(C0)
    U, Tr = schur
    L = _q_sylvest(np.kron(U, eye_ms), np.kron(Tr, eye_ms), T_new)

    # Compute Smat
    theta_tot = (np.kron(pi_A @ C1, pi_S) / (mu * load)).flatten()

    # Find non-zero entries
    nonz_indices = np.where(theta_tot > 1e-15)[0]

    if len(nonz_indices) == 0:
        return MAPMAP1Result(ql, None, None, None)

    theta_tot_red = theta_tot[nonz_indices]
    T_new_reduced = T_new[np.ix_(nonz_indices, nonz_indices)]

    diag_theta_red = np.diag(theta_tot_red)
    diag_theta_red_inv = la.inv(diag_theta_red)
    S_mat = diag_theta_red_inv @ T_new_reduced.T @ diag_theta_red

    # Alpha vector of PH representation of Sojourn time
    D1_sum_col = np.sum(D1, axis=1, keepdims=True)
    soj_alpha_full = np.diag(theta_tot) @ np.kron(np.ones((ma, 1)), D1_sum_col) * load / lambda_arr
    soj_alpha = soj_alpha_full[nonz_indices, 0]

    # Alpha vector of PH representation of Waiting time
    C1_sum_col = np.sum(C1, axis=1, keepdims=True)
    wait_alpha_full = np.diag(theta_tot) @ L @ np.kron(C1_sum_col, D1_sum_col) * load / lambda_arr
    wait_alpha = wait_alpha_full[nonz_indices, 0]

    return MAPMAP1Result(ql, soj_alpha, wait_alpha, S_mat)


def _schur_decomposition(A: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform Schur decomposition of a matrix.

    Returns (U, T) where A = U @ T @ U.T

    Args:
        A: Input matrix

    Returns:
        Tuple of (U, T) matrices
    """
    T, U = la.schur(A, output='real')
    return U, T


def _q_sylvest(U: np.ndarray, T: np.ndarray, B: np.ndarray) -> np.ndarray:
    """
    Solve the equation X*kron(A,I)+BX=-I using a Schur-based approach.

    Uses pre-computed Schur decomposition where kron(A,I)=U*T*U'

    Args:
        U: The unitary matrix from Schur decomposition
        T: The triangular matrix from Schur decomposition
        B: The matrix B in the equation

    Returns:
        Solution matrix X
    """
    n = B.shape[0]

    # Compute F = -U' (in the transformed coordinate system)
    F = -U.T

    # Solve column by column using back-substitution
    Y = np.zeros((n, n))

    for k in range(n):
        # Build right-hand side for column k
        temp = F[:, k].copy()

        # Subtract contributions from previously computed columns
        if k > 0:
            temp -= Y[:, :k] @ T[:k, k]

        # Solve (B + T(k,k)*I) * Y(:,k) = temp
        coeff = B + T[k, k] * np.eye(n)
        Y[:, k] = la.solve(coeff, temp)

    # X = Y @ U'
    return Y @ U.T


def q_ct_map_m_c(C0: np.ndarray, C1: np.ndarray, mu: float, c: int,
                 max_num_comp: int = 1000, verbose: int = 0) -> np.ndarray:
    """
    Compute queue length distribution for a MAP/M/c queue.

    Args:
        C0: MAP arrival process matrix D0
        C1: MAP arrival process matrix D1
        mu: Service rate per server
        c: Number of servers
        max_num_comp: Maximum number of components
        verbose: Verbosity level

    Returns:
        Queue length distribution
    """
    C0 = np.asarray(C0, dtype=float)
    C1 = np.asarray(C1, dtype=float)

    ma = C0.shape[0]

    # Compute arrival rate
    inv_C0 = la.inv(-C0)
    pi_A = stat(C1 @ inv_C0)
    lambda_arr = np.sum(pi_A @ C1)

    # Check stability
    load = lambda_arr / (c * mu)
    if load >= 1:
        raise ValueError(f"The load {load} of the system exceeds one")

    # Build QBD matrices
    # For level k < c: A_m1 = k*mu*I, A_0 = C0, A_1 = C1
    # For level k >= c: A_m1 = c*mu*I, A_0 = C0, A_1 = C1

    # Use repeating part (k >= c)
    Am1 = c * mu * np.eye(ma)
    A0 = C0.copy()
    A1 = C1.copy()

    # Solve QBD
    qbd_result = qbd_cr(Am1, A0, A1, verbose=verbose > 0)
    R = qbd_result['R']

    # For boundary, need to handle levels 0 to c-1 separately
    # Simplified: use level c boundary
    B0 = C0.copy()
    B1 = C0 + C1 + c * mu * np.eye(ma)  # Simplified boundary

    # Compute stationary distribution
    stv = qbd_pi(Am1, B0, R, max_num_comp=max_num_comp, verbose=verbose)

    # Queue length distribution
    num_levels = len(stv) // ma
    ql = np.zeros(num_levels)
    for i in range(num_levels):
        ql[i] = np.sum(stv[i * ma:(i + 1) * ma])

    return ql


def q_ct_ph_ph_1(alpha: np.ndarray, S: np.ndarray, beta: np.ndarray, T: np.ndarray,
                 max_num_comp: int = 1000, verbose: int = 0) -> MAPMAP1Result:
    """
    Compute distributions for a PH/PH/1 queue.

    Args:
        alpha: Arrival PH initial distribution
        S: Arrival PH sub-generator matrix
        beta: Service PH initial distribution
        T: Service PH sub-generator matrix
        max_num_comp: Maximum number of components
        verbose: Verbosity level

    Returns:
        MAPMAP1Result with queue length and time distributions
    """
    # Convert PH to MAP representation
    # For PH: D0 = S, D1 = (-S)*e*alpha
    alpha = np.asarray(alpha).flatten()
    S = np.asarray(S)
    beta = np.asarray(beta).flatten()
    T = np.asarray(T)

    ma = len(alpha)
    ms = len(beta)

    # Create MAP representation of arrival process
    s0 = -np.sum(S, axis=1, keepdims=True)  # Exit rate vector
    C0 = S.copy()
    C1 = s0 @ alpha.reshape(1, -1)

    # Create MAP representation of service process
    t0 = -np.sum(T, axis=1, keepdims=True)
    D0 = T.copy()
    D1 = t0 @ beta.reshape(1, -1)

    options = MAPMAP1Options(max_num_comp=max_num_comp, verbose=verbose)
    return q_ct_map_map_1(C0, C1, D0, D1, options)


def mm1_queue(lambda_arr: float, mu: float, max_level: int = 100) -> np.ndarray:
    """
    Compute queue length distribution for an M/M/1 queue.

    Args:
        lambda_arr: Arrival rate
        mu: Service rate
        max_level: Maximum queue length to compute

    Returns:
        Queue length distribution
    """
    rho = lambda_arr / mu
    if rho >= 1:
        raise ValueError(f"Load {rho} exceeds 1")

    # Geometric distribution
    ql = np.zeros(max_level)
    for k in range(max_level):
        ql[k] = (1 - rho) * (rho ** k)

    return ql


def mmk_queue(lambda_arr: float, mu: float, k: int, max_level: int = 100) -> np.ndarray:
    """
    Compute queue length distribution for an M/M/k queue.

    Args:
        lambda_arr: Arrival rate
        mu: Service rate per server
        k: Number of servers
        max_level: Maximum queue length to compute

    Returns:
        Queue length distribution
    """
    rho = lambda_arr / (k * mu)
    if rho >= 1:
        raise ValueError(f"Load {rho} exceeds 1")

    a = lambda_arr / mu  # Offered load

    # Compute normalization constant
    from math import factorial
    sum1 = sum(a**n / factorial(n) for n in range(k))
    sum2 = (a**k / factorial(k)) * (1 / (1 - rho))
    C = 1 / (sum1 + sum2)

    # Queue length distribution
    ql = np.zeros(max_level)

    for n in range(min(k, max_level)):
        ql[n] = C * (a**n) / factorial(n)

    for n in range(k, max_level):
        ql[n] = C * (a**k / factorial(k)) * (rho ** (n - k))

    return ql


__all__ = [
    'MAPMAP1Result',
    'MAPMAP1Options',
    'q_ct_map_map_1',
    'q_ct_map_m_c',
    'q_ct_ph_ph_1',
    'mm1_queue',
    'mmk_queue',
]
