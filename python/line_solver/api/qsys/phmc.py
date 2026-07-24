"""
PH/M/c queue analysis (phase-type inter-arrivals, exponential service, c servers).

Solves the GI/M/c QBD via Neuts' matrix-geometric method:
    R^2 A2 + R A1 + A0 = 0
where
    A0 = (-T 1) alpha,  A1 = T - c mu I,  A2 = c mu I.
The stationary level-n vector follows pi_n = pi_c R^{n-c} for n >= c, with
boundary states pi_0..pi_c determined by linear balance + normalization.

Returns time-average performance metrics.
"""
from typing import Dict
import numpy as np


def qsys_phmc(alpha, T, mu: float, c: int, max_iter: int = 50000, tol: float = 1e-14) -> Dict:
    """Solve PH/M/c via matrix-geometric.

    Args:
        alpha: PH entry probability vector (length k).
        T: PH sub-generator matrix (k x k), row sums <= 0.
        mu: Exponential service rate per server (>0).
        c: Number of servers (>=1).
        max_iter: Maximum iterations for R fixed-point.
        tol: Convergence tolerance for R.

    Returns:
        dict with mean_queue_length, mean_waiting_queue, mean_waiting_time,
        mean_sojourn_time, utilization.
    """
    alpha_row = np.asarray(alpha, dtype=float).flatten()
    Tm = np.asarray(T, dtype=float)
    if Tm.ndim != 2 or Tm.shape[0] != Tm.shape[1]:
        raise ValueError("T must be square")
    k = Tm.shape[0]
    if alpha_row.size != k:
        raise ValueError("alpha length must match T dimension")
    if mu <= 0:
        raise ValueError("Service rate mu must be positive")
    if c < 1:
        raise ValueError("Number of servers c must be >= 1")

    ones_k = np.ones(k)
    t_vec = -Tm @ ones_k
    D1 = np.outer(t_vec, alpha_row)

    try:
        mean_ia = float(alpha_row @ np.linalg.solve(-Tm, ones_k))
    except np.linalg.LinAlgError as e:
        raise ValueError(f"PH sub-generator T is singular: {e}")
    if mean_ia <= 0:
        raise ValueError(f"Non-positive mean inter-arrival: {mean_ia}")
    lam = 1.0 / mean_ia
    rho = lam / (c * mu)
    if rho >= 1.0 - 1e-12:
        raise ValueError(f"Load rho={rho} must be strictly less than 1")

    A0 = D1
    A1 = Tm - c * mu * np.eye(k)
    A2 = c * mu * np.eye(k)
    R = np.zeros((k, k))
    for _ in range(max_iter):
        try:
            inv_term = np.linalg.inv(A1 + R @ A2)
        except np.linalg.LinAlgError:
            break
        R_new = -A0 @ inv_term
        if np.max(np.abs(R_new - R)) < tol:
            R = R_new
            break
        R = R_new

    # Boundary linear system: variables pi_0, pi_1, ..., pi_c (each row vec dim k)
    n_var = (c + 1) * k
    M = np.zeros((n_var, n_var))

    def block_col(n: int) -> int:
        return n * k

    # Level 0 balance: pi_0 T + pi_1 (mu I) = 0
    for j in range(k):
        for i in range(k):
            M[j, block_col(0) + i] += Tm[i, j]
        M[j, block_col(1) + j] += mu

    # Levels 1 .. c-1 balance: pi_{n-1} D1 + pi_n (T - n mu I) + pi_{n+1} ((n+1) mu I) = 0
    for n in range(1, c):
        eqrow_base = n * k
        A1n = Tm - n * mu * np.eye(k)
        for j in range(k):
            row = eqrow_base + j
            for i in range(k):
                M[row, block_col(n - 1) + i] += D1[i, j]
                M[row, block_col(n) + i] += A1n[i, j]
            M[row, block_col(n + 1) + j] += (n + 1) * mu

    # Level c balance: pi_{c-1} D1 + pi_c (T - c mu I + R (c mu I)) = 0
    A1c_plus_RA2 = Tm - c * mu * np.eye(k) + R * (c * mu)
    eqrow_base = c * k
    for j in range(k):
        row = eqrow_base + j
        for i in range(k):
            M[row, block_col(c - 1) + i] += D1[i, j]
            M[row, block_col(c) + i] += A1c_plus_RA2[i, j]

    # Replace last row with normalization
    IR = np.eye(k) - R
    try:
        sum_geom = np.linalg.solve(IR, ones_k)
    except np.linalg.LinAlgError as e:
        raise ValueError(f"I - R is singular (rho={rho}): {e}")
    norm_row = np.zeros(n_var)
    for n in range(c):
        for i in range(k):
            norm_row[block_col(n) + i] = 1.0
    for i in range(k):
        norm_row[block_col(c) + i] = sum_geom[i]
    M[-1, :] = norm_row
    b = np.zeros(n_var)
    b[-1] = 1.0

    x = np.linalg.solve(M, b)
    pis = [x[n * k:(n + 1) * k] for n in range(c + 1)]
    pi_c = pis[c]

    IR_inv = np.linalg.inv(IR)
    IR_inv2 = IR_inv @ IR_inv
    Lq = float(pi_c @ R @ IR_inv2 @ ones_k)
    L_bulk = float(pi_c @ (c * IR_inv + R @ IR_inv2) @ ones_k)
    L = sum(n * pis[n].sum() for n in range(c)) + L_bulk
    Wq = Lq / lam
    W = Wq + 1.0 / mu

    return {
        "mean_queue_length": L,
        "mean_waiting_queue": Lq,
        "mean_waiting_time": Wq,
        "mean_sojourn_time": W,
        "utilization": rho,
    }
