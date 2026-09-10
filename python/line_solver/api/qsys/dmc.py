"""
D/M/c queue analysis (deterministic interarrivals, exponential service).

Embeds at arrival epochs: X_n = number in system just before n-th arrival.
After arrival, Y = X_n + 1 jobs are present. During the deterministic interval
of length s = 1/lambda, the system evolves as a death-only CTMC with state-
dependent rate min(K, c) * mu.

Returns time-average performance metrics. Time-average L is obtained by
integrating E[N(t)] over a cycle, weighted by the arrival-epoch distribution.
"""
from typing import Dict
import numpy as np
from scipy.linalg import expm


def qsys_dmc(
    lambda_arr: float,
    mu: float,
    c: int,
    truncation: int = -1,
    quad_steps: int = 200,
) -> Dict:
    """Solve D/M/c via embedded DTMC at arrival epochs.

    Args:
        lambda_arr: Arrival rate (deterministic interarrivals of mean 1/lambda).
        mu: Service rate per server.
        c: Number of servers.
        truncation: Optional state-space truncation. Auto-chosen if <=0.
        quad_steps: Number of trapezoidal-rule steps for cycle integration.

    Returns:
        dict with mean_queue_length, mean_waiting_queue, mean_waiting_time,
        mean_sojourn_time, utilization.
    """
    if lambda_arr <= 0:
        raise ValueError("Arrival rate must be positive")
    if mu <= 0:
        raise ValueError("Service rate must be positive")
    if c < 1:
        raise ValueError("Number of servers must be >= 1")

    rho = lambda_arr / (c * mu)
    if rho >= 1 - 1e-12:
        raise ValueError(f"Load rho={rho} must be strictly less than 1")

    s = 1.0 / lambda_arr

    if truncation > 0:
        n_max = truncation
    else:
        n_max = max(200, min(2500, int(15.0 / (1.0 - rho)) + 200))

    n = n_max + 1

    # Death-only sub-generator: rate min(state, c) * mu
    A = np.zeros((n, n))
    for m in range(n):
        rate = min(m, c) * mu
        A[m, m] = -rate
        if m > 0:
            A[m, m - 1] = rate

    expAs = expm(A * s)

    # DTMC at arrival epochs: P[X, e] = expAs[X+1, e]
    P = np.zeros((n, n))
    Y_idx = np.minimum(np.arange(n) + 1, n - 1)
    for X in range(n):
        P[X, :] = expAs[Y_idx[X], :]

    M = (P.T - np.eye(n))
    M[-1, :] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0
    pi_arr = np.linalg.solve(M, b)

    # Time-average via integration of E[N(t)] and E[(N-c)+ at t] over a cycle.
    expAdt = expm(A * (s / quad_steps))
    weights_lq = np.maximum(np.arange(n) - c, 0).astype(float)
    weights_n = np.arange(n, dtype=float)

    Lq_at_t = np.zeros((n, quad_steps + 1))
    N_at_t = np.zeros((n, quad_steps + 1))
    expAt = np.eye(n)
    for k in range(quad_steps + 1):
        if k > 0:
            expAt = expAt @ expAdt
        Lq_at_t[:, k] = expAt @ weights_lq
        N_at_t[:, k] = expAt @ weights_n

    ts = np.linspace(0.0, s, quad_steps + 1)
    Lq_int = np.trapezoid(Lq_at_t, ts, axis=1) / s
    N_int = np.trapezoid(N_at_t, ts, axis=1) / s

    Lq_time = float((pi_arr * Lq_int[Y_idx]).sum())
    N_time = float((pi_arr * N_int[Y_idx]).sum())
    Wq = Lq_time / lambda_arr
    W = Wq + 1.0 / mu

    return {
        "mean_queue_length": N_time,
        "mean_waiting_queue": Lq_time,
        "mean_waiting_time": Wq,
        "mean_sojourn_time": W,
        "utilization": rho,
    }
