"""
M/D/c queue analysis (Poisson arrivals, deterministic service) via Crommelin.

Solves the embedded DTMC at multiples of the service time:
    X_{n+1} = max(0, X_n - c) + A_n,    A_n ~ Poisson(lambda * s)

This gives the time-average queue-length distribution (PASTA) and is exact
for M/D/c with FCFS.
"""
from typing import Dict
import numpy as np
from scipy.special import gammaln


def qsys_mdc_crommelin(
    lambda_arr: float,
    s: float,
    c: int,
    truncation: int = -1,
) -> Dict:
    """Solve M/D/c via Crommelin's embedded DTMC.

    Args:
        lambda_arr: Poisson arrival rate.
        s: Deterministic service time (>0).
        c: Number of servers (>=1).
        truncation: Optional state-space truncation level. If <=0, chosen
            automatically to keep dense LU manageable.

    Returns:
        dict with keys:
            mean_queue_length:  E[N] (number in system)
            mean_waiting_queue: Lq = E[(N-c)+]
            mean_waiting_time:  Wq = Lq / lambda
            mean_sojourn_time:  W = Wq + s
            utilization:        rho = lambda * s / c
    """
    if lambda_arr <= 0:
        raise ValueError("Arrival rate must be positive")
    if s <= 0:
        raise ValueError("Service time must be positive")
    if c < 1:
        raise ValueError("Number of servers must be >= 1")

    a = lambda_arr * s
    rho = a / c
    if rho >= 1 - 1e-12:
        raise ValueError(f"Load rho={rho} must be strictly less than 1")

    if truncation > 0:
        n_max = truncation
    else:
        # Empirically N ~ 10/(1-rho) gives 6+ digits for moderate c. Cap at
        # 2500 to keep dense LU manageable.
        n_max = max(200, min(2500, int(10.0 / (1.0 - rho)) + 200))

    n = n_max + 1
    # Poisson(a) PMF in log-space
    k_arr = np.arange(n)
    log_pmf = -a + k_arr * np.log(a) - gammaln(k_arr + 1)
    pmf = np.exp(log_pmf)

    # Build A = (P^T - I), with last row replaced by all-ones for normalization.
    # P[i, j] = pmf[j - max(0, i - c)] when j >= max(0, i - c), else 0.
    A = np.zeros((n, n))
    for i in range(n):
        base = 0 if i <= c else i - c
        length = n - base
        # Cap length to available pmf entries
        ll = min(length, n)
        A[base:base + ll, i] = pmf[:ll]
    A -= np.eye(n)
    A[-1, :] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0

    pi = np.linalg.solve(A, b)

    idx = np.arange(n)
    mean_n = float((idx * pi).sum())
    Lq = float(np.maximum(idx - c, 0).dot(pi))
    Wq = Lq / lambda_arr
    W = Wq + s

    return {
        "mean_queue_length": mean_n,
        "mean_waiting_queue": Lq,
        "mean_waiting_time": Wq,
        "mean_sojourn_time": W,
        "utilization": rho,
    }
