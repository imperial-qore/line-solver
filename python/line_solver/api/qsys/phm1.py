"""
PH/M/1 queue analysis (phase-type inter-arrivals, exponential service).

Solves the GI/M/1 fundamental equation
    sigma = psi_A(mu*(1 - sigma))
with the LST of a PH distribution
    psi_A(s) = alpha (s I - T)^{-1} (-T 1).

Returns time-average performance metrics:
    L_time  = rho / (1 - sigma)
    Lq_time = rho * sigma / (1 - sigma)
"""
from typing import Dict, Sequence
import numpy as np
from scipy.optimize import brentq


def qsys_phm1(alpha, T, mu: float) -> Dict:
    """Solve PH/M/1 via the GI/M/1 sigma-root.

    Args:
        alpha: PH entry probability vector (length k).
        T: PH sub-generator matrix (k x k), row sums <= 0.
        mu: Exponential service rate (>0).

    Returns:
        dict with mean_queue_length, mean_waiting_queue, mean_waiting_time,
        mean_sojourn_time, utilization, sigma.
    """
    alpha_row = np.asarray(alpha, dtype=float).flatten()
    Tm = np.asarray(T, dtype=float)
    if Tm.ndim != 2 or Tm.shape[0] != Tm.shape[1]:
        raise ValueError("T must be a square matrix")
    k = Tm.shape[0]
    if alpha_row.size != k:
        raise ValueError("alpha length must match T dimension")
    if mu <= 0:
        raise ValueError("Service rate mu must be positive")

    ones_k = np.ones(k)
    t_vec = -Tm @ ones_k  # absorption rates per phase

    # Mean inter-arrival time = alpha * (-T)^{-1} * 1
    try:
        mean_ia = float(alpha_row @ np.linalg.solve(-Tm, ones_k))
    except np.linalg.LinAlgError as e:
        raise ValueError(f"PH sub-generator T is singular: {e}")
    if mean_ia <= 0:
        raise ValueError(f"Non-positive mean inter-arrival time: {mean_ia}")
    lam = 1.0 / mean_ia
    rho = lam / mu
    if rho >= 1.0 - 1e-12:
        raise ValueError(f"Load rho={rho} must be strictly less than 1")

    def lst(s: float) -> float:
        return float(alpha_row @ np.linalg.solve(s * np.eye(k) - Tm, t_vec))

    def f(sigma: float) -> float:
        return sigma - lst(mu * (1.0 - sigma))

    # f(0) = -lst(mu) < 0; f(1-) -> 1 - 1 = 0 (trivial root). Search just below 1.
    f0 = f(1e-12)
    f1 = f(1.0 - 1e-12)
    if f0 * f1 > 0:
        # Fall back to fixed-point iteration if brentq cannot bracket
        sigma = 0.5
        for _ in range(2000):
            sigma_new = lst(mu * (1.0 - sigma))
            if abs(sigma_new - sigma) < 1e-13:
                sigma = sigma_new
                break
            sigma = sigma_new
    else:
        sigma = brentq(f, 1e-12, 1.0 - 1e-12, xtol=1e-12, rtol=1e-12)

    L = rho / (1.0 - sigma)
    Lq = rho * sigma / (1.0 - sigma)
    Wq = Lq / lam
    W = Wq + 1.0 / mu

    return {
        "mean_queue_length": L,
        "mean_waiting_queue": Lq,
        "mean_waiting_time": Wq,
        "mean_sojourn_time": W,
        "utilization": rho,
        "sigma": sigma,
    }
