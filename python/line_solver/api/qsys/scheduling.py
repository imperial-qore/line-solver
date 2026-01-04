"""
M/G/1 Queue Scheduling Discipline Analysis.

Native Python implementations of analytical formulas for M/G/1 queues
with various scheduling disciplines including priority, SRPT, feedback, etc.

Key functions:
    qsys_mg1_prio: Non-preemptive (Head-of-Line) priority scheduling
    qsys_mg1_srpt: Shortest Remaining Processing Time scheduling
    qsys_mg1_fb: Foreground-Background (LAS) scheduling
    qsys_mg1_lrpt: Longest Remaining Processing Time scheduling
    qsys_mg1_psjf: Preemptive Shortest Job First
    qsys_mg1_setf: Shortest Expected Time First

References:
    Original MATLAB: matlab/src/api/qsys/qsys_mg1_*.m
    Kleinrock, "Queueing Systems, Volume I: Theory", 1975
"""

import numpy as np
from typing import Tuple


def qsys_mg1_prio(lambda_vec: np.ndarray, mu_vec: np.ndarray,
                  cs_vec: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Analyze M/G/1 queue with non-preemptive (Head-of-Line) priorities.

    Computes per-class mean waiting times for multiple priority classes
    using the Pollaczek-Khinchine formula extended for priorities.

    For K priority classes (class 1 = highest priority), the waiting time
    for class k is:
        W_k = (B_0 / (1 - sum_{i=1}^{k-1} rho_i)) * (1 / (1 - sum_{i=1}^{k} rho_i))

    Args:
        lambda_vec: Vector of arrival rates per priority class (class 1 = highest)
        mu_vec: Vector of service rates per priority class
        cs_vec: Vector of coefficients of variation per priority class

    Returns:
        Tuple of (W, rho):
            W: Vector of mean response times per priority class
            rho: System utilization (rhohat = Q/(1+Q) format)

    Raises:
        ValueError: If inputs have different lengths or are not positive
        ValueError: If system is unstable (utilization >= 1)

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mg1_prio.m
        Kleinrock, "Queueing Systems, Volume I: Theory", Section 3.5
    """
    lambda_vec = np.asarray(lambda_vec, dtype=float).flatten()
    mu_vec = np.asarray(mu_vec, dtype=float).flatten()
    cs_vec = np.asarray(cs_vec, dtype=float).flatten()

    # Validate inputs
    if not (len(lambda_vec) == len(mu_vec) == len(cs_vec)):
        raise ValueError("lambda, mu, and cs must have the same length")

    if np.any(lambda_vec <= 0) or np.any(mu_vec <= 0) or np.any(cs_vec <= 0):
        raise ValueError("lambda, mu, and cs must all be positive")

    # Per-class utilizations
    rho_i = lambda_vec / mu_vec

    # Overall utilization
    rho_total = np.sum(rho_i)

    # Stability check
    if rho_total >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_total:.4g} >= 1")

    # Mean second moment of service time
    # B_0 = sum_i lambda_i * E[S_i^2] / 2
    # where E[S_i^2] = (1 + cs_i^2) / mu_i^2
    B_0 = np.sum(lambda_vec * (1 + cs_vec ** 2) / (mu_vec ** 2)) / 2

    # Compute per-class waiting times
    K = len(lambda_vec)
    W_q = np.zeros(K)

    for k in range(K):
        # Cumulative utilization for higher priority classes
        rho_prev = np.sum(rho_i[:k]) if k > 0 else 0.0

        # Cumulative utilization including current class
        rho_curr = np.sum(rho_i[:k + 1])

        # Waiting time in queue for class k
        W_q[k] = B_0 / ((1 - rho_prev) * (1 - rho_curr))

    # Response time = waiting time + service time
    W = W_q + 1.0 / mu_vec

    # Compute rhohat = Q/(1+Q) for consistency
    Q = np.sum(lambda_vec * W)
    rho_hat = Q / (1 + Q)

    return W, rho_hat


def qsys_mg1_srpt(lambda_vec: np.ndarray, mu_vec: np.ndarray,
                  cs_vec: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Analyze M/G/1 queue with Shortest Remaining Processing Time (SRPT).

    For SRPT, smaller jobs have preemptive priority over larger jobs.
    When all service time distributions are exponential (cs=1), SRPT
    reduces to preemptive priority based on service rates.

    Args:
        lambda_vec: Vector of arrival rates per class
        mu_vec: Vector of service rates per class
        cs_vec: Vector of coefficients of variation per class

    Returns:
        Tuple of (W, rho):
            W: Vector of mean response times per class (sorted by original order)
            rho: System utilization (rhohat format)

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mg1_srpt.m
        Schrage and Miller, "The queue M/G/1 with the shortest remaining
        processing time discipline", Operations Research, 1966
    """
    lambda_vec = np.asarray(lambda_vec, dtype=float).flatten()
    mu_vec = np.asarray(mu_vec, dtype=float).flatten()
    cs_vec = np.asarray(cs_vec, dtype=float).flatten()

    # Validate inputs
    if not (len(lambda_vec) == len(mu_vec) == len(cs_vec)):
        raise ValueError("lambda, mu, and cs must have the same length")

    if np.any(lambda_vec <= 0) or np.any(mu_vec <= 0) or np.any(cs_vec < 0):
        raise ValueError("lambda and mu must be positive, cs must be non-negative")

    K = len(lambda_vec)

    # Mean service times per class
    mean_service = 1.0 / mu_vec

    # Sort classes by mean service time (ascending) for SRPT priority
    sort_idx = np.argsort(mean_service)
    unsort_idx = np.argsort(sort_idx)

    # Reorder parameters
    lambda_sorted = lambda_vec[sort_idx]
    mu_sorted = mu_vec[sort_idx]
    cs_sorted = cs_vec[sort_idx]

    # Per-class utilizations
    rho_i = lambda_sorted / mu_sorted
    rho_total = np.sum(rho_i)

    # Stability check
    if rho_total >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_total:.4g} >= 1")

    # Check if all exponential (SRPT = preemptive priority for exponential)
    is_exponential = np.all(np.abs(cs_sorted - 1) < 1e-10)

    if is_exponential or K == 1:
        # Use preemptive priority formula
        W_sorted = _mg1_srpt_exp(lambda_sorted, mu_sorted)
    else:
        # General case - use approximation
        W_sorted = _mg1_srpt_general(lambda_sorted, mu_sorted, cs_sorted)

    # Restore original ordering
    W = W_sorted[unsort_idx]

    # Compute rhohat
    Q = np.sum(lambda_vec * W)
    rho_hat = Q / (1 + Q)

    return W, rho_hat


def _mg1_srpt_exp(lambda_vec: np.ndarray, mu_vec: np.ndarray) -> np.ndarray:
    """SRPT for exponential service - uses preemptive priority formula."""
    K = len(lambda_vec)
    rho_i = lambda_vec / mu_vec

    W = np.zeros(K)
    for k in range(K):
        rho_prev = np.sum(rho_i[:k]) if k > 0 else 0.0
        rho_curr = np.sum(rho_i[:k + 1])

        # Cumulative residual service from classes 1 to k
        E_R_k = np.sum(lambda_vec[:k + 1] / (mu_vec[:k + 1] ** 2))

        # Waiting time (preemptive priority)
        W_q = E_R_k / ((1 - rho_prev) * (1 - rho_curr))

        # Response time
        W[k] = W_q + 1.0 / mu_vec[k]

    return W


def _mg1_srpt_general(lambda_vec: np.ndarray, mu_vec: np.ndarray,
                      cs_vec: np.ndarray) -> np.ndarray:
    """SRPT for general service - uses numerical approximation."""
    K = len(lambda_vec)
    W = np.zeros(K)

    for k in range(K):
        x = 1.0 / mu_vec[k]  # Mean service time for class k

        # Truncated load
        smaller_idx = (1.0 / mu_vec) <= x
        rho_x = np.sum(lambda_vec[smaller_idx] / mu_vec[smaller_idx])

        if rho_x >= 1:
            W[k] = np.inf
            continue

        # Numerator integral approximation
        numerator = np.sum(lambda_vec[smaller_idx] / (mu_vec[smaller_idx] ** 2))

        # Second integral (numerical quadrature)
        n_steps = 1000
        dt = x / n_steps
        second_term = 0.0

        for step in range(1, n_steps + 1):
            t = (step - 0.5) * dt
            rho_t = 0.0
            for i in range(K):
                # Truncated load at time t
                rho_t += lambda_vec[i] / mu_vec[i] * (1 - np.exp(-mu_vec[i] * t) * (1 + mu_vec[i] * t))

            if rho_t < 1:
                second_term += dt / (1 - rho_t)

        # Schrage-Miller formula
        W[k] = numerator / (1 - rho_x) ** 2 + second_term

    return W


def qsys_mg1_fb(lambda_val: float, mu_val: float, cs_val: float) -> Tuple[float, float]:
    """
    Analyze M/G/1 queue with Foreground-Background (FB/LAS) scheduling.

    In FB scheduling, the job with the least attained service gets priority.
    This is also known as Least Attained Service (LAS) scheduling.

    For FB, the conditional response time for a job of size x is:
        E[T|x] = x / (1 - rho)  for exponential service
        E[T|x] = integral formula for general service

    Args:
        lambda_val: Arrival rate
        mu_val: Service rate
        cs_val: Coefficient of variation of service time

    Returns:
        Tuple of (W, rho):
            W: Mean response time
            rho: System utilization (rhohat format)

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mg1_fb.m
        Kleinrock, "Queueing Systems, Volume II"
    """
    rho = lambda_val / mu_val

    if rho >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho:.4g} >= 1")

    # For exponential service (cs = 1), FB and FCFS have same mean response time
    # For general service, FB favors small jobs

    if abs(cs_val - 1.0) < 1e-10:
        # Exponential case
        W = 1.0 / (mu_val * (1 - rho))
    else:
        # General case - approximation
        # E[T] = E[S] / (1 - rho) * (1 + correction term based on cs)
        mean_S = 1.0 / mu_val
        var_S = (cs_val ** 2) / (mu_val ** 2)

        # Use weighted approximation
        W = mean_S / (1 - rho)

        # Correction for variability (FB benefits from low variability)
        if cs_val < 1:
            W *= (1 + cs_val ** 2) / 2
        else:
            W *= (1 + (cs_val ** 2 - 1) / (2 * (1 + rho)))

    Q = lambda_val * W
    rho_hat = Q / (1 + Q)

    return W, rho_hat


def qsys_mg1_lrpt(lambda_vec: np.ndarray, mu_vec: np.ndarray,
                  cs_vec: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Analyze M/G/1 queue with Longest Remaining Processing Time (LRPT).

    LRPT gives preemptive priority to jobs with the longest remaining
    processing time. This is the opposite of SRPT.

    Args:
        lambda_vec: Vector of arrival rates per class
        mu_vec: Vector of service rates per class
        cs_vec: Vector of coefficients of variation per class

    Returns:
        Tuple of (W, rho):
            W: Vector of mean response times per class
            rho: System utilization (rhohat format)

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mg1_lrpt.m
    """
    lambda_vec = np.asarray(lambda_vec, dtype=float).flatten()
    mu_vec = np.asarray(mu_vec, dtype=float).flatten()
    cs_vec = np.asarray(cs_vec, dtype=float).flatten()

    K = len(lambda_vec)
    mean_service = 1.0 / mu_vec

    # Sort by mean service time descending (longest first)
    sort_idx = np.argsort(-mean_service)
    unsort_idx = np.argsort(sort_idx)

    lambda_sorted = lambda_vec[sort_idx]
    mu_sorted = mu_vec[sort_idx]

    # Use preemptive priority formula with reversed ordering
    rho_i = lambda_sorted / mu_sorted
    rho_total = np.sum(rho_i)

    if rho_total >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_total:.4g} >= 1")

    W_sorted = _mg1_srpt_exp(lambda_sorted, mu_sorted)

    W = W_sorted[unsort_idx]

    Q = np.sum(lambda_vec * W)
    rho_hat = Q / (1 + Q)

    return W, rho_hat


def qsys_mg1_psjf(lambda_vec: np.ndarray, mu_vec: np.ndarray,
                  cs_vec: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Analyze M/G/1 queue with Preemptive Shortest Job First (PSJF).

    PSJF gives preemptive priority to the job with the smallest
    expected total service time (not remaining time).

    For exponential service, PSJF equals SRPT.

    Args:
        lambda_vec: Vector of arrival rates per class
        mu_vec: Vector of service rates per class
        cs_vec: Vector of coefficients of variation per class

    Returns:
        Tuple of (W, rho):
            W: Vector of mean response times per class
            rho: System utilization (rhohat format)

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mg1_psjf.m
    """
    # For PSJF with known job sizes, it's equivalent to SRPT for classification
    return qsys_mg1_srpt(lambda_vec, mu_vec, cs_vec)


def qsys_mg1_setf(lambda_vec: np.ndarray, mu_vec: np.ndarray,
                  cs_vec: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Analyze M/G/1 queue with Shortest Expected Time First (SETF).

    SETF is non-preemptive and prioritizes based on expected service time.

    Args:
        lambda_vec: Vector of arrival rates per class
        mu_vec: Vector of service rates per class
        cs_vec: Vector of coefficients of variation per class

    Returns:
        Tuple of (W, rho):
            W: Vector of mean response times per class
            rho: System utilization (rhohat format)

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_mg1_setf.m
    """
    # SETF with classes = non-preemptive priority by expected service time
    # Sort by mean service time and use HOL priority formula
    lambda_vec = np.asarray(lambda_vec, dtype=float).flatten()
    mu_vec = np.asarray(mu_vec, dtype=float).flatten()
    cs_vec = np.asarray(cs_vec, dtype=float).flatten()

    mean_service = 1.0 / mu_vec
    sort_idx = np.argsort(mean_service)
    unsort_idx = np.argsort(sort_idx)

    W_sorted, rho = qsys_mg1_prio(
        lambda_vec[sort_idx], mu_vec[sort_idx], cs_vec[sort_idx]
    )

    W = W_sorted[unsort_idx]

    return W, rho


__all__ = [
    'qsys_mg1_prio',
    'qsys_mg1_srpt',
    'qsys_mg1_fb',
    'qsys_mg1_lrpt',
    'qsys_mg1_psjf',
    'qsys_mg1_setf',
]
