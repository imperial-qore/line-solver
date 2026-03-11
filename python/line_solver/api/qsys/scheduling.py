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
    Wierman and Harchol-Balter, SIGMETRICS 2003
    Kleinrock, "Queueing Systems, Volume I: Theory", 1975
"""

import numpy as np
from typing import Tuple
from scipy.integrate import quad


def qsys_mg1_prio(lambda_vec: np.ndarray, mu_vec: np.ndarray,
                  cs_vec: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Analyze M/G/1 queue with non-preemptive (Head-of-Line) priorities.

    Matches MATLAB qsys_mg1_prio.m exactly.

    Args:
        lambda_vec: Vector of arrival rates per priority class (class 1 = highest)
        mu_vec: Vector of service rates per priority class
        cs_vec: Vector of coefficients of variation per priority class

    Returns:
        Tuple of (W, rho):
            W: Vector of mean response times per priority class
            rho: System utilization (rhohat = Q/(1+Q) format)
    """
    lambda_vec = np.asarray(lambda_vec, dtype=float).flatten()
    mu_vec = np.asarray(mu_vec, dtype=float).flatten()
    cs_vec = np.asarray(cs_vec, dtype=float).flatten()

    if not (len(lambda_vec) == len(mu_vec) == len(cs_vec)):
        raise ValueError("lambda, mu, and cs must have the same length")

    if np.any(lambda_vec <= 0) or np.any(mu_vec <= 0) or np.any(cs_vec <= 0):
        raise ValueError("lambda, mu, and cs must all be positive")

    rho_i = lambda_vec / mu_vec
    rho_total = np.sum(rho_i)

    if rho_total >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_total:.4g} >= 1")

    # B_0 = sum_i lambda_i * (1 + cs_i^2) / mu_i^2 / 2
    B_0 = np.sum(lambda_vec * (1 + cs_vec ** 2) / (mu_vec ** 2)) / 2

    K = len(lambda_vec)
    W_q = np.zeros(K)

    for k in range(K):
        rho_prev = np.sum(rho_i[:k]) if k > 0 else 0.0
        rho_curr = np.sum(rho_i[:k + 1])
        W_q[k] = B_0 / ((1 - rho_prev) * (1 - rho_curr))

    W = W_q + 1.0 / mu_vec

    Q = np.sum(lambda_vec * W)
    rho_hat = Q / (1 + Q)

    return W, rho_hat


def qsys_mg1_srpt(lambda_vec: np.ndarray, mu_vec: np.ndarray,
                  cs_vec: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Analyze M/G/1 queue with Shortest Remaining Processing Time (SRPT).

    Matches MATLAB qsys_mg1_srpt.m exactly:
    - Exponential case: preemptive priority with cumulative residual E_R_k
    - General case: Schrage-Miller formula with numerical integration

    Args:
        lambda_vec: Vector of arrival rates per class
        mu_vec: Vector of service rates per class
        cs_vec: Vector of coefficients of variation per class

    Returns:
        Tuple of (W, rho):
            W: Vector of mean response times per class
            rho: System utilization (rhohat format)
    """
    lambda_vec = np.asarray(lambda_vec, dtype=float).flatten()
    mu_vec = np.asarray(mu_vec, dtype=float).flatten()
    cs_vec = np.asarray(cs_vec, dtype=float).flatten()

    if not (len(lambda_vec) == len(mu_vec) == len(cs_vec)):
        raise ValueError("lambda, mu, and cs must have the same length")

    if np.any(lambda_vec <= 0) or np.any(mu_vec <= 0) or np.any(cs_vec < 0):
        raise ValueError("lambda and mu must be positive, cs must be non-negative")

    K = len(lambda_vec)
    mean_service = 1.0 / mu_vec

    # Sort by mean service time ascending
    sort_idx = np.argsort(mean_service)
    unsort_idx = np.argsort(sort_idx)

    lambda_sorted = lambda_vec[sort_idx]
    mu_sorted = mu_vec[sort_idx]
    cs_sorted = cs_vec[sort_idx]

    rho_i = lambda_sorted / mu_sorted
    rho_total = np.sum(rho_i)

    if rho_total >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_total:.4g} >= 1")

    is_exponential = np.all(np.abs(cs_sorted - 1) < 1e-10)

    if is_exponential or K == 1:
        W_sorted = _srpt_exp(lambda_sorted, mu_sorted)
    else:
        W_sorted = _srpt_general(lambda_sorted, mu_sorted, cs_sorted)

    W = W_sorted[unsort_idx]

    Q = np.sum(lambda_vec * W)
    rho_hat = Q / (1 + Q)

    return W, rho_hat


def _srpt_exp(lambda_arr, mu_arr):
    """SRPT for exponential service - preemptive priority formula.
    Matches MATLAB qsys_mg1_srpt_exp."""
    K = len(lambda_arr)
    rho_i = lambda_arr / mu_arr

    W = np.zeros(K)
    for k in range(K):
        rho_prev = np.sum(rho_i[:k]) if k > 0 else 0.0
        rho_curr = np.sum(rho_i[:k + 1])

        # Cumulative residual service from classes 1 to k
        E_R_k = np.sum(lambda_arr[:k + 1] / (mu_arr[:k + 1] ** 2))

        W_q = E_R_k / ((1 - rho_prev) * (1 - rho_curr))
        W[k] = W_q + 1.0 / mu_arr[k]

    return W


def _srpt_general(lambda_arr, mu_arr, cs_arr):
    """SRPT for general service - Schrage-Miller formula.
    Matches MATLAB qsys_mg1_srpt_general."""
    K = len(lambda_arr)
    lambda_total = np.sum(lambda_arr)
    p = lambda_arr / lambda_total

    W = np.zeros(K)
    for k in range(K):
        x = 1.0 / mu_arr[k]

        # Truncated load from classes with smaller mean service times
        smaller_idx = (1.0 / mu_arr) <= x
        rho_x = np.sum(lambda_arr[smaller_idx] / mu_arr[smaller_idx])

        if rho_x >= 1:
            W[k] = np.inf
            continue

        # Numerator
        numerator = 0.0
        for i in range(K):
            if 1.0 / mu_arr[i] <= x:
                numerator += lambda_arr[i] / (mu_arr[i] ** 2)

        # Second integral: integral_0^x dt/(1-rho(t))
        n_steps = 1000
        dt = x / n_steps
        second_term = 0.0
        for step in range(1, n_steps + 1):
            t = (step - 0.5) * dt
            rho_t = 0.0
            for i in range(K):
                rho_t += lambda_arr[i] / mu_arr[i] * (1 - np.exp(-mu_arr[i] * t) * (1 + mu_arr[i] * t))
            if rho_t < 1:
                second_term += dt / (1 - rho_t)

        W[k] = numerator / (1 - rho_x) ** 2 + second_term

    return W


def qsys_mg1_fb(lambda_vec: np.ndarray, mu_vec: np.ndarray,
                cs_vec: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Analyze M/G/1 queue with Foreground-Background (FB/LAS) scheduling.

    Matches MATLAB qsys_mg1_fb.m exactly:
    - Exponential case: numerical integration of E[T(x)] * f_k(x)
    - General case: class-based approximation

    Args:
        lambda_vec: Vector of arrival rates per class
        mu_vec: Vector of service rates per class
        cs_vec: Vector of coefficients of variation per class

    Returns:
        Tuple of (W, rho):
            W: Vector of mean response times per class
            rho: System utilization (rhohat format)
    """
    lambda_vec = np.asarray(lambda_vec, dtype=float).flatten()
    mu_vec = np.asarray(mu_vec, dtype=float).flatten()
    cs_vec = np.asarray(cs_vec, dtype=float).flatten()

    if not (len(lambda_vec) == len(mu_vec) == len(cs_vec)):
        raise ValueError("lambda, mu, and cs must have the same length")

    if np.any(lambda_vec <= 0) or np.any(mu_vec <= 0) or np.any(cs_vec < 0):
        raise ValueError("lambda and mu must be positive, cs must be non-negative")

    K = len(lambda_vec)
    rho_total = np.sum(lambda_vec / mu_vec)

    if rho_total >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_total:.4g} >= 1")

    if np.all(np.abs(cs_vec - 1) < 1e-6):
        W = _fb_exp(lambda_vec, mu_vec)
    else:
        W = _fb_general(lambda_vec, mu_vec, cs_vec)

    Q = np.sum(lambda_vec * W)
    rho_hat = Q / (1 + Q)

    return W, rho_hat


def _compute_fb_response(x, lambda_arr, mu_arr, p, lambda_total):
    """Compute FB/LAS response time E[T(x)] for a job of size x.
    Matches MATLAB compute_fb_response."""
    K = len(lambda_arr)

    # rho_x = lambda * integral_0^x F_bar(t) dt
    rho_x = 0.0
    for i in range(K):
        mu_i = mu_arr[i]
        int_Fbar = (1 - np.exp(-mu_i * x)) / mu_i
        rho_x += p[i] * lambda_total * int_Fbar

    # numerator = lambda * integral_0^x t * F_bar(t) dt
    numerator = 0.0
    for i in range(K):
        mu_i = mu_arr[i]
        int_tFbar = (1 - np.exp(-mu_i * x) * (1 + mu_i * x)) / mu_i**2
        numerator += p[i] * lambda_total * int_tFbar

    if rho_x >= 1:
        return np.inf

    return numerator / (1 - rho_x)**2 + x / (1 - rho_x)


def _fb_exp(lambda_arr, mu_arr):
    """FB for exponential service - numerical integration.
    Matches MATLAB qsys_mg1_fb_exp."""
    K = len(lambda_arr)
    lambda_total = np.sum(lambda_arr)
    p = lambda_arr / lambda_total

    W = np.zeros(K)
    for k in range(K):
        mu_k = mu_arr[k]
        x_max = 20 / mu_k

        def integrand(x):
            T = _compute_fb_response(x, lambda_arr, mu_arr, p, lambda_total)
            f_k = mu_k * np.exp(-mu_k * x)
            return T * f_k

        W[k], _ = quad(integrand, 0, x_max, limit=200, epsrel=1e-8, epsabs=1e-10)

    return W


def _fb_general(lambda_arr, mu_arr, cs_arr):
    """FB for general service - class-based approximation.
    Matches MATLAB qsys_mg1_fb_general."""
    K = len(lambda_arr)
    W = np.zeros(K)

    for k in range(K):
        x = 1.0 / mu_arr[k]

        rho_x = 0.0
        for i in range(K):
            if abs(cs_arr[i] - 1) < 1e-6:
                integral_Fbar = (1 - np.exp(-mu_arr[i] * x)) / mu_arr[i]
            else:
                integral_Fbar = min(x, 1.0 / mu_arr[i])
            rho_x += lambda_arr[i] * integral_Fbar

        numerator = 0.0
        for i in range(K):
            if abs(cs_arr[i] - 1) < 1e-6:
                mu_i = mu_arr[i]
                integral_tFbar = (1 - np.exp(-mu_i * x) * (1 + mu_i * x)) / mu_i**2
            else:
                integral_tFbar = min(x**2 / 2, 1.0 / mu_arr[i]**2)
            numerator += lambda_arr[i] * integral_tFbar

        if rho_x >= 1:
            W[k] = np.inf
        else:
            W[k] = numerator / (1 - rho_x)**2 + x / (1 - rho_x)

    return W


def qsys_mg1_lrpt(lambda_vec: np.ndarray, mu_vec: np.ndarray,
                  cs_vec: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Analyze M/G/1 queue with Longest Remaining Processing Time (LRPT).

    Matches MATLAB qsys_mg1_lrpt.m exactly:
    - Exponential case: numerical integration of E[T(x)] * f_k(x)
    - General case: preemptive priority with descending service time ordering

    Args:
        lambda_vec: Vector of arrival rates per class
        mu_vec: Vector of service rates per class
        cs_vec: Vector of coefficients of variation per class

    Returns:
        Tuple of (W, rho):
            W: Vector of mean response times per class
            rho: System utilization (rhohat format)
    """
    lambda_vec = np.asarray(lambda_vec, dtype=float).flatten()
    mu_vec = np.asarray(mu_vec, dtype=float).flatten()
    cs_vec = np.asarray(cs_vec, dtype=float).flatten()

    if not (len(lambda_vec) == len(mu_vec) == len(cs_vec)):
        raise ValueError("lambda, mu, and cs must have the same length")

    if np.any(lambda_vec <= 0) or np.any(mu_vec <= 0) or np.any(cs_vec < 0):
        raise ValueError("lambda and mu must be positive, cs must be non-negative")

    K = len(lambda_vec)
    rho_total = np.sum(lambda_vec / mu_vec)

    if rho_total >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_total:.4g} >= 1")

    if np.all(np.abs(cs_vec - 1) < 1e-6):
        W = _lrpt_exp(lambda_vec, mu_vec)
    else:
        W = _lrpt_general(lambda_vec, mu_vec, cs_vec)

    Q = np.sum(lambda_vec * W)
    rho_hat = Q / (1 + Q)

    return W, rho_hat


def _lrpt_exp(lambda_arr, mu_arr):
    """LRPT for exponential service - numerical integration.
    Matches MATLAB qsys_mg1_lrpt_exp."""
    K = len(lambda_arr)
    lambda_total = np.sum(lambda_arr)
    rho_total = np.sum(lambda_arr / mu_arr)
    p = lambda_arr / lambda_total

    # E[X^2] for the mixture = sum_i p_i * 2/mu_i^2
    E_X2 = np.sum(p * 2.0 / (mu_arr**2))

    W = np.zeros(K)
    for k in range(K):
        mu_k = mu_arr[k]
        x_max = 20 / mu_k

        def T_of_x(x):
            return x / (1 - rho_total) + lambda_total * E_X2 / (2 * (1 - rho_total)**2)

        def integrand(x):
            f_k = mu_k * np.exp(-mu_k * x)
            return T_of_x(x) * f_k

        W[k], _ = quad(integrand, 0, x_max, limit=200, epsrel=1e-8, epsabs=1e-10)

    return W


def _lrpt_general(lambda_arr, mu_arr, cs_arr):
    """LRPT for general service - preemptive priority with descending ordering.
    Matches MATLAB qsys_mg1_lrpt_general."""
    K = len(lambda_arr)
    mean_service = 1.0 / mu_arr

    # Sort by mean service time DESCENDING for LRPT priority
    sort_idx = np.argsort(-mean_service)
    unsort_idx = np.argsort(sort_idx)

    lambda_sorted = lambda_arr[sort_idx]
    mu_sorted = mu_arr[sort_idx]
    rho_i = lambda_sorted / mu_sorted

    W_sorted = np.zeros(K)
    for k in range(K):
        rho_prev = np.sum(rho_i[:k]) if k > 0 else 0.0
        rho_curr = np.sum(rho_i[:k + 1])
        E_R_k = np.sum(lambda_sorted[:k + 1] / (mu_sorted[:k + 1]**2))
        W_q = E_R_k / ((1 - rho_prev) * (1 - rho_curr))
        W_sorted[k] = W_q + 1.0 / mu_sorted[k]

    return W_sorted[unsort_idx]


def qsys_mg1_psjf(lambda_vec: np.ndarray, mu_vec: np.ndarray,
                  cs_vec: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Analyze M/G/1 queue with Preemptive Shortest Job First (PSJF).

    Matches MATLAB qsys_mg1_psjf.m exactly:
    - Exponential case: numerical integration of E[T(x)] * f_k(x)
    - General case: class-based truncated moment formula

    Args:
        lambda_vec: Vector of arrival rates per class
        mu_vec: Vector of service rates per class
        cs_vec: Vector of coefficients of variation per class

    Returns:
        Tuple of (W, rho):
            W: Vector of mean response times per class
            rho: System utilization (rhohat format)
    """
    lambda_vec = np.asarray(lambda_vec, dtype=float).flatten()
    mu_vec = np.asarray(mu_vec, dtype=float).flatten()
    cs_vec = np.asarray(cs_vec, dtype=float).flatten()

    if not (len(lambda_vec) == len(mu_vec) == len(cs_vec)):
        raise ValueError("lambda, mu, and cs must have the same length")

    if np.any(lambda_vec <= 0) or np.any(mu_vec <= 0) or np.any(cs_vec < 0):
        raise ValueError("lambda and mu must be positive, cs must be non-negative")

    K = len(lambda_vec)
    rho_total = np.sum(lambda_vec / mu_vec)

    if rho_total >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_total:.4g} >= 1")

    if np.all(np.abs(cs_vec - 1) < 1e-6):
        W = _psjf_exp(lambda_vec, mu_vec)
    else:
        W = _psjf_general(lambda_vec, mu_vec, cs_vec)

    Q = np.sum(lambda_vec * W)
    rho_hat = Q / (1 + Q)

    return W, rho_hat


def _compute_psjf_response(x, lambda_arr, mu_arr, p, lambda_total):
    """Compute PSJF response time E[T(x)] for a job of size x.
    Matches MATLAB compute_psjf_response."""
    K = len(lambda_arr)

    # Truncated first moment: integral_0^x t * f(t) dt
    m1_x = 0.0
    for i in range(K):
        mu_i = mu_arr[i]
        int_t = 1.0 / mu_i - (1.0 / mu_i + x) * np.exp(-mu_i * x)
        m1_x += p[i] * int_t

    rho_x = lambda_total * m1_x

    # Truncated second moment: integral_0^x t^2 * f(t) dt
    m2_x = 0.0
    for i in range(K):
        mu_i = mu_arr[i]
        int_t2 = 2.0 / mu_i**2 - (2.0 / mu_i**2 + 2.0 * x / mu_i + x**2) * np.exp(-mu_i * x)
        m2_x += p[i] * int_t2

    m2_x_scaled = lambda_total * m2_x

    if rho_x >= 1:
        return np.inf

    return x / (1 - rho_x) + m2_x_scaled / (2 * (1 - rho_x)**2)


def _psjf_exp(lambda_arr, mu_arr):
    """PSJF for exponential service - numerical integration.
    Matches MATLAB qsys_mg1_psjf_exp."""
    K = len(lambda_arr)
    lambda_total = np.sum(lambda_arr)
    p = lambda_arr / lambda_total

    W = np.zeros(K)
    for k in range(K):
        mu_k = mu_arr[k]
        x_max = 20 / mu_k

        def integrand(x):
            T = _compute_psjf_response(x, lambda_arr, mu_arr, p, lambda_total)
            f_k = mu_k * np.exp(-mu_k * x)
            return T * f_k

        W[k], _ = quad(integrand, 0, x_max, limit=200, epsrel=1e-8, epsabs=1e-10)

    return W


def _psjf_general(lambda_arr, mu_arr, cs_arr):
    """PSJF for general service - class-based formula.
    Matches MATLAB qsys_mg1_psjf_general."""
    K = len(lambda_arr)
    mean_service = 1.0 / mu_arr

    sort_idx = np.argsort(mean_service)
    unsort_idx = np.argsort(sort_idx)

    lambda_sorted = lambda_arr[sort_idx]
    mu_sorted = mu_arr[sort_idx]
    cs_sorted = cs_arr[sort_idx]
    rho_i = lambda_sorted / mu_sorted

    W_sorted = np.zeros(K)
    for k in range(K):
        x = 1.0 / mu_sorted[k]
        rho_x = np.sum(rho_i[:k + 1])

        m2_x = 0.0
        for i in range(k + 1):
            E_S2_i = (1 + cs_sorted[i]**2) / mu_sorted[i]**2
            m2_x += lambda_sorted[i] * E_S2_i

        if rho_x >= 1:
            W_sorted[k] = np.inf
        else:
            waiting_term = m2_x / (2 * (1 - rho_x)**2)
            service_term = x / (1 - rho_x)
            W_sorted[k] = waiting_term + service_term

    return W_sorted[unsort_idx]


def qsys_mg1_setf(lambda_vec: np.ndarray, mu_vec: np.ndarray,
                  cs_vec: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Analyze M/G/1 queue with Shortest Expected Time First (SETF).

    Matches MATLAB qsys_mg1_setf.m exactly:
    SETF = FB/LAS + residual service time penalty (non-preemptive).

    Args:
        lambda_vec: Vector of arrival rates per class
        mu_vec: Vector of service rates per class
        cs_vec: Vector of coefficients of variation per class

    Returns:
        Tuple of (W, rho):
            W: Vector of mean response times per class
            rho: System utilization (rhohat format)
    """
    lambda_vec = np.asarray(lambda_vec, dtype=float).flatten()
    mu_vec = np.asarray(mu_vec, dtype=float).flatten()
    cs_vec = np.asarray(cs_vec, dtype=float).flatten()

    if not (len(lambda_vec) == len(mu_vec) == len(cs_vec)):
        raise ValueError("lambda, mu, and cs must have the same length")

    if np.any(lambda_vec <= 0) or np.any(mu_vec <= 0) or np.any(cs_vec < 0):
        raise ValueError("lambda and mu must be positive, cs must be non-negative")

    K = len(lambda_vec)
    rho_i = lambda_vec / mu_vec
    rho_total = np.sum(rho_i)

    if rho_total >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_total:.4g} >= 1")

    # Mean residual service time for the mixture distribution
    lambda_total = np.sum(lambda_vec)
    E_R = 0.0
    for i in range(K):
        p_i = lambda_vec[i] / lambda_total
        E_S_i = 1.0 / mu_vec[i]
        E_S2_i = (1 + cs_vec[i]**2) / mu_vec[i]**2
        E_R += p_i * E_S2_i / (2 * E_S_i)

    W = np.zeros(K)
    for k in range(K):
        x = 1.0 / mu_vec[k]

        rho_x = 0.0
        for i in range(K):
            if abs(cs_vec[i] - 1) < 1e-6:
                integral_Fbar = (1 - np.exp(-mu_vec[i] * x)) / mu_vec[i]
            else:
                integral_Fbar = min(x, 1.0 / mu_vec[i])
            rho_x += lambda_vec[i] * integral_Fbar

        numerator = 0.0
        for i in range(K):
            if abs(cs_vec[i] - 1) < 1e-6:
                mu_i = mu_vec[i]
                integral_tFbar = (1 - np.exp(-mu_i * x) * (1 + mu_i * x)) / mu_i**2
            else:
                integral_tFbar = min(x**2 / 2, 1.0 / mu_vec[i]**2)
            numerator += lambda_vec[i] * integral_tFbar

        if rho_x >= 1:
            W[k] = np.inf
        else:
            fb_waiting_term = numerator / (1 - rho_x)**2
            fb_service_term = x / (1 - rho_x)
            np_penalty = E_R / (1 - rho_x)
            W[k] = fb_waiting_term + fb_service_term + np_penalty

    Q = np.sum(lambda_vec * W)
    rho_hat = Q / (1 + Q)

    return W, rho_hat


__all__ = [
    'qsys_mg1_prio',
    'qsys_mg1_srpt',
    'qsys_mg1_fb',
    'qsys_mg1_lrpt',
    'qsys_mg1_psjf',
    'qsys_mg1_setf',
]
