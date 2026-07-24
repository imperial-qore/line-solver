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
from scipy.special import gammaln


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

    SRPT is a size-based policy: it always serves the job with the smallest
    remaining processing time, preempting whenever a shorter job arrives. The
    class-conditional mean response time follows the Schrage-Miller formula
    (Bansal-Harchol-Balter, SIGMETRICS 2001, Sec. 4, Eqs (1)-(3), after
    Schrage-Miller 1966). For a job of size x:

        E[T(x)] = E[W(x)] + E[R(x)]
        E[W(x)] = lambda*(m2(x) + x^2*(1-F(x))) / (2*(1-rho(x))^2)
        E[R(x)] = integral_0^x dt/(1-rho(t))

    with f the mixture job-size density, F its CDF, rho(x)=lambda*int_0^x t f(t)dt
    and m2(x)=int_0^x t^2 f(t)dt. The per-class mean is
    E[T_r]=int_0^inf E[T(x)] f_r(x) dx; since E[T(x)] depends only on the job
    size (SRPT is size-based) this is exact. Integrals use cumulative
    trapezoidal quadrature on a common grid. Each class is matched to its
    (mean=1/mu, scv=cs^2): exponential for cs=1, a two-phase balanced
    hyperexponential for cs>1, and a Tijms Erlang-(k-1)/Erlang-k mixture for
    cs<1. The fully exponential case reproduces the exact M/M/1/SRPT result.

    Matches MATLAB qsys_mg1_srpt.m.

    Args:
        lambda_vec: Vector of arrival rates per class
        mu_vec: Vector of service rates per class
        cs_vec: Vector of coefficients of variation per class

    Returns:
        Tuple of (W, rho):
            W: Vector of mean response times per class (original class order)
            rho: System load measure Q/(1+Q) with Q = sum(lambda.*W)
    """
    lambda_vec = np.asarray(lambda_vec, dtype=float).flatten()
    mu_vec = np.asarray(mu_vec, dtype=float).flatten()
    cs_vec = np.asarray(cs_vec, dtype=float).flatten()

    if not (len(lambda_vec) == len(mu_vec) == len(cs_vec)):
        raise ValueError("lambda, mu, and cs must have the same length")

    if np.any(lambda_vec <= 0) or np.any(mu_vec <= 0) or np.any(cs_vec < 0):
        raise ValueError("lambda and mu must be positive, cs must be non-negative")

    K = len(lambda_vec)
    lambda_total = np.sum(lambda_vec)
    p = lambda_vec / lambda_total

    rho_util = np.sum(lambda_vec / mu_vec)
    if rho_util >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_util:.4g} >= 1")

    # Per-class job-size representations and phase-rate bounds
    fits = [_srpt_fit(mu_vec[r], cs_vec[r]) for r in range(K)]
    rate_min = min(f['rate_min'] for f in fits)
    rate_max = max(f['rate_max'] for f in fits)

    # Integration grid: 40 e-foldings of the slowest phase, at least 200
    # points per rate ratio to resolve the fastest phase.
    xmax = 40.0 / rate_min
    N = int(min(2_000_000, max(20000, np.ceil(200.0 * rate_max / rate_min))))
    x = np.linspace(0.0, xmax, N + 1)
    dx = x[1] - x[0]

    fmix = np.zeros(N + 1)
    fbar = np.zeros(N + 1)
    for r in range(K):
        fmix += p[r] * _srpt_pdf(fits[r], x)
        fbar += p[r] * _srpt_tail(fits[r], x)

    def _cumtrap(y):
        return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * dx)])

    rho_x = lambda_total * _cumtrap(x * fmix)
    m2_x = _cumtrap(x * x * fmix)
    denom = np.maximum(1.0 - rho_x, 1e-12)

    wait = lambda_total * (m2_x + x * x * fbar) / (2.0 * denom ** 2)
    res = _cumtrap(1.0 / denom)
    ET = wait + res

    W = np.zeros(K)
    for r in range(K):
        integ = ET * _srpt_pdf(fits[r], x)
        W[r] = np.sum(0.5 * (integ[1:] + integ[:-1]) * dx)

    Q = np.sum(lambda_vec * W)
    rho_hat = Q / (1 + Q)

    return W, rho_hat


def _srpt_fit(mu, cs):
    """Match a job-size distribution to mean 1/mu and scv cs^2.

    Returns a dict with a type tag ('exp', 'h2', 'erlmix'), parameters, and
    the min/max phase rates used to size the integration grid.
    """
    mean_x = 1.0 / mu
    c2 = cs * cs
    if abs(c2 - 1.0) < 1e-9:
        return {'type': 'exp', 'rate': mu, 'rate_min': mu, 'rate_max': mu}
    if c2 > 1.0:
        pr = 0.5 * (1.0 + np.sqrt((c2 - 1.0) / (c2 + 1.0)))
        r1 = 2.0 * pr * mu
        r2 = 2.0 * (1.0 - pr) * mu
        return {'type': 'h2', 'p': pr, 'r1': r1, 'r2': r2,
                'rate_min': min(r1, r2), 'rate_max': max(r1, r2)}
    # c2 < 1: Tijms mixture of Erlang-(k-1) and Erlang-k with common rate
    k = int(np.ceil(1.0 / c2))
    pr = (1.0 / (1.0 + c2)) * (k * c2 - np.sqrt(k * (1.0 + c2) - k * k * c2))
    rate = (k - pr) / mean_x
    return {'type': 'erlmix', 'k': k, 'p': pr, 'rate': rate,
            'rate_min': rate, 'rate_max': rate}


def _srpt_pdf(fit, x):
    """Job-size probability density on the grid x."""
    t = fit['type']
    if t == 'exp':
        return fit['rate'] * np.exp(-fit['rate'] * x)
    if t == 'h2':
        return (fit['p'] * fit['r1'] * np.exp(-fit['r1'] * x)
                + (1.0 - fit['p']) * fit['r2'] * np.exp(-fit['r2'] * x))
    return (fit['p'] * _erlang_pdf(fit['k'] - 1, fit['rate'], x)
            + (1.0 - fit['p']) * _erlang_pdf(fit['k'], fit['rate'], x))


def _srpt_tail(fit, x):
    """Complementary CDF P(X > x) on the grid x."""
    t = fit['type']
    if t == 'exp':
        return np.exp(-fit['rate'] * x)
    if t == 'h2':
        return (fit['p'] * np.exp(-fit['r1'] * x)
                + (1.0 - fit['p']) * np.exp(-fit['r2'] * x))
    return (fit['p'] * _erlang_tail(fit['k'] - 1, fit['rate'], x)
            + (1.0 - fit['p']) * _erlang_tail(fit['k'], fit['rate'], x))


def _erlang_pdf(n, rate, x):
    """Erlang-n density in log space: f(x) = rate * Poisson(n-1; rate*x).
    n=0 is a point mass at 0 (density 0 for x>0)."""
    if n <= 0:
        return np.zeros_like(x)
    t = rate * x
    m = n - 1
    with np.errstate(divide='ignore', invalid='ignore'):
        logp = m * np.log(t) - t - gammaln(m + 1)
    logp = np.where(t > 0, logp, (0.0 if m == 0 else -np.inf))
    return rate * np.exp(logp)


def _erlang_tail(n, rate, x):
    """Erlang-n complementary CDF P(X>x) = sum_{j=0}^{n-1} Poisson(j; rate*x),
    computed in log space. n=0 tail is 0 for x>0."""
    if n <= 0:
        return np.zeros_like(x)
    t = rate * x
    y = np.zeros_like(x)
    for j in range(n):
        with np.errstate(divide='ignore', invalid='ignore'):
            logp = j * np.log(t) - t - gammaln(j + 1)
        logp = np.where(t > 0, logp, (0.0 if j == 0 else -np.inf))
        y = y + np.exp(logp)
    return y


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


def qsys_mm1_dps(lambda_vec: np.ndarray, mu_vec: np.ndarray,
                 w_vec: np.ndarray, tol: float = 1e-10,
                 max_cutoff: int = 2048) -> Tuple[np.ndarray, float]:
    """
    Numerically exact M/M/1 Discriminatory Processor Sharing (DPS) queue.

    Solves the multiclass DPS continuous-time Markov chain on the per-class
    population vector (n_1..n_K): arrivals lambda_k, class-k service completion
    rate mu_k * n_k * w_k / sum_j n_j * w_j. The state space is truncated at a
    total population level chosen from the geometric tail bound (the total-count
    process is stochastically dominated by the M/M/1 with rate min_k mu_k), and
    the truncation level is doubled until the mean queue lengths are stable to
    the requested tolerance -- so the result is exact to solver precision and
    conserves the M/M/1 total for equal service rates by construction.

    Args:
        lambda_vec: Per-class Poisson arrival rates (K,)
        mu_vec: Per-class exponential service rates (K,)
        w_vec: Per-class DPS weights (K,), positive
        tol: Convergence tolerance on the per-class mean counts
        max_cutoff: Hard bound on the total-population truncation level

    Returns:
        Tuple of (T, rho): per-class mean response times (K,) via Little's law,
        and the total utilization sum_k lambda_k/mu_k.
    """
    from itertools import product as _prod
    import scipy.sparse as _sp
    import scipy.sparse.linalg as _spla

    lam = np.asarray(lambda_vec, dtype=float).ravel()
    mu = np.asarray(mu_vec, dtype=float).ravel()
    w = np.asarray(w_vec, dtype=float).ravel()
    K = len(lam)
    if not (len(mu) == K and len(w) == K):
        raise ValueError("lambda, mu, w must have the same length")
    if np.any(lam <= 0) or np.any(mu <= 0) or np.any(w <= 0):
        raise ValueError("lambda, mu, w must all be positive")
    rho = float(np.sum(lam / mu))
    if rho >= 1:
        raise ValueError("System is unstable: utilization rho = %.4g >= 1" % rho)

    def _solve(N):
        # enumerate states with total population <= N
        states = [s for s in _prod(range(N + 1), repeat=K) if sum(s) <= N]
        idx = {s: i for i, s in enumerate(states)}
        n = len(states)
        rows, cols, vals = [], [], []
        for s in states:
            i = idx[s]
            tot = sum(s)
            out = 0.0
            # arrivals
            if tot < N:
                for k in range(K):
                    s2 = list(s); s2[k] += 1
                    j = idx[tuple(s2)]
                    rows.append(i); cols.append(j); vals.append(lam[k])
                    out += lam[k]
            # departures (DPS capacity split)
            if tot > 0:
                denom = sum(s[k] * w[k] for k in range(K))
                for k in range(K):
                    if s[k] > 0:
                        r = mu[k] * s[k] * w[k] / denom
                        s2 = list(s); s2[k] -= 1
                        j = idx[tuple(s2)]
                        rows.append(i); cols.append(j); vals.append(r)
                        out += r
            rows.append(i); cols.append(i); vals.append(-out)
        Q = _sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
        # stationary distribution: solve pi Q = 0, sum pi = 1
        A = Q.T.tolil()
        A[0, :] = 1.0
        b = np.zeros(n); b[0] = 1.0
        pi = _spla.spsolve(A.tocsr(), b)
        EN = np.zeros(K)
        for s in states:
            p = pi[idx[s]]
            for k in range(K):
                EN[k] += p * s[k]
        return EN

    # adaptive truncation: start from a tail-bound estimate, double until stable
    N = max(16, int(np.ceil(np.log(tol) / np.log(rho))) if rho > 0 else 16)
    N = min(N, max_cutoff)
    EN_prev = _solve(N)
    while N < max_cutoff:
        N2 = min(2 * N, max_cutoff)
        EN = _solve(N2)
        if np.max(np.abs(EN - EN_prev)) < tol:
            EN_prev = EN
            break
        EN_prev, N = EN, N2
        if N2 == max_cutoff:
            break
    T = EN_prev / lam
    return T, rho


__all__ = [
    'qsys_mg1_prio',
    'qsys_mg1_srpt',
    'qsys_mg1_fb',
    'qsys_mg1_lrpt',
    'qsys_mg1_psjf',
    'qsys_mg1_setf',
    'qsys_mm1_dps',
]
