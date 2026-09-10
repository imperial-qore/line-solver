"""
Size-based preemptive scheduling policies for M/G/1 queues.

Native Python implementations for SRPT, PSJF, FB/LAS, and LRPT scheduling.
These implementations follow the analytical formulas from:

    A. Wierman and M. Harchol-Balter, "Classifying scheduling policies with
    respect to unfairness in an M/GI/1", SIGMETRICS 2003.

Key algorithms:
    qsys_mg1_srpt: Shortest Remaining Processing Time
    qsys_mg1_psjf: Preemptive Shortest Job First
    qsys_mg1_fb: Feedback / Least Attained Service (LAS)
    qsys_mg1_lrpt: Longest Remaining Processing Time
"""

import numpy as np
from typing import Tuple, Union, List

def qsys_mg1_srpt(
    lambda_vals: Union[np.ndarray, List[float]],
    mu_vals: Union[np.ndarray, List[float]],
    cs_vals: Union[np.ndarray, List[float]]
) -> Tuple[np.ndarray, float]:
    """Compute mean response time for M/G/1/SRPT queue.

    Thin wrapper around the canonical Schrage-Miller implementation in
    :mod:`line_solver.api.qsys.scheduling` (the single source of truth exported
    by the package). SRPT always serves the job with the smallest remaining
    processing time; the per-class mean response time is obtained by
    integrating the Schrage-Miller conditional response time E[T(x)] against
    each class job-size density. See scheduling.qsys_mg1_srpt for details.

    Args:
        lambda_vals: Array of arrival rates per class
        mu_vals: Array of service rates per class
        cs_vals: Array of coefficients of variation per class (cs=1 for exponential)

    Returns:
        Tuple[np.ndarray, float]: (W, rho) with W the per-class mean response
            times (original class order) and rho = Q/(1+Q), Q = sum(lambda.*W).
    """
    from .scheduling import qsys_mg1_srpt as _canonical
    return _canonical(lambda_vals, mu_vals, cs_vals)

def qsys_mg1_psjf(
    lambda_vals: Union[np.ndarray, List[float]],
    mu_vals: Union[np.ndarray, List[float]],
    cs_vals: Union[np.ndarray, List[float]]
) -> Tuple[np.ndarray, float]:
    """
    Compute mean response time for M/G/1/PSJF queue.

    Under PSJF (Preemptive Shortest Job First), priority is based on a job's
    original size (not remaining size). Jobs with smaller original sizes
    always preempt jobs with larger sizes.

    Classification (Wierman-Harchol-Balter 2003):
        PSJF is "Always Unfair" - some job size is treated unfairly under
        all loads and all service distributions.

    For PSJF, the mean response time for a job of size x is:
        E[T(x)]^PSJF = (lambda * integral_0^x t^2*f(t)dt) / (2*(1-rho(x))^2)
                       + x / (1 - rho(x))

    Args:
        lambda_vals: Array of arrival rates per class
        mu_vals: Array of service rates per class
        cs_vals: Array of coefficients of variation per class

    Returns:
        Tuple[np.ndarray, float]: (W, rho) where W is the vector of mean
            response times per class, and rho is the modified utilization.

    References:
        A. Wierman and M. Harchol-Balter, SIGMETRICS 2003, Section 3.2.

    Example:
        >>> W, rho = qsys_mg1_psjf([0.3, 0.2], [1.0, 0.5], [1.0, 1.0])
    """
    # Convert to numpy arrays
    lambda_arr = np.asarray(lambda_vals, dtype=float).flatten()
    mu_arr = np.asarray(mu_vals, dtype=float).flatten()
    cs_arr = np.asarray(cs_vals, dtype=float).flatten()

    # Validate input lengths
    if not (len(lambda_arr) == len(mu_arr) == len(cs_arr)):
        raise ValueError("lambda, mu, and cs must have the same length")

    # Validate positive values
    if np.any(lambda_arr <= 0) or np.any(mu_arr <= 0) or np.any(cs_arr < 0):
        raise ValueError("lambda and mu must be positive, cs must be non-negative")

    K = len(lambda_arr)

    # Compute mean service times per class
    mean_service = 1.0 / mu_arr

    # Sort classes by mean service time (ascending) for priority ordering
    sort_idx = np.argsort(mean_service)

    # Reorder parameters according to sorted service times
    lambda_sorted = lambda_arr[sort_idx]
    mu_sorted = mu_arr[sort_idx]
    cs_sorted = cs_arr[sort_idx]

    # Compute per-class utilizations (sorted order)
    rho_i = lambda_sorted / mu_sorted

    # Overall utilization
    rho_total = np.sum(rho_i)

    # Stability check
    if rho_total >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_total} >= 1")

    # Compute response times using PSJF formula
    W_sorted = np.zeros(K)

    for k in range(K):
        # Mean service time for this class
        x = 1.0 / mu_sorted[k]

        # Truncated load: rho(x) = sum of loads from smaller jobs
        rho_x = np.sum(rho_i[:k+1])

        # Truncated second moment: lambda * integral_0^x t^2*f(t)dt
        m2_x = 0.0
        for i in range(k + 1):
            # E[S_i^2] = (1 + cs_i^2) / mu_i^2
            E_S2_i = (1 + cs_sorted[i] ** 2) / (mu_sorted[i] ** 2)
            m2_x += lambda_sorted[i] * E_S2_i

        # PSJF formula:
        # E[T(x)] = (lambda * m2(x)) / (2*(1-rho(x))^2) + x / (1-rho(x))
        if rho_x >= 1:
            W_sorted[k] = np.inf
        else:
            waiting_term = m2_x / (2 * (1 - rho_x) ** 2)
            service_term = x / (1 - rho_x)
            W_sorted[k] = waiting_term + service_term

    # Restore original class ordering
    unsort_idx = np.argsort(sort_idx)
    W = W_sorted[unsort_idx]

    # Compute rhohat = Q/(1+Q) to match qsys convention
    Q = np.sum(lambda_arr * W)
    rho_out = Q / (1 + Q)

    return W, rho_out

def qsys_mg1_fb(
    lambda_vals: Union[np.ndarray, List[float]],
    mu_vals: Union[np.ndarray, List[float]],
    cs_vals: Union[np.ndarray, List[float]]
) -> Tuple[np.ndarray, float]:
    """
    Compute mean response time for M/G/1/FB (Feedback/LAS) queue.

    Under FB/LAS (Feedback / Least Attained Service), the job with the
    least attained service (smallest age) receives priority. This is an
    age-based policy where priority depends on how much service a job
    has received, not its original or remaining size.

    Also known as Least Attained Service (LAS) or Shortest Elapsed Time (SET).

    Classification (Wierman-Harchol-Balter 2003):
        FB is "Always Unfair" - some job size is treated unfairly under
        all loads and all service distributions. However, FB approximates
        SRPT for heavy-tailed distributions and is practical since job
        sizes need not be known in advance.

    For FB, the mean response time for a job of size x is:
        E[T(x)]^FB = (lambda * integral_0^x t*F_bar(t)dt) / (1-rho_x)^2
                     + x / (1 - rho_x)

    Args:
        lambda_vals: Array of arrival rates per class
        mu_vals: Array of service rates per class
        cs_vals: Array of coefficients of variation per class

    Returns:
        Tuple[np.ndarray, float]: (W, rho) where W is the vector of mean
            response times per class, and rho is the modified utilization.

    References:
        A. Wierman and M. Harchol-Balter, SIGMETRICS 2003, Section 3.3.

    Example:
        >>> W, rho = qsys_mg1_fb([0.3, 0.2], [1.0, 0.5], [1.0, 1.0])
    """
    # Convert to numpy arrays
    lambda_arr = np.asarray(lambda_vals, dtype=float).flatten()
    mu_arr = np.asarray(mu_vals, dtype=float).flatten()
    cs_arr = np.asarray(cs_vals, dtype=float).flatten()

    # Validate input lengths
    if not (len(lambda_arr) == len(mu_arr) == len(cs_arr)):
        raise ValueError("lambda, mu, and cs must have the same length")

    # Validate positive values
    if np.any(lambda_arr <= 0) or np.any(mu_arr <= 0) or np.any(cs_arr < 0):
        raise ValueError("lambda and mu must be positive, cs must be non-negative")

    K = len(lambda_arr)

    # Compute per-class utilizations
    rho_i = lambda_arr / mu_arr

    # Overall utilization
    rho_total = np.sum(rho_i)

    # Stability check
    if rho_total >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_total} >= 1")

    # Compute response times using FB formula
    W = np.zeros(K)

    for k in range(K):
        # Mean service time for this class (job size x)
        x = 1.0 / mu_arr[k]

        # For FB/LAS formula:
        # rho_x = lambda * integral_0^x F_bar(t)dt
        rho_x = 0.0
        for i in range(K):
            if cs_arr[i] == 1:  # Exponential case
                # integral_0^x exp(-mu_i*t)dt = (1 - exp(-mu_i*x)) / mu_i
                integral_Fbar = (1 - np.exp(-mu_arr[i] * x)) / mu_arr[i]
            else:
                # For non-exponential, approximate using mean and variance
                integral_Fbar = min(x, 1.0 / mu_arr[i])  # Bounded approximation
            rho_x += lambda_arr[i] * integral_Fbar

        # Numerator integral: lambda * integral_0^x t*F_bar(t)dt
        numerator = 0.0
        for i in range(K):
            mu_i = mu_arr[i]
            if cs_arr[i] == 1:  # Exponential case
                # integral_0^x t*exp(-mu_i*t)dt
                # = 1/mu_i^2 * (1 - exp(-mu_i*x)*(1 + mu_i*x))
                integral_tFbar = (1 - np.exp(-mu_i * x) * (1 + mu_i * x)) / (mu_i ** 2)
            else:
                # Approximation for non-exponential
                integral_tFbar = min(x ** 2 / 2, 1.0 / (mu_arr[i] ** 2))
            numerator += lambda_arr[i] * integral_tFbar

        # FB formula: E[T(x)] = numerator / (1-rho_x)^2 + x / (1-rho_x)
        if rho_x >= 1:
            W[k] = np.inf
        else:
            waiting_term = numerator / (1 - rho_x) ** 2
            service_term = x / (1 - rho_x)
            W[k] = waiting_term + service_term

    # Compute rhohat = Q/(1+Q) to match qsys convention
    Q = np.sum(lambda_arr * W)
    rho_out = Q / (1 + Q)

    return W, rho_out

def qsys_mg1_setf(
    lambda_vals: Union[np.ndarray, List[float]],
    mu_vals: Union[np.ndarray, List[float]],
    cs_vals: Union[np.ndarray, List[float]]
) -> Tuple[np.ndarray, float]:
    """
    Compute mean response time for M/G/1/SETF queue.

    Under SETF (Shortest Elapsed Time First), priority is based on a job's
    attained service (elapsed processing time), but unlike FB/LAS, jobs
    are not preempted. Once a job starts service, it runs to completion.

    SETF is the non-preemptive version of FB/LAS scheduling.

    For SETF, the mean response time follows a modified FB formula:
        E[T(x)]^SETF = E[T(x)]^FB + E[R] / (1 - rho_x)

    where E[R] is the mean residual service time.

    Args:
        lambda_vals: Array of arrival rates per class
        mu_vals: Array of service rates per class
        cs_vals: Array of coefficients of variation per class

    Returns:
        Tuple[np.ndarray, float]: (W, rho) where W is the vector of mean
            response times per class, and rho is the modified utilization.

    References:
        - M. Nuyens and A. Wierman, "The Foreground-Background queue: A survey",
          Performance Evaluation, 2008.
        - A. Wierman and M. Harchol-Balter, SIGMETRICS 2003.

    Example:
        >>> W, rho = qsys_mg1_setf([0.3, 0.2], [1.0, 0.5], [1.0, 1.0])
    """
    # Convert to numpy arrays
    lambda_arr = np.asarray(lambda_vals, dtype=float).flatten()
    mu_arr = np.asarray(mu_vals, dtype=float).flatten()
    cs_arr = np.asarray(cs_vals, dtype=float).flatten()

    # Validate input lengths
    if not (len(lambda_arr) == len(mu_arr) == len(cs_arr)):
        raise ValueError("lambda, mu, and cs must have the same length")

    # Validate positive values
    if np.any(lambda_arr <= 0) or np.any(mu_arr <= 0) or np.any(cs_arr < 0):
        raise ValueError("lambda and mu must be positive, cs must be non-negative")

    K = len(lambda_arr)

    # Compute per-class utilizations
    rho_i = lambda_arr / mu_arr

    # Overall utilization
    rho_total = np.sum(rho_i)

    # Stability check
    if rho_total >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_total} >= 1")

    # Compute mean residual service time for the mixture distribution
    # E[R] = sum_i (lambda_i / lambda_total) * E[S_i^2] / (2 * E[S_i])
    lambda_total = np.sum(lambda_arr)
    mean_residual = 0.0
    for i in range(K):
        p_i = lambda_arr[i] / lambda_total
        mean_s = 1.0 / mu_arr[i]
        mean_s2 = (1.0 + cs_arr[i] ** 2) / (mu_arr[i] ** 2)
        mean_residual += p_i * mean_s2 / (2.0 * mean_s)

    # Compute response times using SETF formula
    W = np.zeros(K)

    for k in range(K):
        # Mean service time for this class (job size x)
        x = 1.0 / mu_arr[k]

        # For SETF formula (similar to FB but with residual):
        # rho_x = lambda * integral_0^x F_bar(t)dt
        rho_x = 0.0
        for i in range(K):
            if cs_arr[i] == 1:  # Exponential case
                # integral_0^x exp(-mu_i*t)dt = (1 - exp(-mu_i*x)) / mu_i
                integral_Fbar = (1 - np.exp(-mu_arr[i] * x)) / mu_arr[i]
            else:
                # For non-exponential, approximate using mean and variance
                integral_Fbar = min(x, 1.0 / mu_arr[i])  # Bounded approximation
            rho_x += lambda_arr[i] * integral_Fbar

        # Numerator integral: lambda * integral_0^x t*F_bar(t)dt
        numerator = 0.0
        for i in range(K):
            mu_i = mu_arr[i]
            if cs_arr[i] == 1:  # Exponential case
                # integral_0^x t*exp(-mu_i*t)dt
                # = 1/mu_i^2 * (1 - exp(-mu_i*x)*(1 + mu_i*x))
                integral_tFbar = (1 - np.exp(-mu_i * x) * (1 + mu_i * x)) / (mu_i ** 2)
            else:
                # Approximation for non-exponential
                integral_tFbar = min(x ** 2 / 2, 1.0 / (mu_arr[i] ** 2))
            numerator += lambda_arr[i] * integral_tFbar

        # SETF formula: E[T(x)]^SETF = E[T(x)]^FB + E[R] / (1 - rho_x)
        # where E[T(x)]^FB = numerator / (1-rho_x)^2 + x / (1-rho_x)
        if rho_x >= 1:
            W[k] = np.inf
        else:
            # FB waiting term
            fb_waiting_term = numerator / (1 - rho_x) ** 2
            # FB service term (slowdown)
            fb_service_term = x / (1 - rho_x)
            # Non-preemptive penalty: residual service time adjusted
            np_penalty = mean_residual / (1 - rho_x)

            W[k] = fb_waiting_term + fb_service_term + np_penalty

    # Compute rhohat = Q/(1+Q) to match qsys convention
    Q = np.sum(lambda_arr * W)
    rho_out = Q / (1 + Q)

    return W, rho_out

def qsys_mg1_lrpt(
    lambda_vals: Union[np.ndarray, List[float]],
    mu_vals: Union[np.ndarray, List[float]],
    cs_vals: Union[np.ndarray, List[float]]
) -> Tuple[np.ndarray, float]:
    """
    Compute mean response time for M/G/1/LRPT queue.

    Under LRPT (Longest Remaining Processing Time), at any given point,
    the processor is shared evenly among all jobs with the longest
    remaining processing time. This is a remaining-size based policy
    that prioritizes large jobs.

    Classification (Wierman-Harchol-Balter 2003):
        LRPT is "Always Unfair" - for all finite job sizes y,
        E[S(y)]^LRPT > 1/(1-rho) under any service distribution.

    For LRPT, the expected slowdown for a job of size x is:
        E[S(x)]^LRPT = 1/(1-rho) + lambda*E[X^2] / (2*x*(1-rho)^2)

    Therefore the expected response time is:
        E[T(x)]^LRPT = x/(1-rho) + lambda*E[X^2] / (2*(1-rho)^2)

    Args:
        lambda_vals: Array of arrival rates per class
        mu_vals: Array of service rates per class
        cs_vals: Array of coefficients of variation per class

    Returns:
        Tuple[np.ndarray, float]: (W, rho) where W is the vector of mean
            response times per class, and rho is the modified utilization.

    References:
        A. Wierman and M. Harchol-Balter, SIGMETRICS 2003, Section 3.2.

    Example:
        >>> W, rho = qsys_mg1_lrpt([0.3, 0.2], [1.0, 0.5], [1.0, 1.0])
    """
    # Convert to numpy arrays
    lambda_arr = np.asarray(lambda_vals, dtype=float).flatten()
    mu_arr = np.asarray(mu_vals, dtype=float).flatten()
    cs_arr = np.asarray(cs_vals, dtype=float).flatten()

    # Validate input lengths
    if not (len(lambda_arr) == len(mu_arr) == len(cs_arr)):
        raise ValueError("lambda, mu, and cs must have the same length")

    # Validate positive values
    if np.any(lambda_arr <= 0) or np.any(mu_arr <= 0) or np.any(cs_arr < 0):
        raise ValueError("lambda and mu must be positive, cs must be non-negative")

    K = len(lambda_arr)

    # Compute per-class utilizations
    rho_i = lambda_arr / mu_arr

    # Overall utilization
    rho_total = np.sum(rho_i)

    # Stability check
    if rho_total >= 1:
        raise ValueError(f"System is unstable: utilization rho = {rho_total} >= 1")

    # Compute overall second moment of service time
    # E[X^2] = sum_i (lambda_i / lambda_total) * E[S_i^2]
    lambda_total = np.sum(lambda_arr)
    p = lambda_arr / lambda_total  # mixture probabilities

    E_X2 = 0.0
    for i in range(K):
        E_S2_i = (1 + cs_arr[i] ** 2) / (mu_arr[i] ** 2)
        E_X2 += p[i] * E_S2_i

    # The reference has TWO branches and only the exponential one was ported
    # here, so a class with cs != 1 was answered with the exponential formula:
    # on lambda=[0.3,0.2], mu=[1,0.5], cs=[1,sqrt(0.5)] that read
    # W = [13.33, 5.33] against MATLAB's [7.111, 3.333].
    if np.all(np.abs(cs_arr - 1.0) < 1e-6):
        # E[T(x)] = x/(1-rho) + lambda_total*E[X^2] / (2*(1-rho)^2)
        W = np.zeros(K)
        term2 = lambda_total * E_X2 / (2 * (1 - rho_total) ** 2)
        for k in range(K):
            W[k] = (1.0 / mu_arr[k]) / (1 - rho_total) + term2
    else:
        W = _lrpt_general(lambda_arr, mu_arr)

    # Compute rhohat = Q/(1+Q) to match qsys convention
    Q = np.sum(lambda_arr * W)
    rho_out = Q / (1 + Q)

    return W, rho_out


def _lrpt_general(lambda_arr: np.ndarray, mu_arr: np.ndarray) -> np.ndarray:
    """General (non-exponential) LRPT branch of MATLAB qsys_mg1_lrpt.

    The class-based preemptive priority surrogate, with the classes ordered by
    DECREASING mean service time, which is the order LRPT serves them in::

        W_q(k) = (sum_{j<=k} lambda_j/mu_j^2) / ((1 - rho_{<k})(1 - rho_{<=k}))
        W(k)   = W_q(k) + 1/mu_k

    The sort must be STABLE, as MATLAB's descending sort is, or two classes of
    equal mean size swap places and the cumulative sums change.
    """
    K = len(lambda_arr)
    mean_service = 1.0 / mu_arr
    order = sorted(range(K), key=lambda i: -mean_service[i])
    W = np.zeros(K)
    rho_prev = 0.0
    er_k = 0.0
    for r in order:
        rho_curr = rho_prev + lambda_arr[r] / mu_arr[r]
        er_k += lambda_arr[r] / (mu_arr[r] ** 2)
        W[r] = er_k / ((1 - rho_prev) * (1 - rho_curr)) + 1.0 / mu_arr[r]
        rho_prev = rho_curr
    return W
