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

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    def njit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator


def qsys_mg1_srpt(
    lambda_vals: Union[np.ndarray, List[float]],
    mu_vals: Union[np.ndarray, List[float]],
    cs_vals: Union[np.ndarray, List[float]]
) -> Tuple[np.ndarray, float]:
    """
    Compute mean response time for M/G/1/SRPT queue.

    Under SRPT (Shortest Remaining Processing Time), the job with the
    smallest remaining processing time receives priority. This is the
    optimal policy for minimizing mean response time.

    For exponential service (M/M/1/SRPT), SRPT reduces to preemptive priority
    queueing since the memoryless property makes remaining time distribution
    equal to the original distribution.

    Classification (Wierman-Harchol-Balter 2003):
        SRPT is "Sometimes Unfair" - it is unfair only for some loads
        or some service distributions.

    Args:
        lambda_vals: Array of arrival rates per class
        mu_vals: Array of service rates per class
        cs_vals: Array of coefficients of variation per class (cs=1 for exponential)

    Returns:
        Tuple[np.ndarray, float]: (W, rho) where W is the vector of mean
            response times per class, and rho is the modified utilization.

    References:
        - L. E. Schrage and L. W. Miller, "The queue M/G/1 with the shortest
          remaining processing time discipline", Operations Research, 1966.
        - A. Wierman and M. Harchol-Balter, SIGMETRICS 2003.

    Example:
        >>> W, rho = qsys_mg1_srpt([0.3, 0.2], [1.0, 0.5], [1.0, 1.0])
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

    # Sort classes by mean service time (ascending) for SRPT priority
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

    # Check if all cs values are approximately 1 (exponential case)
    is_exponential = np.all(np.abs(cs_sorted - 1) < 1e-10)

    if is_exponential or K == 1:
        # Exponential case or single class: use preemptive priority formula
        W_sorted = _qsys_mg1_srpt_exp(lambda_sorted, mu_sorted)
    else:
        # General case: use numerical integration for each class
        W_sorted = _qsys_mg1_srpt_general(lambda_sorted, mu_sorted, cs_sorted)

    # Restore original class ordering
    unsort_idx = np.argsort(sort_idx)
    W = W_sorted[unsort_idx]

    # Compute rhohat = Q/(1+Q) to match qsys convention
    Q = np.sum(lambda_arr * W)
    rho_out = Q / (1 + Q)

    return W, rho_out


def _qsys_mg1_srpt_exp(lambda_arr: np.ndarray, mu_arr: np.ndarray) -> np.ndarray:
    """SRPT for exponential service - uses preemptive priority formula."""
    K = len(lambda_arr)
    rho_i = lambda_arr / mu_arr

    # Compute mean residual service time
    # E[R] = sum_i lambda_i * E[S_i^2] / 2 = sum_i lambda_i / mu_i^2
    # For exponential: E[S^2] = 2/mu^2
    E_R = np.sum(lambda_arr / (mu_arr ** 2))

    # Compute response times using preemptive priority formula
    W = np.zeros(K)

    for k in range(K):
        # Cumulative utilization up to class k-1 (higher priority classes)
        rho_prev = np.sum(rho_i[:k]) if k > 0 else 0.0

        # Cumulative utilization up to class k
        rho_curr = np.sum(rho_i[:k+1])

        # Mean waiting time for class k (preemptive priority)
        W_q = E_R / ((1 - rho_prev) * (1 - rho_curr))

        # Response time = waiting time + service time
        W[k] = W_q + 1.0 / mu_arr[k]

    return W


def _qsys_mg1_srpt_general(
    lambda_arr: np.ndarray,
    mu_arr: np.ndarray,
    cs_arr: np.ndarray
) -> np.ndarray:
    """SRPT for general service distributions - uses numerical integration."""
    K = len(lambda_arr)
    W = np.zeros(K)

    for k in range(K):
        x = 1.0 / mu_arr[k]  # Mean service time for class k

        # Compute truncated load rho(x) - sum from smaller classes
        smaller_idx = (1.0 / mu_arr) <= x
        rho_x = np.sum(lambda_arr[smaller_idx] / mu_arr[smaller_idx])

        if rho_x >= 1:
            W[k] = np.inf
            continue

        # Numerator integral: lambda * integral_0^x t*F_bar(t)dt
        numerator = 0.0
        for i in range(K):
            if 1.0 / mu_arr[i] <= x:
                numerator += lambda_arr[i] / (mu_arr[i] ** 2)

        # Second integral: integral_0^x dt/(1-rho(t))
        n_steps = 100
        dt = x / n_steps
        second_term = 0.0
        for step in range(n_steps):
            t = (step + 0.5) * dt
            # Compute rho(t) at this point
            rho_t = 0.0
            for i in range(K):
                # For exponential: integral_0^t s*mu*exp(-mu*s)ds
                mu_i = mu_arr[i]
                rho_t += (lambda_arr[i] / mu_i) * (1 - np.exp(-mu_i * t) * (1 + mu_i * t))
            if rho_t < 1:
                second_term += dt / (1 - rho_t)

        # Schrage-Miller formula
        W[k] = numerator / (1 - rho_x) ** 2 + second_term

    return W


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

    # Compute response times using LRPT formula
    # E[T(x)] = x/(1-rho) + lambda_total*E[X^2] / (2*(1-rho)^2)
    W = np.zeros(K)

    for k in range(K):
        # Mean service time for this class
        x = 1.0 / mu_arr[k]

        # LRPT formula
        term1 = x / (1 - rho_total)
        term2 = lambda_total * E_X2 / (2 * (1 - rho_total) ** 2)
        W[k] = term1 + term2

    # Compute rhohat = Q/(1+Q) to match qsys convention
    Q = np.sum(lambda_arr * W)
    rho_out = Q / (1 + Q)

    return W, rho_out
