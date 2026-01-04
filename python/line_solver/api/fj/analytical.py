"""
Fork-Join analytical functions.

Native Python implementations of analytical methods for Fork-Join queueing
systems, including response time approximations, expected maximum calculations,
and order statistics.

References:
    A. Thomasian, "Analysis of Fork/Join and Related Queueing Systems",
    ACM Computing Surveys, Vol. 47, No. 2, Article 17, July 2014.
"""

import numpy as np
from typing import Tuple, Optional, Callable, Dict, Union
from dataclasses import dataclass
from scipy.special import comb, factorial
from scipy.integrate import quad
from scipy.optimize import brentq


# =============================================================================
# Basic Functions
# =============================================================================

def fj_harmonic(K: int) -> float:
    """
    Compute Harmonic sum H_K = sum(1/k) for k=1 to K.

    The K-th Harmonic number is fundamental in Fork-Join analysis.
    For large K, H_K ~ ln(K) + gamma, where gamma ~ 0.57721 is the
    Euler-Mascheroni constant.

    Args:
        K: Number of parallel servers (positive integer)

    Returns:
        Harmonic sum H_K = 1 + 1/2 + 1/3 + ... + 1/K

    References:
        Thomasian, Table I and throughout.
    """
    if K < 1:
        raise ValueError(f'K must be a positive integer. Got K={K}.')
    return np.sum(1.0 / np.arange(1, K + 1))


def fj_synch_delay(lambda_val: float, mu: float) -> float:
    """
    Synchronization delay S_2(rho) for two-way Fork-Join.

    Computes the mean synchronization delay - the average time tasks spend
    waiting at synchronization queues for their siblings to complete.

    S_2(rho) = (1/2) * (1 - rho/4) * R(rho)

    The F/J response time decomposes as:
        R_2^{F/J}(rho) = R(rho) + S_2(rho)

    Args:
        lambda_val: Arrival rate
        mu: Service rate (mu > lambda for stability)

    Returns:
        Mean synchronization delay S_2(rho)

    References:
        Thomasian, Eq. (8) on page 17:11.
    """
    rho = lambda_val / mu
    if rho >= 1:
        raise ValueError(f'System is unstable: rho = {rho:.4f} >= 1.')

    # M/M/1 mean response time
    R_rho = 1.0 / (mu - lambda_val)

    # Synchronization delay (Eq. 8)
    S = 0.5 * (1 - rho / 4) * R_rho
    return S


# =============================================================================
# Response Time Approximations
# =============================================================================

def fj_respt_2way(lambda_val: float, mu: float) -> float:
    """
    Exact two-way Fork-Join response time R_2^{F/J}(rho).

    Computes the exact mean response time for a 2-way (K=2) Fork-Join
    queueing system with Poisson arrivals and exponential service times.

    R_2^{F/J}(rho) = (H_2 - rho/8) * R(rho) = (12 - rho)/8 * R(rho)

    where H_2 = 1.5 and R(rho) = (mu - lambda)^{-1}.

    Args:
        lambda_val: Arrival rate
        mu: Service rate (mu > lambda for stability)

    Returns:
        Exact 2-way F/J response time

    References:
        Thomasian, Eq. (6) on page 17:10.
        Flatto and Hahn, SIAM J. Appl. Math., 1984.
    """
    rho = lambda_val / mu
    if rho >= 1:
        raise ValueError(f'System is unstable: rho = {rho:.4f} >= 1.')

    R_rho = 1.0 / (mu - lambda_val)
    R = (12 - rho) / 8 * R_rho
    return R


def fj_respt_nt(K: int, lambda_val: float, mu: float) -> float:
    """
    Nelson-Tantawi approximation for K-way F/J response time.

    Valid for 2 <= K <= 32.

    R_K^{F/J}(rho) ~ [H_K/H_2 + (1 - H_K/H_2) * 4*rho/11] * (1.5 - rho/8) * R(rho)

    Args:
        K: Number of parallel servers (2 <= K <= 32)
        lambda_val: Arrival rate
        mu: Service rate (mu > lambda for stability)

    Returns:
        Approximate K-way F/J response time

    References:
        Thomasian, Eq. (3) on page 17:10.
        Nelson and Tantawi, IEEE Trans. Computers, 1988.
    """
    if K < 2:
        raise ValueError(f'Nelson-Tantawi requires K >= 2. Got K={K}.')

    rho = lambda_val / mu
    if rho >= 1:
        raise ValueError(f'System is unstable: rho = {rho:.4f} >= 1.')

    H_K = fj_harmonic(K)
    H_2 = fj_harmonic(2)  # H_2 = 1.5

    S_K = H_K / H_2 + (1 - H_K / H_2) * (4 * rho / 11)
    R_2_factor = 1.5 - rho / 8
    R_rho = 1.0 / (mu - lambda_val)

    return S_K * R_2_factor * R_rho


def fj_respt_vm(K: int, lambda_val: float, mu: float) -> float:
    """
    Varma-Makowski approximation for K-way F/J response time.

    R_K^{F/J}(rho) ~ [H_K + (A_K - H_K) * rho] * R(rho)

    where A_K = sum_{i=1}^{K} C(K,i) * (-1)^{i-1} * sum_{m=1}^{i} C(i,m) * (m-1)! / i^{m+1}

    Args:
        K: Number of parallel servers (positive integer)
        lambda_val: Arrival rate
        mu: Service rate (mu > lambda for stability)

    Returns:
        Approximate K-way F/J response time

    References:
        Thomasian, Eq. (4) on page 17:10.
        Varma and Makowski, Performance Evaluation, 1994.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')

    rho = lambda_val / mu
    if rho >= 1:
        raise ValueError(f'System is unstable: rho = {rho:.4f} >= 1.')

    H_K = fj_harmonic(K)

    # Compute A_K
    A_K = 0.0
    for i in range(1, K + 1):
        inner_sum = 0.0
        for m in range(1, i + 1):
            inner_sum += comb(i, m, exact=True) * factorial(m - 1) / (i ** (m + 1))
        A_K += comb(K, i, exact=True) * ((-1) ** (i - 1)) * inner_sum

    R_rho = 1.0 / (mu - lambda_val)
    return (H_K + (A_K - H_K) * rho) * R_rho


def fj_respt_varki(K: int, lambda_val: float, mu: float) -> float:
    """
    Varki et al. approximation for K-way F/J response time.

    Mean of pessimistic (upper) and optimistic (lower) bounds:

    R_K^{F/J}(rho) ~ (1/mu) * [H_K + (rho/(2*(1-rho))) * (S1 + (1-2*rho)*S2)]

    where:
        S1 = sum_{i=1}^{K} 1/(i - rho)
        S2 = sum_{i=1}^{K} 1/(i*(i - rho))

    Args:
        K: Number of parallel servers (positive integer)
        lambda_val: Arrival rate
        mu: Service rate (mu > lambda for stability)

    Returns:
        Approximate K-way F/J response time

    References:
        Thomasian, Eq. (5) on page 17:10.
        Varki et al., IEEE TPDS, 2014.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')

    rho = lambda_val / mu
    if rho >= 1:
        raise ValueError(f'System is unstable: rho = {rho:.4f} >= 1.')

    H_K = fj_harmonic(K)
    indices = np.arange(1, K + 1)

    S1 = np.sum(1.0 / (indices - rho))
    S2 = np.sum(1.0 / (indices * (indices - rho)))

    R = (1.0 / mu) * (H_K + (rho / (2 * (1 - rho))) * (S1 + (1 - 2 * rho) * S2))
    return R


def fj_rmax(K: int, lambda_val: float, mu: float) -> float:
    """
    Maximum response time R_K^max(rho) for K M/M/1 queues.

    Upper bound for K-way Fork-Join response time:

    R_K^max(rho) = H_K * R(rho) = H_K / (mu - lambda)

    Args:
        K: Number of parallel servers (positive integer)
        lambda_val: Arrival rate
        mu: Service rate (mu > lambda for stability)

    Returns:
        Expected maximum response time R_K^max(rho)

    References:
        Thomasian, Eq. (1) on page 17:8.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')

    rho = lambda_val / mu
    if rho >= 1:
        raise ValueError(f'System is unstable: rho = {rho:.4f} >= 1.')

    H_K = fj_harmonic(K)
    R_rho = 1.0 / (mu - lambda_val)
    return H_K * R_rho


def fj_bounds(K: int, lambda_val: float, mu: float) -> Tuple[float, float]:
    """
    Upper and lower bounds for K-way F/J response time.

    Upper bound (pessimistic):
        R_K^max(rho) = H_K / (mu * (1 - rho))

    Lower bound (optimistic):
        R_K^{F/J(opt)}(rho) = (1/mu) * [H_K + S_{K(K-rho)}]
        where S_{K(K-rho)} = sum_{j=1}^{K} (1/j) * (rho/(j - rho))

    Args:
        K: Number of parallel servers (positive integer)
        lambda_val: Arrival rate
        mu: Service rate (mu > lambda for stability)

    Returns:
        Tuple of (Rmax, Rmin)

    References:
        Thomasian, Eq. (1) and (2) on page 17:9.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')

    rho = lambda_val / mu
    if rho >= 1:
        raise ValueError(f'System is unstable: rho = {rho:.4f} >= 1.')

    H_K = fj_harmonic(K)

    # Upper bound
    Rmax = H_K / (mu * (1 - rho))

    # Lower bound
    indices = np.arange(1, K + 1)
    S_K = np.sum((1.0 / indices) * (rho / (indices - rho)))
    Rmin = (1.0 / mu) * (H_K + S_K)

    return Rmax, Rmin


def fj_rmax_evd(K: int, R: float, sigma_R: float, calibrated: bool = False) -> float:
    """
    Maximum response time using Extreme Value Distribution approximation.

    R_K^max(rho) = R(rho) + (sqrt(6)*ln(K)/pi) * sigma_R(rho)

    Args:
        K: Number of parallel servers (positive integer)
        R: Mean response time
        sigma_R: Standard deviation of response time
        calibrated: Use calibrated formula with 1.27 divisor (default: False)

    Returns:
        Approximate maximum response time

    References:
        Thomasian, Eq. (43) on page 17:21.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')
    if R <= 0:
        raise ValueError('Mean response time R must be positive.')
    if sigma_R < 0:
        raise ValueError('Standard deviation sigma_R must be non-negative.')

    correction_factor = np.sqrt(6) * np.log(K) / np.pi
    if calibrated:
        correction_factor /= 1.27

    return R + correction_factor * sigma_R


def fj_rmax_erlang(K: int, k: int, lambda_val: float, mu: float) -> float:
    """
    Maximum response time R_K^max for Erlang service times.

    For K M/E_k/1 queues (Poisson arrivals, Erlang-k service times).

    Args:
        K: Number of parallel servers (positive integer)
        k: Number of Erlang stages (positive integer)
        lambda_val: Arrival rate
        mu: Service rate per Erlang stage (mean service = k/mu)

    Returns:
        Expected maximum response time

    References:
        Thomasian, Eq. (33) and (34) on page 17:18.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')
    if k < 1:
        raise ValueError(f'k (stages) must be positive. Got k={k}.')

    mean_service = k / mu
    rho = lambda_val * mean_service
    if rho >= 1:
        raise ValueError(f'System is unstable: rho = {rho:.4f} >= 1.')

    # M/E_k/1 response time using Pollaczek-Khinchin
    cv2 = 1.0 / k
    R_single = mean_service * (1 + rho * (1 + cv2) / (2 * (1 - rho)))

    if K == 2 and k <= 10:
        # Closed-form for K=2
        mu_resp = k / R_single
        correction = 0.0
        for m in range(k):
            for n in range(k):
                correction += comb(m + n, m, exact=True) * (mu_resp ** m) * (mu_resp ** n) / ((2 * mu_resp) ** (m + n + 1))
        return 2 * R_single - correction
    else:
        # Numerical integration for general K
        def erlang_cdf(t):
            if t <= 0:
                return 0.0
            k_resp = max(1, int(np.ceil(k / (1 + rho))))
            mu_resp = k_resp / R_single
            S = sum((mu_resp * t) ** j / factorial(j) for j in range(k_resp))
            return 1 - np.exp(-mu_resp * t) * S

        def integrand(t):
            return 1 - erlang_cdf(t) ** K

        result, _ = quad(integrand, 0, R_single * 20, limit=100)
        return result


# =============================================================================
# Expected Maximum Functions (X_K^max)
# =============================================================================

def fj_xmax_exp(K: int, mu: float) -> float:
    """
    Expected maximum of K i.i.d. exponential random variables.

    X_K^max = H_K / mu = H_K * x

    where H_K is the K-th Harmonic number.

    Args:
        K: Number of parallel servers (positive integer)
        mu: Service rate (mean service time is 1/mu)

    Returns:
        Expected maximum service time

    References:
        Thomasian, Eq. (10) on page 17:11.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')
    if mu <= 0:
        raise ValueError(f'Rate mu must be positive. Got mu={mu}.')

    return fj_harmonic(K) / mu


def fj_xmax_2(lambda1: float, lambda2: Optional[float] = None) -> float:
    """
    Expected maximum of 2 exponential random variables.

    Y_2^max = 1/lambda1 + 1/lambda2 - 1/(lambda1 + lambda2)

    For identical rates: Y_2^max = 1.5/lambda = H_2/lambda

    Args:
        lambda1: Rate of first exponential
        lambda2: Rate of second exponential (default: same as lambda1)

    Returns:
        Expected maximum Y_2^max

    References:
        Thomasian, Eq. (27) on page 17:17.
    """
    if lambda2 is None:
        lambda2 = lambda1

    if lambda1 <= 0 or lambda2 <= 0:
        raise ValueError(f'Rates must be positive. Got lambda1={lambda1}, lambda2={lambda2}.')

    return 1 / lambda1 + 1 / lambda2 - 1 / (lambda1 + lambda2)


def fj_xmax_erlang(K: int, k: int, mu: float) -> float:
    """
    Expected maximum of K i.i.d. Erlang-k random variables.

    Args:
        K: Number of parallel servers (positive integer)
        k: Number of Erlang stages (positive integer, default=2)
        mu: Rate parameter per stage (mean service time = k/mu)

    Returns:
        Expected maximum service time

    References:
        Thomasian, page 17:13.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')
    if k < 1:
        raise ValueError(f'k (stages) must be positive. Got k={k}.')
    if mu <= 0:
        raise ValueError(f'Rate mu must be positive.')

    if k == 2:
        # Closed-form for Erlang-2
        outer_sum = 0.0
        for n in range(1, K + 1):
            inner_sum = 0.0
            for m in range(1, n + 1):
                inner_sum += comb(n, m, exact=True) * factorial(m) / (2 * n ** (m + 1))
            outer_sum += comb(K, n, exact=True) * ((-1) ** (n - 1)) * inner_sum
        return outer_sum / mu
    else:
        # Numerical integration
        def erlang_cdf(t):
            if t <= 0:
                return 0.0
            S = sum((mu * t) ** j / factorial(j) for j in range(k))
            return 1 - np.exp(-mu * t) * S

        def integrand(t):
            return 1 - erlang_cdf(t) ** K

        mean_erlang = k / mu
        upper_limit = mean_erlang * 10 + 10 * np.sqrt(k) / mu
        result, _ = quad(integrand, 0, upper_limit, limit=100)
        return result


def fj_xmax_hyperexp(K: int, p1: float, mu1: float, mu2: float) -> float:
    """
    Expected maximum of K i.i.d. Hyperexponential-2 random variables.

    X_K^max = sum_{n=1}^{K} (-1)^{n+1} * sum_{m=0}^{n} C(n,m) * p1^m * p2^{n-m} / (m*mu1 + (n-m)*mu2)

    Args:
        K: Number of parallel servers (positive integer)
        p1: Probability of branch 1 (0 < p1 < 1)
        mu1: Rate of branch 1
        mu2: Rate of branch 2

    Returns:
        Expected maximum service time

    References:
        Thomasian, Eq. (15) on page 17:13.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')
    if not 0 < p1 < 1:
        raise ValueError(f'p1 must be in (0,1). Got p1={p1}.')
    if mu1 <= 0 or mu2 <= 0:
        raise ValueError('Rates mu1 and mu2 must be positive.')

    p2 = 1 - p1
    Xmax = 0.0
    for n in range(1, K + 1):
        inner_sum = 0.0
        for m in range(n + 1):
            denom = m * mu1 + (n - m) * mu2
            if denom > 0:
                inner_sum += comb(n, m, exact=True) * (p1 ** m) * (p2 ** (n - m)) / denom
        Xmax += ((-1) ** (n + 1)) * inner_sum

    return Xmax


def fj_xmax_approx(
    K: int,
    mu_X: float,
    sigma_X: float,
    dist_type: str = 'exp'
) -> Tuple[float, float]:
    """
    General approximation for expected maximum of K random variables.

    X_K^max ~ mu_X + sigma_X * G(K)

    where G(K) depends on distribution type:
        - Exponential: G(K) = H_K - 1
        - Uniform:     G(K) = sqrt(3) * (K-1) / (K+1)
        - EVD:         G(K) = sqrt(6) * ln(K) / pi
        - General:     G(K) <= (K-1) / sqrt(2K-1)  (upper bound)

    Args:
        K: Number of random variables (positive integer)
        mu_X: Mean of the distribution
        sigma_X: Standard deviation of the distribution
        dist_type: Distribution type ('exp', 'uniform', 'evd', 'bound')

    Returns:
        Tuple of (Xmax, GK)

    References:
        Thomasian, Eq. (22) and (23) on page 17:16.
        David, "Order Statistics", Wiley, 1970.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')
    if sigma_X < 0:
        raise ValueError('Standard deviation sigma_X must be non-negative.')

    dist_type = dist_type.lower()
    if dist_type == 'exp':
        GK = fj_harmonic(K) - 1
    elif dist_type == 'uniform':
        GK = np.sqrt(3) * (K - 1) / (K + 1)
    elif dist_type == 'evd':
        GK = np.sqrt(6) * np.log(K) / np.pi
    elif dist_type == 'bound':
        GK = (K - 1) / np.sqrt(2 * K - 1)
    else:
        raise ValueError(f'Unknown distribution type: {dist_type}')

    Xmax = mu_X + sigma_X * GK
    return Xmax, GK


def fj_xmax_emma(K: int, param: Union[float, Callable], dist_type: str = 'auto') -> float:
    """
    Expected maximum using EMMA (Expected Maximum from Marginal Approximation).

    Based on: [F_X(E[Y_K])]^K ~ phi = 0.570376

    Therefore: E[Y_K] ~ F^{-1}(phi^{1/K})

    For exponential with rate mu: E[Y_K] = -(1/mu) * ln(1 - phi^{1/K})

    Args:
        K: Number of random variables (positive integer)
        param: Rate mu for exponential, or inverse CDF function handle
        dist_type: 'exp', 'general', or 'auto' (detect from param type)

    Returns:
        Approximate expected maximum

    References:
        Thomasian, Eq. (44) on page 17:22.
    """
    phi = 0.570376  # exp(-exp(-gamma)) where gamma is Euler-Mascheroni

    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')

    if dist_type == 'auto':
        dist_type = 'general' if callable(param) else 'exp'

    if dist_type == 'exp':
        mu = param
        if mu <= 0:
            raise ValueError('Rate mu must be positive.')
        return -(1 / mu) * np.log(1 - phi ** (1 / K))
    else:
        F_inv = param
        return F_inv(phi ** (1 / K))


def fj_xmax_normal(
    K: int,
    mu: float,
    sigma: float,
    method: str = 'johnson'
) -> Tuple[float, float]:
    """
    Expected maximum for normal distribution.

    Args:
        K: Number of random variables (K >= 2)
        mu: Mean of the normal distribution
        sigma: Standard deviation
        method: 'johnson' (default), 'arnold', or 'corrected'

    Returns:
        Tuple of (Xmax, Vmax) - expected value and variance of maximum

    References:
        Thomasian, Eq. (45) on page 17:23.
    """
    if K < 2:
        raise ValueError(f'K must be >= 2 for normal approximation. Got K={K}.')
    if sigma < 0:
        raise ValueError('Standard deviation must be non-negative.')

    gamma_em = 0.5772156649015329  # Euler-Mascheroni
    sqrt_2lnK = np.sqrt(2 * np.log(K))

    method = method.lower()
    if method == 'arnold':
        GK = sqrt_2lnK
    elif method == 'johnson':
        correction = (np.log(np.log(K)) - np.log(4 * np.pi) + 2 * gamma_em) / (2 * sqrt_2lnK)
        GK = sqrt_2lnK - correction
    elif method == 'corrected':
        correction = (np.log(np.log(K)) - np.log(4 * np.pi) + 2 * gamma_em) / (2 * sqrt_2lnK)
        delta_K = 0.1727 * K ** (-0.2750)
        GK = sqrt_2lnK - correction - delta_K
    else:
        raise ValueError(f'Unknown method: {method}')

    Xmax = mu + sigma * GK
    Vmax = 1.64492 * sigma ** 2 / (2 * np.log(K))
    return Xmax, Vmax


def fj_xmax_pareto(K: int, beta: float, k: Optional[float] = None) -> Tuple[float, float]:
    """
    Expected maximum for Pareto distribution.

    Args:
        K: Number of random variables (positive integer)
        beta: Shape parameter (beta > 2 for finite moments)
        k: Scale parameter (default: beta - 1 for mean = 1)

    Returns:
        Tuple of (Xmax, MK) - expected maximum and characteristic maximum

    References:
        Thomasian, page 17:24.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')
    if beta <= 2:
        raise ValueError(f'beta must be > 2 for finite moments. Got beta={beta}.')

    if k is None:
        k = beta - 1

    if k <= 0:
        raise ValueError('Scale parameter k must be positive.')

    # m_K = k * (K^(1/beta) - 1)
    mK = k * (K ** (1 / beta) - 1)

    # Numerical integration for E[max]
    def integrand(x):
        F = 1 - (k / (k + max(x, 0))) ** beta
        return 1 - F ** K

    upper_limit = k * K ** (2 / beta) * 10
    Xmax, _ = quad(integrand, 0, upper_limit, limit=100)

    # Characteristic maximum
    if beta > 1:
        tail_integral = k ** beta * (k + mK) ** (1 - beta) / (beta - 1)
        MK = mK + K * tail_integral
    else:
        MK = np.inf

    return Xmax, MK


# =============================================================================
# Order Statistics and Bounds
# =============================================================================

def fj_order_stat(
    y: Union[float, np.ndarray],
    k: int,
    K: int,
    F_X: Callable[[float], float]
) -> Tuple[np.ndarray, Optional[float]]:
    """
    CDF and expected value of k-th order statistic.

    For maximum (k=K): F_{Y_K}(y) = [F_X(y)]^K

    For k-th order statistic:
        F_{Y_k}(y) = sum_{j=k}^{K} C(K,j) * F_X(y)^j * (1-F_X(y))^{K-j}

    Args:
        y: Value(s) at which to evaluate CDF
        k: Order of the statistic (1 = minimum, K = maximum)
        K: Total number of random variables
        F_X: CDF function handle

    Returns:
        Tuple of (F_Yk, E_Yk) - CDF values and expected value

    References:
        Thomasian, Eq. (18) and (19) on page 17:15.
    """
    if k < 1 or k > K:
        raise ValueError(f'k must satisfy 1 <= k <= K. Got k={k}, K={K}.')

    y = np.atleast_1d(y)
    F_y = np.array([F_X(yi) for yi in y])

    if k == K:
        F_Yk = F_y ** K
    else:
        F_Yk = np.zeros_like(F_y)
        for j in range(k, K + 1):
            F_Yk += comb(K, j, exact=True) * (F_y ** j) * ((1 - F_y) ** (K - j))

    # Compute expected value
    def find_upper_limit(F_X, target=0.999999):
        u = 1.0
        for _ in range(100):
            if F_X(u) >= target:
                break
            u *= 2
        return u

    if k == K:
        integrand = lambda t: 1 - F_X(t) ** K
    elif k == 1:
        integrand = lambda t: (1 - F_X(t)) ** K
    else:
        eps_val = 1e-8

        def f_X(t):
            return (F_X(t + eps_val) - F_X(t - eps_val)) / (2 * eps_val)

        coeff = K * comb(K - 1, k - 1, exact=True)
        integrand = lambda t: t * coeff * f_X(t) * (F_X(t) ** (k - 1)) * ((1 - F_X(t)) ** (K - k))

    upper_limit = find_upper_limit(F_X)
    E_Yk, _ = quad(integrand, 0, upper_limit, limit=100)

    return F_Yk, E_Yk


def fj_gk_bound(K: int, gk_type: str = 'all') -> Union[float, Dict[str, float]]:
    """
    Compute G(K) factors for expected maximum approximation.

    G(K) factors for X_K^max ~ mu_X + sigma_X * G(K):
        - Exponential: G(K) = H_K - 1
        - Uniform:     G(K) = sqrt(3) * (K-1) / (K+1)
        - EVD:         G(K) = sqrt(6) * ln(K) / pi
        - Upper bound: G(K) <= (K-1) / sqrt(2K-1)

    Args:
        K: Number of random variables (positive integer)
        gk_type: 'exp', 'uniform', 'evd', 'bound', or 'all' (default)

    Returns:
        G(K) value or dict with all values

    References:
        Thomasian, Eq. (22) and (23) on page 17:16.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')

    gk_type = gk_type.lower()

    if gk_type == 'exp':
        return fj_harmonic(K) - 1
    elif gk_type == 'uniform':
        return np.sqrt(3) * (K - 1) / (K + 1)
    elif gk_type == 'evd':
        return np.sqrt(6) * np.log(K) / np.pi
    elif gk_type == 'bound':
        return (K - 1) / np.sqrt(2 * K - 1)
    elif gk_type == 'all':
        return {
            'K': K,
            'exponential': fj_harmonic(K) - 1,
            'uniform': np.sqrt(3) * (K - 1) / (K + 1),
            'evd': np.sqrt(6) * np.log(K) / np.pi,
            'upper_bound': (K - 1) / np.sqrt(2 * K - 1),
        }
    else:
        raise ValueError(f'Unknown type: {gk_type}')


def fj_char_max(
    K: int,
    param: Union[float, Tuple[int, float], Callable],
    dist_type: str = 'exp'
) -> Tuple[float, float]:
    """
    Characteristic maximum M_K for order statistics.

    M_K provides a bound for the expected maximum of K i.i.d. random variables.

    For exponential with rate mu:
        m_K = ln(K) / mu
        M_K = H_K / mu

    Args:
        K: Number of random variables (positive integer)
        param: Rate mu for exponential, [k, mu] for Erlang, or survival function
        dist_type: 'exp' (default), 'erlang', or 'general'

    Returns:
        Tuple of (MK, mK) - characteristic maximum and threshold

    References:
        Thomasian, Eq. (47) on page 17:24.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')

    dist_type = dist_type.lower()

    if dist_type == 'exp':
        mu = param
        if mu <= 0:
            raise ValueError('Rate mu must be positive.')
        mK = np.log(K) / mu
        MK = fj_harmonic(K) / mu

    elif dist_type == 'erlang':
        k_stages, mu = param
        if k_stages < 1:
            raise ValueError('Erlang stages must be positive.')
        if mu <= 0:
            raise ValueError('Rate mu must be positive.')

        # Erlang survival function
        def S_erlang(x):
            if x <= 0:
                return 1.0
            S = sum((mu * x) ** i / factorial(i) for i in range(k_stages))
            return np.exp(-mu * x) * S

        # Find m_K where S(m_K) = 1/K
        initial_guess = k_stages / mu + np.log(K) / mu
        mK = brentq(lambda x: S_erlang(x) - 1 / K, 0, initial_guess * 10)
        MK = (k_stages / mu) * (1 + K * np.exp(-mu * mK) * (mu * mK) ** k_stages / factorial(k_stages))

    elif dist_type == 'general':
        S_X = param  # Survival function
        # Find m_K where S_X(m_K) = 1/K
        upper = 1.0
        while S_X(upper) > 1 / K:
            upper *= 2
        mK = brentq(lambda x: S_X(x) - 1 / K, 0, upper)

        # Find upper limit for integration
        u = mK
        while S_X(u) > 1e-10:
            u *= 2

        tail_integral, _ = quad(S_X, mK, u)
        MK = mK + K * tail_integral

    else:
        raise ValueError(f'Unknown distribution type: {dist_type}')

    return MK, mK


def fj_quantile(
    K: int,
    q: Union[float, np.ndarray],
    F_inv: Optional[Callable[[float], float]] = None
) -> Union[float, np.ndarray]:
    """
    Quantile approximation for maximum of K random variables.

    Standard Gumbel approximation: x(K,q) ~ ln(K) - ln(ln(1/q))

    For general distributions: x(K,q) = F^{-1}(q^{1/K})

    Args:
        K: Number of random variables (positive integer)
        q: Quantile probability (0 < q < 1)
        F_inv: Optional inverse CDF for general distributions

    Returns:
        q-th quantile of the maximum

    References:
        Thomasian, page 17:22.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')

    q = np.atleast_1d(q)
    if np.any(q <= 0) or np.any(q >= 1):
        raise ValueError('Quantile q must satisfy 0 < q < 1.')

    if F_inv is None:
        # Standard Gumbel approximation
        result = np.log(K) - np.log(np.log(1.0 / q))
    else:
        # General distribution
        result = np.array([F_inv(qi ** (1 / K)) for qi in q])

    return result if len(result) > 1 else result[0]


def fj_sm_tput(K: int, mu: float) -> float:
    """
    Maximum throughput for Split-Merge queueing system.

    In an SM system, all tasks must complete before the next request.

    lambda_K^SM = 1 / X_K^max = mu / H_K

    For comparison, the F/J maximum throughput is mu (factor H_K higher).

    Args:
        K: Number of parallel servers (positive integer)
        mu: Service rate at each server

    Returns:
        Maximum throughput lambda_K^SM

    References:
        Thomasian, page 17:2.
    """
    if K < 1:
        raise ValueError(f'K must be positive. Got K={K}.')
    if mu <= 0:
        raise ValueError(f'Service rate mu must be positive. Got mu={mu}.')

    return mu / fj_harmonic(K)
