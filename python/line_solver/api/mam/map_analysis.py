"""
Markovian Arrival Process (MAP) analysis algorithms.

Native Python implementations for analyzing MAPs, including:
- Steady-state distributions (map_piq, map_pie)
- Arrival rate and moments (map_lambda, map_mean, map_var, map_scv)
- MAP transformations and utilities

References:
    Neuts, M.F. "Matrix-Geometric Solutions in Stochastic Models:
    An Algorithmic Approach." Dover Publications, 1994.
"""

import numpy as np
from numpy.linalg import LinAlgError
from scipy import linalg
from typing import Tuple, Optional, Union

# Import CTMC solver from mc module
from ..mc import ctmc_solve


def map_infgen(D0: np.ndarray, D1: np.ndarray) -> np.ndarray:
    """
    Compute the infinitesimal generator of a MAP.

    The generator Q = D0 + D1 represents the underlying CTMC.

    Args:
        D0: Hidden transition matrix (non-arrival transitions)
        D1: Visible transition matrix (arrival transitions)

    Returns:
        Infinitesimal generator matrix Q = D0 + D1
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    return D0 + D1


def map_piq(D0: np.ndarray, D1: np.ndarray = None) -> np.ndarray:
    """
    Compute steady-state distribution of the underlying CTMC of a MAP.

    Solves πQ = 0 where Q = D0 + D1 is the generator.

    Args:
        D0: Hidden transition matrix, or stacked [D0, D1] if D1 is None
        D1: Visible transition matrix (optional)

    Returns:
        Steady-state probability vector π
    """
    if D1 is None:
        # D0 is actually [D0, D1] stacked
        D0_arr = np.asarray(D0)
        if D0_arr.ndim == 3 and D0_arr.shape[0] == 2:
            D0, D1 = D0_arr[0], D0_arr[1]
        else:
            raise ValueError("D0 must be a (2, n, n) array when D1 is not provided")

    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    Q = map_infgen(D0, D1)
    pi = ctmc_solve(Q)

    return pi


def map_lambda(D0: np.ndarray, D1: np.ndarray = None) -> float:
    """
    Compute the arrival rate (λ) of a MAP.

    The arrival rate is λ = π * D1 * e where π is the steady-state
    and e is the column vector of ones.

    Args:
        D0: Hidden transition matrix, or stacked [D0, D1] if D1 is None
        D1: Visible transition matrix (optional)

    Returns:
        Arrival rate λ
    """
    if D1 is None:
        D0_arr = np.asarray(D0)
        if D0_arr.ndim == 3 and D0_arr.shape[0] == 2:
            D0, D1 = D0_arr[0], D0_arr[1]
        else:
            raise ValueError("D0 must be a (2, n, n) array when D1 is not provided")

    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    n = D0.shape[0]
    e = np.ones(n)

    pi = map_piq(D0, D1)
    lambda_rate = pi @ D1 @ e

    return float(lambda_rate)


def map_pie(D0: np.ndarray, D1: np.ndarray = None) -> np.ndarray:
    """
    Compute equilibrium distribution of embedded DTMC.

    The embedded DTMC has transition matrix P = (-D0)^{-1} * D1.
    Its steady-state is π_e = π * D1 / (π * D1 * e).

    Args:
        D0: Hidden transition matrix, or stacked [D0, D1] if D1 is None
        D1: Visible transition matrix (optional)

    Returns:
        Equilibrium distribution of embedded DTMC
    """
    if D1 is None:
        D0_arr = np.asarray(D0)
        if D0_arr.ndim == 3 and D0_arr.shape[0] == 2:
            D0, D1 = D0_arr[0], D0_arr[1]
        else:
            raise ValueError("D0 must be a (2, n, n) array when D1 is not provided")

    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    n = D0.shape[0]
    e = np.ones(n)

    pi = map_piq(D0, D1)
    A = pi @ D1  # Row vector

    # Normalize
    normalizer = A @ e
    if normalizer > 0:
        pie = A / normalizer
    else:
        pie = np.ones(n) / n

    return pie


def map_mean(D0: np.ndarray, D1: np.ndarray = None) -> float:
    """
    Compute mean inter-arrival time of a MAP.

    The mean is 1/λ where λ is the arrival rate.

    Args:
        D0: Hidden transition matrix, or stacked [D0, D1] if D1 is None
        D1: Visible transition matrix (optional)

    Returns:
        Mean inter-arrival time
    """
    lambda_rate = map_lambda(D0, D1)
    if lambda_rate > 0:
        return 1.0 / lambda_rate
    else:
        return float('inf')


def map_var(D0: np.ndarray, D1: np.ndarray = None) -> float:
    """
    Compute variance of inter-arrival times of a MAP.

    Var[X] = E[X²] - E[X]²

    Uses map_moment for consistency with Kotlin implementation.

    Args:
        D0: Hidden transition matrix, or stacked [D0, D1] if D1 is None
        D1: Visible transition matrix (optional)

    Returns:
        Variance of inter-arrival times
    """
    if D1 is None:
        D0_arr = np.asarray(D0)
        if D0_arr.ndim == 3 and D0_arr.shape[0] == 2:
            D0, D1 = D0_arr[0], D0_arr[1]
        else:
            raise ValueError("D0 must be a (2, n, n) array when D1 is not provided")

    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    # Variance = E[X²] - E[X]² using map_moment and map_mean
    mean = map_mean(D0, D1)
    moment2 = map_moment(D0, D1, 2)
    variance = moment2 - mean**2

    return float(max(0, variance))


def map_scv(D0: np.ndarray, D1: np.ndarray = None) -> float:
    """
    Compute squared coefficient of variation (SCV) of a MAP.

    SCV = Var[X] / E[X]² = (E[X²] - E[X]²) / E[X]²

    Args:
        D0: Hidden transition matrix, or stacked [D0, D1] if D1 is None
        D1: Visible transition matrix (optional)

    Returns:
        Squared coefficient of variation
    """
    if D1 is None:
        D0_arr = np.asarray(D0)
        if D0_arr.ndim == 3 and D0_arr.shape[0] == 2:
            D0, D1 = D0_arr[0], D0_arr[1]
        else:
            raise ValueError("D0 must be a (2, n, n) array when D1 is not provided")

    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    # Match Kotlin: compute from moments directly
    e1 = map_moment(D0, D1, 1)
    e2 = map_moment(D0, D1, 2)

    if e1 > 0:
        var = e2 - e1 * e1
        scv = var / (e1 * e1)
        return max(0, scv)
    else:
        return 1.0


def map_moment(D0: np.ndarray, D1: np.ndarray, k: int) -> float:
    """
    Compute the k-th moment of inter-arrival time distribution.

    E[X^k] = k! * π_e * (-D0)^{-k} * e

    where π_e is the embedded DTMC steady-state.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        k: Moment order (k >= 1)

    Returns:
        k-th raw moment
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    if k < 1:
        raise ValueError("Moment order must be >= 1")

    # Check if D0 is zero matrix or singular (matching Kotlin implementation)
    if np.all(np.abs(D0) < 1e-14) or np.abs(linalg.det(D0)) < 1e-12:
        return 0.0

    n = D0.shape[0]
    e = np.ones(n)

    # Use embedded DTMC steady-state (map_pie), not CTMC steady-state
    pie = map_pie(D0, D1)

    try:
        D0_inv = linalg.inv(-D0)
    except LinAlgError:
        D0_inv = linalg.pinv(-D0)

    # Compute k! * (-D0)^{-k} incrementally as in Kotlin
    # Start with (-D0)^{-1}, then multiply by (-D0)^{-1} * i for factorial
    D0_inv_k = D0_inv.copy()
    for i in range(2, k + 1):
        D0_inv_k = D0_inv_k @ D0_inv
        D0_inv_k *= i  # Incorporate factorial incrementally

    # E[X^k] = k! * π_e * (-D0)^{-k} * e (factorial already incorporated)
    moment = pie @ D0_inv_k @ e

    return float(moment)


def map_scale(D0: np.ndarray, D1: np.ndarray, factor: float
              ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Scale a MAP by a given factor.

    Scaling changes the time scale: D0' = factor * D0, D1' = factor * D1.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        factor: Scaling factor (> 0)

    Returns:
        Tuple of scaled (D0', D1')
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    if factor <= 0:
        raise ValueError("Scaling factor must be positive")

    return factor * D0, factor * D1


def map_normalize(D0: np.ndarray, D1: np.ndarray
                  ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Normalize a MAP to have unit mean inter-arrival time.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix

    Returns:
        Tuple of normalized (D0', D1')
    """
    mean = map_mean(D0, D1)
    if mean > 0 and np.isfinite(mean):
        return map_scale(D0, D1, mean)
    else:
        return D0.copy(), D1.copy()


def map_isfeasible(D0: np.ndarray, D1: np.ndarray,
                   tolerance: float = 1e-10) -> bool:
    """
    Check if (D0, D1) form a valid MAP.

    A valid MAP requires:
    - D0 has non-positive diagonal and non-negative off-diagonal
    - D1 has non-negative elements
    - D0 + D1 is a valid generator (row sums = 0)

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        tolerance: Numerical tolerance

    Returns:
        True if valid MAP
    """
    D0 = np.asarray(D0)
    D1 = np.asarray(D1)

    n = D0.shape[0]
    if D0.shape != (n, n) or D1.shape != (n, n):
        return False

    # Check D0 diagonal non-positive
    if np.any(np.diag(D0) > tolerance):
        return False

    # Check D0 off-diagonal non-negative
    D0_offdiag = D0.copy()
    np.fill_diagonal(D0_offdiag, 0)
    if np.any(D0_offdiag < -tolerance):
        return False

    # Check D1 non-negative
    if np.any(D1 < -tolerance):
        return False

    # Check row sums of Q = D0 + D1 are zero
    Q = D0 + D1
    row_sums = Q.sum(axis=1)
    if np.any(np.abs(row_sums) > tolerance):
        return False

    return True


def exp_map(lambda_rate: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create a MAP representation of an exponential distribution.

    Args:
        lambda_rate: Rate parameter (λ > 0)

    Returns:
        Tuple of (D0, D1) matrices representing Exp(λ)
    """
    if lambda_rate <= 0:
        raise ValueError("Rate must be positive")

    D0 = np.array([[-lambda_rate]])
    D1 = np.array([[lambda_rate]])

    return D0, D1


def erlang_map(k: int, lambda_rate: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create a MAP representation of an Erlang-k distribution.

    Args:
        k: Number of phases (k >= 1)
        lambda_rate: Overall rate parameter

    Returns:
        Tuple of (D0, D1) matrices representing Erlang(k, k*λ)
    """
    if k < 1:
        raise ValueError("Number of phases must be >= 1")
    if lambda_rate <= 0:
        raise ValueError("Rate must be positive")

    # Phase rate = k * lambda to get mean = 1/lambda
    mu = k * lambda_rate

    # Build D0 and D1
    D0 = np.zeros((k, k))
    D1 = np.zeros((k, k))

    for i in range(k):
        D0[i, i] = -mu
        if i < k - 1:
            D0[i, i + 1] = mu

    # Arrival from last phase
    D1[k - 1, 0] = mu

    return D0, D1


def hyperexp_map(probs: np.ndarray, rates: np.ndarray
                 ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create a MAP representation of a hyperexponential distribution.

    Args:
        probs: Probability vector for choosing each phase
        rates: Rate parameters for each phase

    Returns:
        Tuple of (D0, D1) matrices representing hyperexponential
    """
    probs = np.asarray(probs, dtype=np.float64)
    rates = np.asarray(rates, dtype=np.float64)

    if len(probs) != len(rates):
        raise ValueError("probs and rates must have same length")

    k = len(probs)

    # D0: diagonal with -rates
    D0 = np.diag(-rates)

    # D1: arrivals return to initial state with probability probs
    D1 = np.zeros((k, k))
    for i in range(k):
        D1[i, :] = rates[i] * probs

    return D0, D1


# MATLAB-style wrapper functions for compatibility
def map_exponential(mean: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create a MAP representation of an exponential distribution (MATLAB-style).

    This is a MATLAB-compatible wrapper that takes mean instead of rate.

    Args:
        mean: Mean inter-arrival time (= 1/λ)

    Returns:
        Tuple of (D0, D1) matrices representing Exp(1/mean)

    Examples:
        >>> D0, D1 = map_exponential(2)  # Poisson process with rate λ=0.5
    """
    if mean <= 0:
        raise ValueError("Mean must be positive")

    mu = 1.0 / mean
    D0 = np.array([[-mu]])
    D1 = np.array([[mu]])
    return D0, D1


def map_erlang(mean: float, k: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create a MAP representation of an Erlang-k distribution (MATLAB-style).

    This is a MATLAB-compatible wrapper that takes mean as first argument.

    Args:
        mean: Mean inter-arrival time
        k: Number of phases

    Returns:
        Tuple of (D0, D1) matrices representing Erlang-k with given mean

    Examples:
        >>> D0, D1 = map_erlang(2, 3)  # Erlang-3 with mean 2
    """
    if mean <= 0:
        raise ValueError("Mean must be positive")
    if k < 1:
        raise ValueError("Number of phases must be >= 1")

    mu = k / mean

    # Build D0 and D1
    D0 = np.zeros((k, k))
    D1 = np.zeros((k, k))

    # D0: transitions between phases
    for i in range(k - 1):
        D0[i, i + 1] = mu

    # D1: arrival from last phase
    D1[k - 1, 0] = mu

    # Normalize D0 diagonal to make it a valid generator
    D0, D1 = _map_normalize_generator(D0, D1)

    return D0, D1


def _map_normalize_generator(D0: np.ndarray, D1: np.ndarray
                              ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Normalize MAP to have proper generator structure.

    Sets D0 diagonal so that (D0 + D1) has zero row sums.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix

    Returns:
        Tuple of normalized (D0, D1)
    """
    D0 = np.asarray(D0, dtype=np.float64).copy()
    D1 = np.asarray(D1, dtype=np.float64).copy()
    n = D0.shape[0]

    # Set non-negative off-diagonal elements
    for i in range(n):
        for j in range(n):
            D0[i, j] = np.real(D0[i, j])
            D1[i, j] = np.real(D1[i, j])
            if i != j and D0[i, j] < 0:
                D0[i, j] = 0
            if D1[i, j] < 0:
                D1[i, j] = 0

    # Normalize diagonal
    for i in range(n):
        D0[i, i] = 0
        row_sum = np.sum(D0[i, :]) + np.sum(D1[i, :])
        D0[i, i] = -row_sum

    return D0, D1


def map_sumind(maps: list) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the sum of independent MAPs.

    Creates a MAP representing the sum (concatenation) of independent
    random variables represented by the input MAPs.

    Args:
        maps: List of MAPs, each as (D0, D1) tuple

    Returns:
        Tuple of (D0, D1) matrices representing the sum

    Examples:
        >>> # Sum of exponential and Erlang-2
        >>> MAP1 = map_exponential(1.0)
        >>> MAP2 = map_erlang(1.0, 2)
        >>> D0, D1 = map_sumind([MAP1, MAP2])
    """
    n = len(maps)
    if n == 0:
        raise ValueError("At least one MAP required")
    if n == 1:
        return maps[0]

    # Get orders of each MAP
    orders = [map[0].shape[0] for map in maps]
    total_order = sum(orders)

    D0 = np.zeros((total_order, total_order))
    D1 = np.zeros((total_order, total_order))

    curpos = 0
    for i in range(n):
        D0_i, D1_i = maps[i][0], maps[i][1]
        order_i = orders[i]

        # Set diagonal block from D0_i
        D0[curpos:curpos+order_i, curpos:curpos+order_i] = D0_i

        if i < n - 1:
            # Transition to next MAP
            order_next = orders[i + 1]
            pie_next = map_pie(maps[i + 1][0], maps[i + 1][1])
            e_i = np.ones((order_i, 1))

            # D0 transition: D1_i * e * pie_next
            D0[curpos:curpos+order_i, curpos+order_i:curpos+order_i+order_next] = \
                D1_i @ e_i @ pie_next.reshape(1, -1)
        else:
            # Last MAP: transition back to first
            order_first = orders[0]
            pie_first = map_pie(maps[0][0], maps[0][1])
            e_i = np.ones((order_i, 1))

            # D1 transition back to first MAP
            D1[curpos:curpos+order_i, 0:order_first] = \
                D1_i @ e_i @ pie_first.reshape(1, -1)

        curpos += order_i

    return D0, D1


def map_cdf(D0: np.ndarray, D1: np.ndarray, points: np.ndarray) -> np.ndarray:
    """
    Compute cumulative distribution function of inter-arrival times.

    F(t) = 1 - π_e * exp(D0*t) * e

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        points: Time points at which to evaluate CDF

    Returns:
        CDF values at specified points

    Examples:
        >>> map_cdf(D0, D1, 1.0)  # Returns P(T <= 1)
        >>> map_cdf(D0, D1, [1.0, 5.0])  # Returns [P(T<=1), P(T<=5)]
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    points = np.asarray(points, dtype=np.float64).ravel()

    n = D0.shape[0]
    e = np.ones(n)
    pie = map_pie(D0, D1)

    cdf_vals = np.zeros(len(points))
    for i, t in enumerate(points):
        if t <= 0:
            cdf_vals[i] = 0.0
        else:
            exp_D0t = linalg.expm(D0 * t)
            cdf_vals[i] = 1.0 - pie @ exp_D0t @ e

    return cdf_vals


def map_pdf(D0: np.ndarray, D1: np.ndarray, points: np.ndarray) -> np.ndarray:
    """
    Compute probability density function of inter-arrival times.

    f(t) = π_e * exp(D0*t) * (-D0) * e

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        points: Time points at which to evaluate PDF

    Returns:
        PDF values at specified points

    Examples:
        >>> map_pdf(D0, D1, [0.5, 1.0, 2.0])
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    points = np.asarray(points, dtype=np.float64).ravel()

    n = D0.shape[0]
    e = np.ones(n)
    pie = map_pie(D0, D1)
    neg_D0 = -D0

    pdf_vals = np.zeros(len(points))
    for i, t in enumerate(points):
        if t < 0:
            pdf_vals[i] = 0.0
        else:
            exp_D0t = linalg.expm(D0 * t)
            pdf_vals[i] = pie @ exp_D0t @ neg_D0 @ e

    return pdf_vals


def map_hyperexp(probs: np.ndarray, means: np.ndarray
                 ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create a MAP representation of a hyperexponential distribution (MATLAB-style).

    This is a MATLAB-compatible wrapper that takes means instead of rates.

    Args:
        probs: Probability vector for choosing each phase
        means: Mean values for each phase (1/rates)

    Returns:
        Tuple of (D0, D1) matrices representing hyperexponential
    """
    probs = np.asarray(probs, dtype=np.float64)
    means = np.asarray(means, dtype=np.float64)
    rates = 1.0 / means
    return hyperexp_map(probs, rates)


def map_gamma(mean: float, n: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create a MAP approximation of a Gamma distribution.

    Uses Erlang-n with matching mean as an approximation.

    Args:
        mean: Mean of the distribution
        n: Number of phases (shape parameter)

    Returns:
        Tuple of (D0, D1) matrices
    """
    return map_erlang(mean, n)


def map_sample(D0: np.ndarray, D1: np.ndarray, n_samples: int,
               rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """
    Generate random samples from a MAP distribution.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        n_samples: Number of samples to generate
        rng: Optional random number generator

    Returns:
        Array of inter-arrival times
    """
    if rng is None:
        rng = np.random.default_rng()

    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    n = D0.shape[0]

    # Get initial state distribution
    pie = map_pie(D0, D1)

    samples = np.zeros(n_samples)

    for i in range(n_samples):
        # Choose initial state
        state = rng.choice(n, p=pie)
        t = 0.0

        while True:
            # Total rate out of current state
            rate_out = -D0[state, state]
            if rate_out <= 0:
                break

            # Sample sojourn time
            t += rng.exponential(1.0 / rate_out)

            # Transition probabilities
            trans_rates = np.zeros(2 * n)
            # Hidden transitions
            for j in range(n):
                if j != state:
                    trans_rates[j] = D0[state, j]
            # Visible transitions
            for j in range(n):
                trans_rates[n + j] = D1[state, j]

            trans_probs = trans_rates / np.sum(trans_rates)

            # Choose next transition
            next_trans = rng.choice(2 * n, p=trans_probs)

            if next_trans >= n:
                # Visible transition (arrival)
                samples[i] = t
                break
            else:
                # Hidden transition
                state = next_trans

    return samples


# ============================================================================
# Statistics Functions
# ============================================================================

def map_skew(D0: np.ndarray, D1: np.ndarray) -> float:
    """
    Compute skewness of inter-arrival times.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix

    Returns:
        Skewness of inter-arrival times
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    m1 = map_moment(D0, D1, 1)
    m2 = map_moment(D0, D1, 2)
    m3 = map_moment(D0, D1, 3)

    # M3 = E[X^3] - 3*E[X^2]*E[X] + 2*E[X]^3 (third central moment)
    M3 = m3 - 3 * m2 * m1 + 2 * m1 ** 3

    scv = map_scv(D0, D1)
    if scv > 0:
        skewness = M3 / (np.sqrt(scv) * m1) ** 3
    else:
        skewness = 0.0

    return float(skewness)


def map_kurt(D0: np.ndarray, D1: np.ndarray) -> float:
    """
    Compute kurtosis of inter-arrival times.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix

    Returns:
        Kurtosis of inter-arrival times
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    m1 = map_moment(D0, D1, 1)
    m2 = map_moment(D0, D1, 2)
    m3 = map_moment(D0, D1, 3)
    m4 = map_moment(D0, D1, 4)

    # Fourth central moment / variance^2
    var = map_var(D0, D1)
    if var > 0:
        kurt = (m4 - 4 * m3 * m1 + 6 * m2 * m1 ** 2 - 3 * m1 ** 4) / var ** 2
    else:
        kurt = 0.0

    return float(kurt)


def map_acf(D0: np.ndarray, D1: np.ndarray,
            lags: Union[int, np.ndarray] = 1) -> np.ndarray:
    """
    Compute autocorrelation coefficients of inter-arrival times.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        lags: Lag(s) at which to compute autocorrelation (default: 1)

    Returns:
        Array of autocorrelation coefficients at specified lags

    Examples:
        >>> map_acf(D0, D1)  # lag-1 autocorrelation
        >>> map_acf(D0, D1, np.arange(1, 11))  # first 10 autocorrelations
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    if isinstance(lags, (int, np.integer)):
        lags = np.array([lags])
    else:
        lags = np.asarray(lags, dtype=int)

    n = D0.shape[0]
    e = np.ones(n)

    P = map_embedded(D0, D1)
    lam = map_lambda(D0, D1)
    piq = map_piq(D0, D1)
    x = lam * piq

    try:
        invD0 = linalg.inv(-D0)
    except LinAlgError:
        invD0 = linalg.pinv(-D0)

    y = invD0 @ e

    acf_coeffs = np.zeros(len(lags))
    for i, lag in enumerate(lags):
        P_lag = np.linalg.matrix_power(P, lag)
        acf_coeffs[i] = x @ P_lag @ y

    scv = map_scv(D0, D1)
    if scv > 0:
        acf_coeffs = (acf_coeffs - 1) / scv
    else:
        acf_coeffs = np.zeros(len(lags))

    return acf_coeffs


def map_acfc(D0: np.ndarray, D1: np.ndarray,
             kset: np.ndarray, u: float) -> np.ndarray:
    """
    Compute autocorrelation of counting process at given lags.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        kset: Set of lags
        u: Length of timeslot (time scale)

    Returns:
        Autocorrelation coefficients at specified lags
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    kset = np.asarray(kset, dtype=int)

    n = D0.shape[0]
    Q = map_infgen(D0, D1)
    I = np.eye(n)
    e = np.ones(n)
    piq = map_piq(D0, D1)

    exp_Qu = linalg.expm(Q * u)
    I_minus_expQu = I - exp_Qu

    # PRE = piq * D1 * (I - exp(Q*u))
    PRE = piq @ D1 @ I_minus_expQu

    # (e * piq - Q)^(-1)
    try:
        inv_term = linalg.inv(np.outer(e, piq) - Q)
    except LinAlgError:
        inv_term = linalg.pinv(np.outer(e, piq) - Q)

    # POST = (I - exp(Q*u)) * (e*piq - Q)^(-2) * D1 * e
    inv_term_sq = inv_term @ inv_term
    POST = I_minus_expQu @ inv_term_sq @ D1 @ e

    vart = map_varcount(D0, D1, np.array([u]))[0]

    acfc = np.zeros(len(kset))
    for j, k in enumerate(kset):
        exp_Qku = linalg.expm(Q * (k - 1) * u)
        acfc[j] = PRE @ exp_Qku @ POST / vart if vart > 0 else 0.0

    return acfc


def map_idc(D0: np.ndarray, D1: np.ndarray) -> float:
    """
    Compute the asymptotic index of dispersion.

    I = SCV * (1 + 2 * sum_{k=1}^{inf} rho_k)

    where SCV is the squared coefficient of variation and rho_k is the
    lag-k autocorrelation coefficient.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix

    Returns:
        Asymptotic index of dispersion
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    n = D0.shape[0]
    e = np.ones(n)

    lam = map_lambda(D0, D1)
    pie = map_pie(D0, D1)
    piq = map_piq(D0, D1)
    Q = map_infgen(D0, D1)

    try:
        inv_term = linalg.inv(Q + np.outer(e, piq))
    except LinAlgError:
        inv_term = linalg.pinv(Q + np.outer(e, piq))

    I = 1 + 2 * (lam - pie @ inv_term @ D1 @ e)

    return float(I)


# ============================================================================
# Count Process Functions
# ============================================================================

def map_count_mean(D0: np.ndarray, D1: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    Compute mean of counting process at resolution t.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        t: Time period(s) for counting

    Returns:
        Mean arrivals in (0, t]
    """
    t = np.asarray(t, dtype=np.float64).ravel()
    lam = map_lambda(D0, D1)
    return lam * t


def map_count_var(D0: np.ndarray, D1: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    Compute variance of counting process at resolution t.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        t: Time period(s) for counting

    Returns:
        Variance of arrivals in (0, t]

    Reference:
        He and Neuts, "Markov chains with marked transitions", 1998
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64).ravel()

    n = D0.shape[0]
    D = D0 + D1
    I = np.eye(n)
    e = np.ones(n)

    theta = map_piq(D0, D1)

    try:
        tmp = linalg.inv(np.outer(e, theta) - D)
    except LinAlgError:
        tmp = linalg.pinv(np.outer(e, theta) - D)

    lam = theta @ D1 @ e
    c = theta @ D1 @ tmp
    d = tmp @ D1 @ e
    ll = theta @ D1 @ e

    v = np.zeros(len(t))
    for i, ti in enumerate(t):
        exp_Dt = linalg.expm(D * ti)
        v[i] = (ll - 2 * lam ** 2 + 2 * c @ D1 @ e) * ti - 2 * c @ (I - exp_Dt) @ d

    return v


def map_varcount(D0: np.ndarray, D1: np.ndarray, tset: np.ndarray) -> np.ndarray:
    """
    Compute variance of counting process (alternative implementation).

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        tset: Set of time points

    Returns:
        Variance at each time point
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    tset = np.asarray(tset, dtype=np.float64).ravel()

    n = D0.shape[0]
    Q = map_infgen(D0, D1)
    I = np.eye(n)
    e = np.ones(n)

    piq = map_piq(D0, D1)
    lam = 1.0 / map_mean(D0, D1)

    try:
        inv_term = linalg.inv(np.outer(e, piq) - Q)
    except LinAlgError:
        inv_term = linalg.pinv(np.outer(e, piq) - Q)

    PRE = lam - 2 * lam ** 2 + 2 * piq @ D1 @ inv_term @ D1 @ e
    inv_term_sq = inv_term @ inv_term

    varc = np.zeros(len(tset))
    for j, t in enumerate(tset):
        exp_Qt = linalg.expm(Q * t)
        POST = -2 * piq @ D1 @ (I - exp_Qt) @ inv_term_sq @ D1 @ e
        varc[j] = PRE * t + POST

    return varc


def map_count_moment(D0: np.ndarray, D1: np.ndarray, t: float,
                     orders: np.ndarray) -> np.ndarray:
    """
    Compute power moments of counts at resolution t.

    Uses numerical differentiation of the moment generating function.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        t: Resolution (time period)
        orders: Orders of moments to compute

    Returns:
        Power moments of counts
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    orders = np.asarray(orders, dtype=int).ravel()

    n = D0.shape[0]
    e = np.ones(n)
    theta = map_piq(D0, D1)

    def mgfunc(z):
        return theta @ linalg.expm(D0 * t + D1 * np.exp(z) * t) @ e

    # Numerical differentiation
    from scipy.misc import derivative

    M = np.zeros(len(orders))
    for i, order in enumerate(orders):
        M[i] = derivative(mgfunc, 0, n=order, dx=1e-6)

    return M


# ============================================================================
# Constructor Functions
# ============================================================================

def map_mmpp2(mean: float, scv: float, skew: float = -1,
              acf1: float = -1) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit an MMPP(2) as a MAP.

    Args:
        mean: Mean inter-arrival time
        scv: Squared coefficient of variation
        skew: Skewness (-1 for automatic minimization)
        acf1: Lag-1 autocorrelation (-1 for maximum feasible)

    Returns:
        Tuple of (D0, D1) matrices

    Examples:
        >>> D0, D1 = map_mmpp2(1, 2, -1, 0.2)  # Minimal skewness, ACF=0.2
        >>> D0, D1 = map_mmpp2(1, 2, -1, -1)   # Minimal skewness, max ACF
    """
    E1 = mean
    E2 = (1 + scv) * E1 ** 2
    E3 = -(2 * E1 ** 3 - 3 * E1 * E2 - skew * (E2 - E1 ** 2) ** (3 / 2))

    FEASTOL = map_feastol()
    if acf1 == -1:
        G2 = 1 - 10 * 10 ** (-FEASTOL)
    else:
        G2 = acf1 / (1 - 1 / scv) / 0.5 if scv != 1 else 0

    if skew == -1 and scv > 1:
        E3 = (3 / 2 + 0.001) * E2 ** 2 / E1

    scv_recomputed = (E2 - E1 ** 2) / E1 ** 2

    # Simplified MMPP(2) construction
    if G2 < 1e-6:
        # Renewal case (no autocorrelation)
        denom = 6 * E1 ** 3 * scv_recomputed + 3 * E1 ** 3 * scv_recomputed ** 2 + 3 * E1 ** 3 - 2 * E3
        if abs(denom) < 1e-14:
            denom = 1e-14
        mu00 = 2 * (6 * E1 ** 3 * scv_recomputed - E3) / E1 / denom
        mu11 = 0
        q01 = 9 * E1 ** 5 * (scv_recomputed - 1) * (scv_recomputed ** 2 - 2 * scv_recomputed + 1) / max(1e-14, (6 * E1 ** 3 * scv_recomputed - E3)) / denom
        q10 = -3 * (scv_recomputed - 1) * E1 ** 2 / max(1e-14, (6 * E1 ** 3 * scv_recomputed - E3))
    else:
        # Use hyperexponential for general case
        if scv_recomputed >= 1:
            # H2 approximation
            mu1 = E1 - 0.5 * np.sqrt(max(0, -4 * E1 ** 2 + 2 * E2))
            mu2 = E1 + 0.5 * np.sqrt(max(0, -4 * E1 ** 2 + 2 * E2))
            p = 0.5 - 0.5 * G2
            mu1 = max(1e-10, np.real(mu1))
            mu2 = max(1e-10, np.real(mu2))
            mu00 = 1 / mu1
            mu11 = 1 / mu2
            q01 = mu00 * p
            q10 = mu11 * (1 - p)
        else:
            # Erlang approximation for SCV < 1
            mu00 = 2 / E1
            mu11 = 2 / E1
            q01 = mu00 * 0.5
            q10 = mu11 * 0.5

    D0 = np.array([[-mu00 - q01, q01],
                   [q10, -mu11 - q10]])
    D1 = np.array([[mu00, 0],
                   [0, mu11]])

    # Normalize to ensure valid generator
    D0, D1 = _map_normalize_generator(D0, D1)

    return D0, D1


def map_gamma2(D0: np.ndarray, D1: np.ndarray) -> float:
    """
    Compute the second largest eigenvalue of embedded DTMC.

    This is the autocorrelation decay rate.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix

    Returns:
        Second largest eigenvalue (gamma_2)
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    P = map_embedded(D0, D1)
    eigvals = linalg.eigvals(P)

    # Sort by absolute value descending
    sorted_idx = np.argsort(np.abs(eigvals))[::-1]
    sorted_eigvals = eigvals[sorted_idx]

    # Return second largest
    if len(sorted_eigvals) > 1:
        return float(np.real(sorted_eigvals[1]))
    else:
        return 0.0


def map_rand(k: int = 2) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate a random MAP of order k.

    Args:
        k: Order of the MAP (default: 2)

    Returns:
        Tuple of (D0, D1) matrices
    """
    D0 = np.random.rand(k, k)
    D1 = np.random.rand(k, k)
    return _map_normalize_generator(D0, D1)


def map_randn(k: int, mu: Tuple[float, float] = (1.0, 1.0),
              sigma: Tuple[float, float] = (0.5, 0.5)
              ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate a random MAP with normally distributed elements.

    Args:
        k: Order of the MAP
        mu: Mean for (D0, D1) elements
        sigma: Standard deviation for (D0, D1) elements

    Returns:
        Tuple of (D0, D1) matrices
    """
    D0 = np.abs(np.random.normal(mu[0], sigma[0], (k, k)))
    D1 = np.abs(np.random.normal(mu[1], sigma[1], (k, k)))
    return _map_normalize_generator(D0, D1)


def map_renewal(D0: np.ndarray, D1: np.ndarray
                ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Remove all correlations from a MAP, creating a renewal process.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix

    Returns:
        Renewal MAP with same CDF but no correlations
    """
    D0 = np.asarray(D0, dtype=np.float64).copy()
    D1 = np.asarray(D1, dtype=np.float64).copy()

    n = D0.shape[0]
    e = np.ones((n, 1))
    pie = map_pie(D0, D1).reshape(1, -1)

    # D1_new = D1 * e * pie
    D1_new = D1 @ e @ pie

    return D0, D1_new


# ============================================================================
# Operation Functions
# ============================================================================

def map_embedded(D0: np.ndarray, D1: np.ndarray) -> np.ndarray:
    """
    Compute embedded discrete-time transition matrix.

    P = (-D0)^{-1} * D1

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix

    Returns:
        Transition matrix of embedded DTMC
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    try:
        P = linalg.inv(-D0) @ D1
    except LinAlgError:
        P = linalg.pinv(-D0) @ D1

    return P


def map_sum(D0: np.ndarray, D1: np.ndarray, n: int
            ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create MAP for sum of n IID random variables.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        n: Number of variables to sum

    Returns:
        Tuple of (D0_new, D1_new) for the sum
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    order = D0.shape[0]
    D0_new = np.zeros((n * order, n * order))
    D1_new = np.zeros((n * order, n * order))

    curpos = 0
    for i in range(n):
        D0_new[curpos:curpos + order, curpos:curpos + order] = D0
        if i < n - 1:
            D0_new[curpos:curpos + order,
                   curpos + order:curpos + 2 * order] = D1
        else:
            D1_new[curpos:curpos + order, 0:order] = D1
        curpos += order

    return D0_new, D1_new


def _krons(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Kronecker sum: A oplus B = A otimes I + I otimes B."""
    return np.kron(A, np.eye(B.shape[0])) + np.kron(np.eye(A.shape[0]), B)


def map_super(D0_a: np.ndarray, D1_a: np.ndarray,
              D0_b: np.ndarray, D1_b: np.ndarray
              ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create superposition of two MAPs.

    Args:
        D0_a, D1_a: First MAP
        D0_b, D1_b: Second MAP

    Returns:
        Superposed MAP (D0, D1)
    """
    D0_new = _krons(D0_a, D0_b)
    D1_new = _krons(D1_a, D1_b)
    return _map_normalize_generator(D0_new, D1_new)


def map_mixture(alpha: np.ndarray, maps: list
                ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create probabilistic mixture of MAPs.

    Args:
        alpha: Probability vector for choosing each MAP
        maps: List of MAPs as (D0, D1) tuples

    Returns:
        Mixture MAP (D0, D1)
    """
    alpha = np.asarray(alpha, dtype=np.float64)
    n_maps = len(maps)

    if len(alpha) != n_maps:
        raise ValueError("Alpha length must match number of MAPs")

    # Build block diagonal D0
    D0_blocks = [maps[i][0] for i in range(n_maps)]
    D0 = linalg.block_diag(*D0_blocks)

    # Build D1 with mixture transitions
    orders = [maps[i][0].shape[0] for i in range(n_maps)]
    total_order = sum(orders)
    D1 = np.zeros((total_order, total_order))

    row_offset = 0
    for i in range(n_maps):
        D1_i = maps[i][1]
        e_i = np.ones((orders[i], 1))

        col_offset = 0
        for j in range(n_maps):
            pie_j = map_pie(maps[j][0], maps[j][1]).reshape(1, -1)
            D1[row_offset:row_offset + orders[i],
               col_offset:col_offset + orders[j]] = \
                D1_i @ e_i * alpha[j] @ pie_j
            col_offset += orders[j]
        row_offset += orders[i]

    return _map_normalize_generator(D0, D1)


def map_max(D0_a: np.ndarray, D1_a: np.ndarray,
            D0_b: np.ndarray, D1_b: np.ndarray
            ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Create MAP for max(X, Y) where X ~ MAP_A and Y ~ MAP_B.

    Args:
        D0_a, D1_a: First MAP
        D0_b, D1_b: Second MAP

    Returns:
        MAP for the maximum (D0, D1)
    """
    D0_a = np.asarray(D0_a, dtype=np.float64)
    D1_a = np.asarray(D1_a, dtype=np.float64)
    D0_b = np.asarray(D0_b, dtype=np.float64)
    D1_b = np.asarray(D1_b, dtype=np.float64)

    na = D0_a.shape[0]
    nb = D0_b.shape[0]

    a = -D0_a @ np.ones((na, 1))
    b = -D0_b @ np.ones((nb, 1))

    # Build D0
    D0 = np.zeros((na * nb + na + nb, na * nb + na + nb))

    # Block (1,1): krons(A, B)
    D0[:na * nb, :na * nb] = _krons(D0_a, D0_b)
    # Block (1,2): kron(a, I_na)
    D0[:na * nb, na * nb:na * nb + na] = np.kron(a, np.eye(na))
    # Block (1,3): kron(b, I_nb)
    D0[:na * nb, na * nb + na:] = np.kron(b, np.eye(nb))
    # Block (2,2): D0_b
    D0[na * nb:na * nb + nb, na * nb:na * nb + nb] = D0_b
    # Block (3,3): D0_a
    D0[na * nb + nb:, na * nb + nb:] = D0_a

    # Build D1
    pie_a = map_pie(D0_a, D1_a)
    pie_b = map_pie(D0_b, D1_b)
    pie = np.concatenate([np.kron(pie_a, pie_b), np.zeros(na), np.zeros(nb)])
    D1 = np.outer(pie, np.ones(D0.shape[0]))

    return D0, D1


def map_timereverse(D0: np.ndarray, D1: np.ndarray
                    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute time-reversed MAP.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix

    Returns:
        Time-reversed MAP (D0_r, D1_r)
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    piq = map_piq(D0, D1)
    D = np.diag(piq)

    try:
        D_inv = linalg.inv(D)
    except LinAlgError:
        D_inv = linalg.pinv(D)

    D0_r = D_inv @ D0.T @ D
    D1_r = D_inv @ D1.T @ D

    return D0_r, D1_r


def map_mark(D0: np.ndarray, D1: np.ndarray, prob: np.ndarray
             ) -> list:
    """
    Mark arrivals from a MAP according to given probabilities.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        prob: Marking probabilities (prob[k] = P(mark type k))

    Returns:
        MMAP as list [D0, D1_class1, D1_class2, ...]
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    prob = np.asarray(prob, dtype=np.float64)

    if abs(np.sum(prob) - 1.0) > 1e-6:
        prob = prob / np.sum(prob)

    mmap = [D0]
    for k in range(len(prob)):
        mmap.append(prob[k] * D1)

    return mmap


def map_stochcomp(D0: np.ndarray, D1: np.ndarray,
                  retain_idx: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Stochastic complementation to reduce MAP order.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        retain_idx: Indices of states to retain (0-indexed)

    Returns:
        Reduced MAP (D0_new, D1_new)
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    retain_idx = np.asarray(retain_idx, dtype=int)

    Q = D0 + D1
    n = Q.shape[0]
    eliminated_idx = np.setdiff1d(np.arange(n), retain_idx)

    Q_RE = Q[np.ix_(retain_idx, eliminated_idx)]
    Q_EE = Q[np.ix_(eliminated_idx, eliminated_idx)]
    Q_RR = Q[np.ix_(retain_idx, retain_idx)]
    Q_ER = Q[np.ix_(eliminated_idx, retain_idx)]

    try:
        Q_EE_inv = linalg.inv(-Q_EE)
    except LinAlgError:
        Q_EE_inv = linalg.pinv(-Q_EE)

    Q_new = Q_RR + Q_RE @ Q_EE_inv @ Q_ER

    D1_RR = D1[np.ix_(retain_idx, retain_idx)]
    D1_ER = D1[np.ix_(eliminated_idx, retain_idx)]

    D0_new = Q_new - D1_RR
    D1_new = D1_RR + Q_RE @ Q_EE_inv @ D1_ER

    return _map_normalize_generator(D0_new, D1_new)


# ============================================================================
# Fitting Functions
# ============================================================================

def map_kpc(maps: list) -> Tuple[np.ndarray, np.ndarray]:
    """
    Kronecker product composition of MAPs.

    Args:
        maps: List of MAPs as (D0, D1) tuples

    Returns:
        Composed MAP (D0, D1)
    """
    if len(maps) == 0:
        raise ValueError("At least one MAP required")

    if len(maps) == 1:
        return maps[0]

    D0, D1 = maps[0]
    for i in range(1, len(maps)):
        D0_i, D1_i = maps[i]
        D0 = -np.kron(D0, D0_i)
        D1 = np.kron(D1, D1_i)

    return D0, D1


def map_bernstein(f, n: int = 20) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert distribution to MAP via Bernstein approximation.

    Args:
        f: PDF function handle
        n: Number of phases (default: 20)

    Returns:
        MAP representation (D0, D1)
    """
    # Bernstein approximation
    c = 0.0
    for i in range(1, n + 1):
        xi = -np.log(i / n)
        fi = f(xi)
        if np.isfinite(fi) and fi > 0:
            c += fi / i

    if c <= 0 or not np.isfinite(c):
        return map_erlang(1.0, n)

    # Build subgenerator T
    T = np.diag(-np.arange(1, n + 1)) + np.diag(np.arange(1, n), 1)

    # Build initial probability vector alpha
    alpha = np.zeros(n)
    for i in range(1, n + 1):
        xi = -np.log(i / n)
        fi = f(xi)
        if np.isfinite(fi) and fi > 0:
            alpha[i - 1] = fi / (i * c)

    if np.sum(alpha) > 0:
        alpha = alpha / np.sum(alpha)
    else:
        alpha[0] = 1

    # Convert to MAP
    P = np.tile(alpha, (n, 1))
    D0 = T
    D1 = -T @ P

    return D0, D1


def map_pntiter(D0: np.ndarray, D1: np.ndarray, na: int, t: float,
                M: Optional[int] = None) -> np.ndarray:
    """
    Compute probability of na arrivals in interval [0, t].

    Uses iterative bisection method (Neuts and Li).

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        na: Number of arrivals
        t: Time interval length
        M: Bisection depth (auto-computed if None)

    Returns:
        Matrix of probabilities
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    if M is None:
        mean = map_mean(D0, D1)
        M = max(0, int(np.ceil(np.log2(t * 100 / mean)))) if mean > 0 else 0

    def br(tau, t, r):
        """Poisson probability."""
        return np.exp(-tau * t) * (tau * t) ** r / np.math.factorial(r)

    def pnt_bisect(na, t):
        """Base case computation."""
        n_states = D0.shape[0]
        tau = np.max(-np.diag(D0))
        I = np.eye(n_states)
        K = D0 / tau + I
        K1 = D1 / tau

        # Find N for convergence
        epsilon = np.finfo(float).eps
        Nmax = 100
        N = 1
        for N in range(1, Nmax):
            S = sum(br(tau, t, nn) for nn in range(N + 1, Nmax))
            if S < epsilon:
                break

        # Initialize
        V = [[np.zeros_like(I) for _ in range(N + 1)] for _ in range(na + 1)]
        P = [np.zeros_like(I) for _ in range(na + 1)]

        V[0][0] = I
        P[0] = V[0][0] * br(tau, t, 0)

        for n in range(1, na + 1):
            for k in range(1, N + 1):
                V[n][k] = V[n][k - 1] @ K + V[n - 1][k - 1] @ K1
                P[n] = P[n] + V[n][k] * br(tau, t, n)

        return P

    if M <= 0:
        P = pnt_bisect(na, t)
    else:
        P = pnt_bisect(na, t / (2 ** M))
        for _ in range(M):
            Pold = P
            P = [np.zeros_like(Pold[0]) for _ in range(na + 1)]
            for n in range(na + 1):
                for j in range(n + 1):
                    P[n] = P[n] + Pold[j] @ Pold[n - j]

    return P[na]


def map_pntquad(D0: np.ndarray, D1: np.ndarray, na: int, t: float
                ) -> np.ndarray:
    """
    Compute probability of na arrivals using ODE solver.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        na: Number of arrivals
        t: Time interval length

    Returns:
        Matrix P_n(t)
    """
    from scipy.integrate import solve_ivp

    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    Ki = D0.shape[0]

    def pnt_ode(t_val, P):
        dP = np.zeros_like(P)
        # n = 0 case
        P0 = P[:Ki ** 2].reshape(Ki, Ki)
        dP[:Ki ** 2] = (P0 @ D0).flatten()

        # n >= 1 cases
        for n in range(1, na + 1):
            Pn_1 = P[((n - 1) * Ki ** 2):(n * Ki ** 2)].reshape(Ki, Ki)
            Pn = P[(n * Ki ** 2):((n + 1) * Ki ** 2)].reshape(Ki, Ki)
            dP[(n * Ki ** 2):((n + 1) * Ki ** 2)] = (Pn @ D0 + Pn_1 @ D1).flatten()

        return dP

    # Initial condition
    P0 = np.zeros((na + 1) * Ki ** 2)
    P0[:Ki ** 2] = np.eye(Ki).flatten()

    sol = solve_ivp(pnt_ode, [0, t], P0, method='RK45',
                    rtol=1e-10, atol=1e-10)

    P_final = sol.y[:, -1]
    Pnt = P_final[(na * Ki ** 2):((na + 1) * Ki ** 2)].reshape(Ki, Ki)

    return Pnt


# ============================================================================
# Utility Functions
# ============================================================================

def map_joint(D0: np.ndarray, D1: np.ndarray,
              a: np.ndarray, i: np.ndarray) -> float:
    """
    Compute joint moments of a MAP.

    E[(X_{a1})^{i1} * (X_{a1+a2})^{i2} * ...]

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        a: Vector of inter-arrival lags
        i: Vector of powers

    Returns:
        Joint moment
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    a = np.asarray(a, dtype=int)
    i = np.asarray(i, dtype=int)

    a_cum = np.cumsum(a)
    P = map_embedded(D0, D1)

    try:
        invD0 = linalg.inv(-D0)
    except LinAlgError:
        invD0 = linalg.pinv(-D0)

    K = len(a)
    JM = np.eye(D0.shape[0])

    for k in range(K - 1):
        P_power = np.linalg.matrix_power(P, a_cum[k + 1] - a_cum[k])
        invD0_power = np.linalg.matrix_power(invD0, i[k])
        JM = JM @ (np.math.factorial(i[k]) * invD0_power) @ P_power

    # Final term
    pie = map_pie(D0, D1)
    invD0_power_final = np.linalg.matrix_power(invD0, i[K - 1])
    e = np.ones(D0.shape[0])

    result = pie @ JM @ (np.math.factorial(i[K - 1]) * invD0_power_final) @ e

    return float(result)


def map_issym(D0: np.ndarray, D1: np.ndarray = None) -> bool:
    """
    Check if MAP contains symbolic elements.

    In Python, this always returns False as we use numpy arrays.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix (optional)

    Returns:
        False (Python arrays are always numeric)
    """
    return False


def map_feastol() -> int:
    """
    Get tolerance for feasibility checks.

    Returns:
        Tolerance exponent (10^(-feastol))
    """
    return 8


def map_feasblock(E1: float, E2: float, E3: float, G2: float,
                  opt: str = '') -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit the most similar feasible MAP(2).

    Args:
        E1: First moment (mean)
        E2: Second moment (or SCV if opt='scv')
        E3: Third moment
        G2: Autocorrelation decay rate
        opt: 'scv' if E2 is SCV instead of second moment

    Returns:
        Feasible MAP (D0, D1)
    """
    if opt.lower() == 'scv':
        E2 = (1 + E2) * E1 ** 2

    # Exponential case
    if abs(E2 - 2 * E1 ** 2) < 1e-10:
        D0 = np.array([[-1.0, 0], [0, -1.0]])
        D1 = np.array([[0.5, 0.5], [0.5, 0.5]])
        return map_scale(D0, D1, 1.0 / E1)

    KPC_TOL = 1e-10

    if E2 <= 2 * E1 ** 2:
        E2 = (1 + KPC_TOL) * E1 ** 2

    if E3 <= (3 / 2) * E2 ** 2 / E1:
        E3 = (3 / 2 + KPC_TOL) * E2 ** 2 / E1

    return map_block(E1, E2, E3, G2)


def map_block(E1: float, E2: float, E3: float, G2: float,
              opt: str = '') -> Tuple[np.ndarray, np.ndarray]:
    """
    Construct a MAP(2) from moments and autocorrelation.

    Args:
        E1: First moment
        E2: Second moment (or SCV if opt='scv')
        E3: Third moment
        G2: Autocorrelation decay rate
        opt: 'scv' if E2 is SCV

    Returns:
        MAP (D0, D1)
    """
    if opt.lower() == 'scv':
        E2 = (1 + E2) * E1 ** 2

    SCV = (E2 - E1 ** 2) / E1 ** 2

    if SCV >= 1:
        # Hyperexponential case
        mu1 = E1 - 0.5 * np.sqrt(max(0, 2 * E2 - 4 * E1 ** 2))
        mu2 = E1 + 0.5 * np.sqrt(max(0, 2 * E2 - 4 * E1 ** 2))
        mu1 = max(1e-10, np.real(mu1))
        mu2 = max(1e-10, np.real(mu2))
        p = 0.5 - 0.5 * G2

        D0 = np.array([[-1 / mu1, 0], [0, -1 / mu2]])
        D1 = np.array([[(1 - p) / mu1, p / mu1], [p / mu2, (1 - p) / mu2]])
    else:
        # Erlang-like case for SCV < 1
        mu = 2 / E1
        D0 = np.array([[-mu, mu * 0.5], [0, -mu]])
        D1 = np.array([[0, 0], [mu, 0]])

    return _map_normalize_generator(D0, D1)


def map_largemap() -> int:
    """
    Get threshold for "large" MAP where exact computation is expensive.

    Returns:
        Order threshold (default: 100)
    """
    return 100


def map2_fit(e1: float, e2: float, e3: float = -1.0, g2: float = 0.0
             ) -> Tuple[Optional[Tuple[np.ndarray, np.ndarray]], int]:
    """
    Fit a MAP(2) distribution to moments and autocorrelation.

    Based on: A. Heindl, G. Horvath, K. Gross "Explicit inverse characterization
    of acyclic MAPs of second order"

    Args:
        e1: First moment E[X]
        e2: Second moment E[X^2]
        e3: Third moment E[X^3], or -1 for automatic selection,
            -2 for minimum, -3 for maximum, or negative fraction for interpolation
        g2: Autocorrelation decay rate (gamma_2)

    Returns:
        Tuple of (MAP, error_code) where:
        - MAP: Tuple (D0, D1) if successful, None if failed
        - error_code: 0 for success, >0 for various errors

    Error codes:
        0: Success
        10: Mean out of bounds
        20: Correlated exponential
        30: h2 out of bounds
        40: h3 out of bounds
        51-54: g2 out of bounds
    """
    TOL = 1e-10

    r1 = e1
    r2 = e2 / 2
    h2 = (r2 - r1 ** 2) / r1 ** 2

    # Handle special cases for e3
    scv = (e2 - e1 ** 2) / e1 ** 2

    if e3 == -1:
        # Select e3 that maximizes the range of correlations
        if 1 <= scv < 3:
            if g2 < 0:
                h3 = h2 - h2 ** 2
                e3 = 12 * e1 ** 3 * h2 + 6 * e1 ** 3 * h3 + 6 * e1 ** 3 * (1 + h2 ** 2)
            else:
                e3 = (3 / 2 + 1e-3) * e2 ** 2 / e1
        elif scv >= 3:
            e3 = (3 / 2 + 1e-3) * e2 ** 2 / e1
        elif 0 < scv < 1:
            e3 = (1 + TOL) * (12 * e1 ** 3 * h2 + 6 * e1 ** 3 * (
                h2 * (1 - h2 - 2 * np.sqrt(-h2))) + 6 * e1 ** 3 * (1 + h2 ** 2))
    elif e3 == -2:
        # Select minimum e3
        if scv >= 1:
            e3 = (3 / 2 + 1e-6) * e2 ** 2 / e1
        elif 0 < scv < 1:
            h3 = h2 * (1 - h2 - 2 * np.sqrt(-h2))
            e3 = 6 * e1 ** 3 * (h2 ** 2 + h3)
    elif e3 == -3:
        # Select maximum e3
        if scv >= 1:
            e3 = 1e6
        elif 0 < scv < 1:
            h3 = (-h2) ** 2
            e3 = 6 * e1 ** 3 * (h2 ** 2 + h3)
    elif e3 == -4:
        # Select random e3
        r = np.random.rand()
        if scv >= 1:
            e3 = r * (3 / 2 + 1e-6) * e2 ** 2 / e1 + (1 - r) * 1e6
        elif 0 < scv < 1:
            h3 = r * (-h2) ** 2 + (1 - r) * h2 * (1 - h2 - 2 * np.sqrt(-h2))
            e3 = 6 * e1 ** 3 * (h2 ** 2 + h3)
    elif -1 < e3 < 0:
        # Use a custom random e3
        r = abs(e3)
        if scv >= 1:
            e3 = r * (3 / 2 + 1e-6) * e2 ** 2 / e1 + (1 - r) * 1e6
        elif 0 < scv < 1:
            h3 = r * h2 * (1 - h2 - 2 * np.sqrt(-h2)) + (1 - r) * (-h2) ** 2
            e3 = 6 * e1 ** 3 * (h2 ** 2 + h3)

    r3 = e3 / 6
    h3 = (r3 * r1 - r2 ** 2) / r1 ** 4
    b = h3 + h2 ** 2 - h2
    c = np.sqrt(b ** 2 + 4 * h2 ** 3 + 0j)  # Complex sqrt

    if r1 <= 0:
        return None, 10  # Mean out of bounds

    if h2 == 0:
        if h3 == 0 and g2 == 0:
            return map_exponential(e1), 0
        else:
            return None, 20  # Correlated exponential

    if h2 > 0 and h3 > 0:
        # Hyperexponential case
        if np.real(b) >= 0:
            bc_ratio = np.real(b / c) if np.abs(c) > TOL else 0
            lower_bound = (b - c) / (b + c) if np.abs(b + c) > TOL else -1
            if np.real(lower_bound) <= g2 < 1:
                D0 = (1 / (2 * r1 * h3)) * np.array([
                    [-(2 * h2 + b - c), 0],
                    [0, -(2 * h2 + b + c)]
                ])
                D1 = (1 / (4 * r1 * h3)) * np.array([
                    [(2 * h2 + b - c) * (1 - bc_ratio + g2 * (1 + bc_ratio)),
                     (2 * h2 + b - c) * (1 + bc_ratio) * (1 - g2)],
                    [(2 * h2 + b + c) * (1 - bc_ratio) * (1 - g2),
                     (2 * h2 + b + c) * (1 + bc_ratio + g2 * (1 - bc_ratio))]
                ])
                D0 = np.real(D0)
                D1 = np.real(D1)
                return (D0, D1), 0
            else:
                return None, 51  # g2 out of bounds
        else:  # b < 0
            if 0 <= g2 < 1:
                bc_ratio = np.real(b / c) if np.abs(c) > TOL else 0
                D0 = (1 / (2 * r1 * h3)) * np.array([
                    [-(2 * h2 + b - c), 0],
                    [0, -(2 * h2 + b + c)]
                ])
                D1 = (1 / (4 * r1 * h3)) * np.array([
                    [(2 * h2 + b - c) * (1 - bc_ratio + g2 * (1 + bc_ratio)),
                     (2 * h2 + b - c) * (1 + bc_ratio) * (1 - g2)],
                    [(2 * h2 + b + c) * (1 - bc_ratio) * (1 - g2),
                     (2 * h2 + b + c) * (1 + bc_ratio + g2 * (1 - bc_ratio))]
                ])
                D0 = np.real(D0)
                D1 = np.real(D1)
                return (D0, D1), 0
            elif -(h3 + h2 ** 2) / h2 <= g2 < 0:
                a = (h3 + h2 ** 2) / h2
                d1 = ((1 - a) * (2 * h2 * g2 + b - c) + g2 * (b + c) - (b - c)) / (
                    (1 - a) * (2 * h2 + b - c) + 2 * c)
                d2 = ((g2 - 1) * (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c)
                D0 = (1 / (2 * r1 * h3)) * np.array([
                    [-(2 * h2 + b - c), (2 * h2 + b - c) * (1 - a)],
                    [0, -(2 * h2 + b + c)]
                ])
                D1 = (1 / (2 * r1 * h3)) * np.array([
                    [(2 * h2 + b - c) * d1, (2 * h2 + b - c) * (a - d1)],
                    [(2 * h2 + b + c) * d2, (2 * h2 + b + c) * (1 - d2)]
                ])
                D0 = np.real(D0)
                D1 = np.real(D1)
                return (D0, D1), 0
            else:
                return None, 52  # g2 out of bounds

    elif -1 / 4 <= h2 < 0 and h2 * (1 - h2 - 2 * np.sqrt(-h2)) <= h3 <= -h2 ** 2:
        # Hypoexponential case
        c = -c  # Flip sign for hypo case
        if g2 >= 0:
            upper_bound = -(h2 + np.sqrt(-h3)) ** 2 / h2 if h2 != 0 else np.inf
            if g2 <= np.real(upper_bound):
                a = (2 * h2 + b - c) * (h2 + np.sqrt(-h3)) / (2 * h2 * np.sqrt(-h3))
                d1 = ((1 - a) * (2 * h2 * g2 + b - c) + g2 * (b + c) - (b - c)) / (
                    (1 - a) * (2 * h2 + b - c) + 2 * c)
                d2 = ((g2 - 1) * (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c)
                D0 = (1 / (2 * r1 * h3)) * np.array([
                    [-(2 * h2 + b - c), (2 * h2 + b - c) * (1 - a)],
                    [0, -(2 * h2 + b + c)]
                ])
                D1 = (1 / (2 * r1 * h3)) * np.array([
                    [(2 * h2 + b - c) * d1, (2 * h2 + b - c) * (a - d1)],
                    [(2 * h2 + b + c) * d2, (2 * h2 + b + c) * (1 - d2)]
                ])
                D0 = np.real(D0)
                D1 = np.real(D1)
                return (D0, D1), 0
            else:
                return None, 53  # g2 out of bounds
        else:  # g2 < 0
            lower_bound = -(h3 + h2 ** 2) / h2 if h2 != 0 else -np.inf
            if g2 >= np.real(lower_bound):
                a = (h3 + h2 ** 2) / h2
                d1 = ((1 - a) * (2 * h2 * g2 + b - c) + g2 * (b + c) - (b - c)) / (
                    (1 - a) * (2 * h2 + b - c) + 2 * c)
                d2 = ((g2 - 1) * (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c)
                D0 = (1 / (2 * r1 * h3)) * np.array([
                    [-(2 * h2 + b - c), (2 * h2 + b - c) * (1 - a)],
                    [0, -(2 * h2 + b + c)]
                ])
                D1 = (1 / (2 * r1 * h3)) * np.array([
                    [(2 * h2 + b - c) * d1, (2 * h2 + b - c) * (a - d1)],
                    [(2 * h2 + b + c) * d2, (2 * h2 + b + c) * (1 - d2)]
                ])
                D0 = np.real(D0)
                D1 = np.real(D1)
                return (D0, D1), 0
            else:
                return None, 54  # g2 out of bounds
    else:
        if not (-1 / 4 <= h2 < 0):
            return None, 30  # h2 out of bounds
        else:
            return None, 40  # h3 out of bounds


__all__ = [
    'map_infgen',
    'map_piq',
    'map_pie',
    'map_lambda',
    'map_mean',
    'map_var',
    'map_scv',
    'map_moment',
    'map_scale',
    'map_normalize',
    'map_isfeasible',
    # MAP constructors (rate-based API)
    'exp_map',
    'erlang_map',
    'hyperexp_map',
    # MAP constructors (MATLAB-style mean-based API)
    'map_exponential',
    'map_erlang',
    'map_hyperexp',
    'map_gamma',
    # MAP operations
    'map_sumind',
    'map_cdf',
    'map_pdf',
    'map_sample',
    # Statistics functions
    'map_skew',
    'map_kurt',
    'map_acf',
    'map_acfc',
    'map_idc',
    # Count process functions
    'map_count_mean',
    'map_count_var',
    'map_varcount',
    'map_count_moment',
    # Constructor functions
    'map_mmpp2',
    'map_gamma2',
    'map_rand',
    'map_randn',
    'map_renewal',
    # Operation functions
    'map_embedded',
    'map_sum',
    'map_super',
    'map_mixture',
    'map_max',
    'map_timereverse',
    'map_mark',
    'map_stochcomp',
    # Fitting functions
    'map_kpc',
    'map_bernstein',
    'map_pntiter',
    'map_pntquad',
    'map2_fit',
    # Utility functions
    'map_joint',
    'map_issym',
    'map_feastol',
    'map_feasblock',
    'map_block',
    'map_largemap',
]
