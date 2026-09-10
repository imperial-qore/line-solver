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

import math
import numpy as np
from numpy.linalg import LinAlgError
from scipy import linalg, optimize
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


def map_prob(D0: np.ndarray, D1: np.ndarray = None) -> np.ndarray:
    """Stationary distribution of the underlying CTMC of a MAP (alias of
    :func:`map_piq`, matching the MATLAB ``map_prob`` name)."""
    return map_piq(D0, D1)


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

    Uses map_moment for consistency with the JAR implementation.

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

    # Compute from moments directly
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

    # MATLAB map_moment.m guards exactly one degenerate case, `if MAP{1}==0`.
    # Do NOT reinstate a det(D0) test here: det scales as rate^n, so any fixed
    # threshold misfires on a perfectly healthy slow MAP -- order 4 with rates
    # around 1e-3, or order 16 with rates around 0.1, both land under 1e-12 --
    # and the moment comes back silently as 0.0. Singularity is scale-free only
    # through the condition number.
    if not np.any(D0):
        return 0.0
    if not np.isfinite(D0).all():
        return float('nan')
    with np.errstate(divide='ignore', invalid='ignore'):
        rcond = 1.0 / np.linalg.cond(D0)
    if not np.isfinite(rcond) or rcond < np.finfo(np.float64).eps:
        return 0.0

    n = D0.shape[0]
    e = np.ones(n)

    # Use embedded DTMC steady-state (map_pie), not CTMC steady-state
    pie = map_pie(D0, D1)

    try:
        D0_inv = linalg.inv(-D0)
    except LinAlgError:
        D0_inv = linalg.pinv(-D0)

    # Compute k! * (-D0)^{-k} incrementally
    # Start with (-D0)^{-1}, then multiply by (-D0)^{-1} * i for factorial
    D0_inv_k = D0_inv.copy()
    for i in range(2, k + 1):
        D0_inv_k = D0_inv_k @ D0_inv
        D0_inv_k *= i  # Incorporate factorial incrementally

    # E[X^k] = k! * π_e * (-D0)^{-k} * e (factorial already incorporated)
    moment = pie @ D0_inv_k @ e

    return float(moment)


def map_scale(D0: np.ndarray, D1: np.ndarray, new_mean: float
              ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Rescale a MAP to a given MEAN inter-arrival time.

    Port of map_scale.m: the third argument is the TARGET MEAN, not a
    multiplier. The rates are scaled by mean/new_mean, which leaves every
    normalized moment and every autocorrelation alone and moves only the first
    moment, and the result is passed through map_normalize (the feasibility
    repair) as the reference does.

    THIS ARGUMENT USED TO BE A FACTOR here, and nowhere else: MATLAB, the JAR
    and the C++ port all take the new mean, and this module's own private
    helper `_map_scale` in api/solvers/mam/mmap_fj.py already did too. The two
    conventions are each other's reciprocal-ish (a factor c gives mean/c), so
    a call written for one and read by the other produces a MAP with the wrong
    rate and the right shape, which no moment check on the SCV would catch.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        new_mean: target mean inter-arrival time (> 0)

    Returns:
        Tuple of rescaled (D0', D1')
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)

    if new_mean <= 0:
        raise ValueError("The target mean must be positive")

    ratio = map_mean(D0, D1) / new_mean
    return map_normalize(ratio * np.asarray(D0, dtype=np.float64),
                         ratio * np.asarray(D1, dtype=np.float64))


def map_normalize(D0: np.ndarray, D1: np.ndarray
                  ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Make a MAP feasible again: port of map_normalize.m.

    Takes real parts, clips negative entries to zero and rebuilds D0's diagonal
    so that (D0 + D1) e = 0. MATLAB, the JAR (`Map_normalize.java`) and the C++
    port (`map_transform.h`) all do exactly this, and every caller here wants
    it: `map_scale` closes with it, `kpcfit`'s rescaling helper closes with it,
    and `mmpp_rand` uses it to turn two random matrices into a generator pair.

    IT USED TO RESCALE THE MEAN TO ONE, which is a different operation
    altogether and left every one of those callers wrong in a way no moment
    check would show: `mmpp_rand` returned a D0 whose diagonal had never been
    repaired, so the pair was not a generator at all, and the two rescaling
    helpers had the mean they had just set pulled straight back to one. Use
    `map_scale(D0, D1, 1.0)` where unit mean is what is wanted.

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix

    Returns:
        Tuple of (D0', D1') satisfying the generator condition
    """
    A = np.real(np.asarray(D0, dtype=np.float64)).copy()
    B = np.real(np.asarray(D1, dtype=np.float64)).copy()
    A[A < 0] = 0.0
    B[B < 0] = 0.0
    for n in range(A.shape[0]):
        A[n, n] = 0.0
        A[n, n] = -(np.sum(A[n, :]) + np.sum(B[n, :]))
    return A, B


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
    The same repair as `map_normalize`, under the name the fitting code uses.

    Kept as one implementation rather than two: the pair had drifted, with the
    public name rescaling the mean and this one repairing the generator.
    """
    return map_normalize(D0, D1)


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


def map_hyperexp(probs, means=None, p: float = 0.99):
    """
    Two-phase hyperexponential process as a MAP.

    Two call forms, mirroring the two conventions in the codebase:

    - ``map_hyperexp(MEAN, SCV, p)`` with scalar arguments fits a two-phase
      hyperexponential to the given mean and squared coefficient of variation,
      selecting phase 1 with probability ``p`` (default 0.99). This is the
      MATLAB ``map_hyperexp.m`` signature. Returns ``None`` when the request is
      outside the feasible set, as MATLAB returns ``{}``: with a fixed ``p``
      the reachable SCV is bounded, e.g. SCV <= 3 at p = 0.5.
    - ``map_hyperexp(probs, means)`` with array arguments builds the MAP of a
      hyperexponential with the given phase probabilities and phase means.

    Returns:
        Tuple of (D0, D1) matrices, or None if the moment fit is infeasible
    """
    if means is not None and np.isscalar(probs) and np.isscalar(means):
        return _map_hyperexp_moments(float(probs), float(means), float(p))
    probs = np.asarray(probs, dtype=np.float64)
    means = np.asarray(means, dtype=np.float64)
    rates = 1.0 / means
    return hyperexp_map(probs, rates)


def _map_hyperexp_moments(MEAN: float, SCV: float, p: float):
    """Moment fit of map_hyperexp.m, including its fallback to a smaller p."""
    E2 = (1 + SCV) * MEAN ** 2
    Delta = -4 * p * MEAN ** 2 + 4 * p ** 2 * MEAN ** 2 + 2 * E2 * p - 2 * E2 * p ** 2
    if Delta < 0:
        return None
    root = np.sqrt(Delta)
    for sign in (1.0, -1.0):
        denom = E2 * p - 2 * MEAN ** 2
        if denom == 0:
            continue
        mu2 = (-2 * MEAN + 2 * p * MEAN + sign * root) / denom
        den1 = p - 1 + MEAN * mu2
        if den1 == 0:
            continue
        mu1 = mu2 * p / den1
        D0 = np.array([[-mu1, 0.0], [0.0, -mu2]])
        D1 = np.array([[mu1 * p, mu1 * (1 - p)], [mu2 * p, mu2 * (1 - p)]])
        if map_isfeasible(D0, D1):
            return (D0, D1)
    if p > 1e-6:
        return _map_hyperexp_moments(MEAN, SCV, p / 10.0)
    return None


def map_gamma(D0: np.ndarray, D1: np.ndarray, limit: int = 1000) -> float:
    """
    Estimate the autocorrelation decay rate of a MAP.

    Mirrors MATLAB map_gamma. For MAPs of order higher than 2 the ACF is not
    geometric, so the decay rate is obtained by fitting rho_k = RHO0*gamma^k in
    the least-squares sense, with RHO0 = (1 - 1/SCV)/2 held fixed.

    This is NOT a Gamma-distribution constructor; use map_erlang for that. It is
    also distinct from map_gamma2, which returns the second largest eigenvalue
    of the embedded DTMC (the two agree for order 2 only).

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        limit: Maximum lag considered when fitting (default: 1000)

    Returns:
        Autocorrelation decay rate
    """
    D0 = np.asarray(D0, dtype=np.float64)
    D1 = np.asarray(D1, dtype=np.float64)
    n = D0.shape[0]

    if n == 1:
        # Poisson process: no correlation
        return 0.0

    if n == 2:
        # second-order MAP: geometric ACF
        acf1 = float(np.ravel(map_acf(D0, D1, 1))[0])
        if abs(acf1) < 1e-8:
            # phase-type
            return 0.0
        return float(np.ravel(map_acf(D0, D1, 2))[0]) / acf1

    # higher-order MAP: non-geometric, fit the ACF curve
    lag = np.arange(1, limit + 1, max(1, limit // 10))
    m1 = map_mean(D0, D1)
    m2 = map_moment(D0, D1, 2)
    scv = (m2 - m1 ** 2) / m1 ** 2
    rho0 = 0.5 * (1.0 - 1.0 / scv)
    rho = np.ravel(map_acf(D0, D1, lag))

    return _fit_geometric_acf(lag, rho, rho0)


def _fit_geometric_acf(lag: np.ndarray, rho: np.ndarray, rho0: float,
                       start: float = 0.99) -> float:
    """
    Fit rho_k = rho0*gamma^k for gamma by robust nonlinear least squares.

    Mirrors MATLAB nlinfit with RobustWgtFun='fair': an ordinary least-squares
    fit, then iteratively reweighted fits with the fair weight w=1/(1+|r|),
    tuning constant 1.4, residuals adjusted by the leverage of the
    least-squares Jacobian and scaled by a MAD estimate of sigma.

    Args:
        lag: Lags at which rho was evaluated
        rho: Autocorrelation at those lags
        rho0: Fixed lag-0 coefficient (1 - 1/SCV)/2
        start: Initial guess for gamma

    Returns:
        The fitted decay rate gamma
    """
    lag = np.asarray(lag, dtype=np.float64)
    rho = np.asarray(rho, dtype=np.float64)

    def model(gamma):
        return rho0 * np.power(gamma, lag)

    def jacobian(gamma):
        return (rho0 * lag * np.power(gamma, lag - 1.0)).reshape(-1, 1)

    def lsq(x0, weights=None):
        if weights is None:
            fun = lambda p: model(p[0]) - rho
            jac = lambda p: jacobian(p[0])
        else:
            sw = np.sqrt(weights)
            fun = lambda p: sw * (model(p[0]) - rho)
            jac = lambda p: sw.reshape(-1, 1) * jacobian(p[0])
        fit = optimize.least_squares(fun, np.array([x0]), jac=jac, method='lm',
                                     xtol=1e-14, ftol=1e-14, gtol=1e-14,
                                     max_nfev=100000)
        return float(fit.x[0])

    gamma = lsq(start)

    # Leverage of the least-squares Jacobian, as advised by DuMouchel & O'Brien.
    # Held fixed across the reweighting, as nlinfit does.
    q, _ = np.linalg.qr(jacobian(gamma))
    h = np.minimum(0.9999, np.sum(q * q, axis=1))
    adjust = 1.0 / np.sqrt(1.0 - h)

    # A near-perfect fit would drive the MAD estimate of sigma to zero and make
    # every point an outlier, so floor it against the spread of the response
    tiny_s = 1e-6 * float(np.std(rho, ddof=1))
    if tiny_s == 0.0:
        tiny_s = 1.0

    TUNE = 1.4  # fair
    delta = np.sqrt(np.finfo(np.float64).eps)
    for _ in range(200):
        previous = gamma
        radj = (rho - model(gamma)) * adjust
        # one parameter, so no residual is dropped from the MAD
        sigma = float(np.median(np.abs(radj))) / 0.6745
        weights = 1.0 / (1.0 + np.abs(radj / (max(sigma, tiny_s) * TUNE)))
        gamma = lsq(previous, weights)
        if abs(gamma - previous) < delta * max(abs(gamma), abs(previous)):
            break

    return gamma


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


def map_count_idc(D0: np.ndarray, D1: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    Index of dispersion for counts (IDC) of a MAP at time point(s) t.

    The IDC of the counting process A(t) associated to the MAP is
    I_a(t) = Var(A(t)) / E[A(t)], t > 0, i.e. the scaled variance-time curve.
    It interpolates between I_a(0+) = SCV of the interarrival time (renewal MAP)
    and the asymptotic value I_a(inf) = map_idc(MAP).

    Reference:
        W. Whitt and W. You, "A Robust Queueing Network Analyzer Based on
        Indices of Dispersion", eq. (1).

    Args:
        D0: Hidden transition matrix
        D1: Visible transition matrix
        t: Time period(s)

    Returns:
        Column vector of IDC values, one per element of t.
    """
    t = np.asarray(t, dtype=np.float64).ravel()
    m = map_count_mean(D0, D1, t)
    v = map_count_var(D0, D1, t)
    I = np.ones(len(t))
    nz = m > 0
    I[nz] = v[nz] / m[nz]
    # For an orderly point process the counting IDC tends to 1 as t->0
    # (locally Poisson: Var(N(t)) ~ E[N(t)]). Only hit at t==0 exactly.
    return I


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
        return float(theta @ linalg.expm(D0 * t + D1 * np.exp(z) * t) @ e)

    # see _kb/03-api-layer.md for rationale
    def _nth_derivative(f, x0, n, h=1e-2):
        def stencil(hh):
            if n == 1:
                return (f(x0 + hh) - f(x0 - hh)) / (2 * hh)
            if n == 2:
                return (f(x0 + hh) - 2 * f(x0) + f(x0 - hh)) / hh ** 2
            if n == 3:
                return (f(x0 + 2 * hh) - 2 * f(x0 + hh) + 2 * f(x0 - hh)
                        - f(x0 - 2 * hh)) / (2 * hh ** 3)
            if n == 4:
                return (f(x0 + 2 * hh) - 4 * f(x0 + hh) + 6 * f(x0)
                        - 4 * f(x0 - hh) + f(x0 - 2 * hh)) / hh ** 4
            raise ValueError("map_count_moment supports orders 1..4 numerically")
        d1 = stencil(h)
        d2 = stencil(h / 2.0)
        return (4.0 * d2 - d1) / 3.0

    M = np.zeros(len(orders))
    for i, order in enumerate(orders):
        M[i] = _nth_derivative(mgfunc, 0.0, int(order))

    return M


# ============================================================================
# Constructor Functions
# ============================================================================

def map_mmpp2(mean: float, scv: float, skew: float = -1,
              acf1: float = -1) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fit an MMPP(2) as a MAP.

    Matches the requested mean, SCV, skewness and lag-1 autocorrelation
    exactly. Raises ValueError when the request lies outside the MMPP(2)
    feasible set rather than returning a non-MAP.

    Args:
        mean: Mean inter-arrival time
        scv: Squared coefficient of variation (>= 1)
        skew: Skewness (-1 for automatic minimization)
        acf1: Lag-1 autocorrelation (-1 for maximum feasible)

    Returns:
        Tuple of (D0, D1) matrices

    Examples:
        >>> D0, D1 = map_mmpp2(1, 2, -1, 0.2)  # Minimal skewness, ACF=0.2
        >>> D0, D1 = map_mmpp2(1, 2, -1, -1)   # Minimal skewness, max ACF
    """
    # see _kb/03-api-layer.md for rationale
    from line_solver.lib.kpctoolbox.map import map_mmpp2 as _kpc_map_mmpp2
    D0, D1 = _kpc_map_mmpp2(mean, scv, skew, acf1)
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

    The phase space is ordered as [(i,j) pairs, B-only phases, A-only phases]:
    in the first block both A and B are still running, in the second block A has
    already completed and B is awaited, in the third block B has completed and A
    is awaited. An arrival is recorded when the second of the two completes,
    i.e. only out of the last two blocks.

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
    npair = na * nb

    # completion-rate vectors of each process out of each of its phases
    a = -D0_a @ np.ones((na, 1))
    b = -D0_b @ np.ones((nb, 1))

    # Build D0
    D0 = np.zeros((npair + nb + na, npair + nb + na))

    # Block (1,1): both still running; state (i,j) has index i*nb+j
    D0[:npair, :npair] = _krons(D0_a, D0_b)
    # Block (1,2): A completes out of phase i, B stays in phase j -> B-only j
    D0[:npair, npair:npair + nb] = np.kron(a, np.eye(nb))
    # Block (1,3): B completes out of phase j, A stays in phase i -> A-only i
    D0[:npair, npair + nb:] = np.kron(np.eye(na), b)
    # Block (2,2): B evolves on its own until it completes
    D0[npair:npair + nb, npair:npair + nb] = D0_b
    # Block (3,3): A evolves on its own until it completes
    D0[npair + nb:, npair + nb:] = D0_a

    # Build D1: both processes restart from their embedded equilibrium
    pie_a = map_pie(D0_a, D1_a)
    pie_b = map_pie(D0_b, D1_b)
    pie = np.concatenate([np.kron(pie_a, pie_b), np.zeros(nb), np.zeros(na)])
    # rate at which the second of the two processes completes; it is zero in the
    # pair block, where only the first of the two can still complete
    d = np.concatenate([np.zeros(npair), b.ravel(), a.ravel()])
    D1 = np.outer(d, pie)

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
        return np.exp(-tau * t) * (tau * t) ** r / math.factorial(r)

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

        # Uniformization: V(n,k) is the sub-stochastic matrix of paths making
        # exactly n arrivals in k STEPS of the uniformized chain, so the Poisson
        # weight mixing the terms is that of the STEP count k, not of the
        # arrival count n, and V(0,k) = V(0,k-1) @ K carries real mass for every
        # k rather than being zero past k = 0.
        #   V(0,0) = I,  V(0,k) = V(0,k-1) @ K,
        #   V(n,k) = V(n,k-1) @ K + V(n-1,k-1) @ K1,
        #   P_n(t) = sum_k br(tau,t,k) * V(n,k)
        #
        # Both faults were previously present and NEITHER IS VISIBLE ON A
        # POISSON PROCESS: there K = D0/tau + I = 0, V(n,k) collapses to
        # delta(n,k), and the two weights coincide on the only surviving term.
        # The identities that do see them are P_0(t) = expm(D0 t) and
        # sum_n P_n(t) = expm((D0 + D1) t).
        V[0][0] = I
        P[0] = V[0][0] * br(tau, t, 0)
        for k in range(1, N + 1):
            V[0][k] = V[0][k - 1] @ K
            P[0] = P[0] + V[0][k] * br(tau, t, k)

        for n in range(1, na + 1):
            for k in range(1, N + 1):
                V[n][k] = V[n][k - 1] @ K + V[n - 1][k - 1] @ K1
                P[n] = P[n] + V[n][k] * br(tau, t, k)

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
        JM = JM @ (math.factorial(i[k]) * invD0_power) @ P_power

    # Final term
    pie = map_pie(D0, D1)
    invD0_power_final = np.linalg.matrix_power(invD0, i[K - 1])
    e = np.ones(D0.shape[0])

    result = pie @ JM @ (math.factorial(i[K - 1]) * invD0_power_final) @ e

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
    Get the feasibility tolerance exponent for MAPs.

    This is the exponent k of the toolbox feasibility tolerance 10^-k, so the
    tolerance itself is 10^-8. It is NOT the tolerance: map_feastol() == 8.

    Returns:
        Tolerance exponent k, to be used as 10**(-map_feastol())
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
        return map_scale(D0, D1, E1)

    KPC_TOL = 1e-10

    if E2 <= 2 * E1 ** 2:
        E2 = (2 + KPC_TOL) * E1 ** 2

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


def map2_fit_idc(e1: float, e2: float, e3: float, idc: float
                 ) -> Tuple[Tuple[np.ndarray, np.ndarray], int]:
    """
    Fit a MAP(2) matching the first three moments and the index of dispersion.

    A MAP(2) has a geometrically decaying autocorrelation, so its index of
    dispersion obeys I = SCV + (SCV-1)*g2/(1-g2), as reported in Section 5.2.2 of
    Casale, Mi, Cherkasova and Smirni, IEEE Trans. Soft. Eng. 37(5), 2011. The
    relation is inverted in closed form as g2 = (I-SCV)/(I-1) and the decay rate is
    passed to map2_fit. A third moment outside the feasible region is replaced by
    its lower limit (3/2)*e2^2/e1, the largest heavy-tail decay a MAP(2) admits.

    The paper returns an exponential whenever SCV <= 1 or I < SCV. The rule does
    more than avoid an infeasible fit and must not be relaxed: a flow-equivalent
    server whose service is exponential and load dependent is exact for a
    product-form subnetwork by Norton's theorem, whereas any MAP(2) fitted to the
    marginal inter-departure statistics is not, because the departure stream of the
    subnetwork is not independent of the rest of the model.

    Args:
        e1: first moment
        e2: second moment
        e3: third moment
        idc: asymptotic index of dispersion

    Returns:
        Tuple (MAP, status) with status 0 all four descriptors matched,
        1 exponential as burstiness is not representable, 2 third moment clamped,
        3 third moment selected automatically, 4 fit failed and an exponential is
        returned
    """
    scv = (e2 - e1 ** 2) / e1 ** 2

    if scv <= 1 + 1e-8 or idc < scv:
        return map_exponential(e1), 1

    g2 = (idc - scv) / (idc - 1)

    fit, err = map2_fit(e1, e2, e3, g2)
    if err == 0:
        return fit, 0

    e3min = (3.0 / 2 + 1e-6) * e2 ** 2 / e1
    if e3 < e3min:
        fit, err = map2_fit(e1, e2, e3min, g2)
        if err == 0:
            return fit, 2

    fit, err = map2_fit(e1, e2, -1, g2)
    if err == 0:
        return fit, 3

    return map_exponential(e1), 4


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
