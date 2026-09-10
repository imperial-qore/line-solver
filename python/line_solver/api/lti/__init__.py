"""
LTI: Laplace Transform Inversion algorithms.

Native Python implementations of numerical Laplace transform inversion methods:
- CME method (Concentrated Matrix Exponential, based on Abate-Whitt with pre-computed coefficients)
- Euler method (based on Abate-Whitt acceleration)
- Talbot method (contour integration)
- Gaver-Stehfest method

These methods are used for computing distributions and performance metrics
from their Laplace transforms in queueing theory applications.
"""

import json
import os
import numpy as np
from numpy.typing import NDArray
from typing import Callable, Optional, Tuple, List
from scipy.special import comb
from math import factorial, log, tan, pi, exp, floor

def euler_get_alpha(n: int) -> np.ndarray:
    """
    Get alpha coefficients for Euler method.

    Args:
        n: Number of terms (should be odd)

    Returns:
        Complex array of alpha values
    """
    result = np.zeros(n, dtype=complex)
    for i in range(n):
        result[i] = complex((n - 1) * log(10.0) / 6, pi * i)
    return result

def euler_get_eta(n: int) -> np.ndarray:
    """
    Get eta coefficients for Euler method.

    Args:
        n: Number of terms (should be odd)

    Returns:
        Array of eta values
    """
    res = np.zeros(n)
    res[0] = 0.5

    # Euler defined only for odd n
    for i in range(1, (n + 1) // 2):
        res[i] = 1.0

    res[n - 1] = 1.0 / (2.0 ** ((n - 1) / 2.0))

    for i in range(1, (n - 1) // 2):
        res[n - i - 1] = res[n - i] + (2.0 ** ((1 - n) / 2.0)) * comb((n - 1) // 2, i, exact=True)

    return res

def euler_get_omega(n: int) -> np.ndarray:
    """
    Get omega coefficients for Euler method.

    Args:
        n: Number of terms (should be odd)

    Returns:
        Complex array of omega values
    """
    eta = euler_get_eta(n)

    res = np.zeros(n, dtype=complex)

    for i in range(1, n + 1):
        res[i - 1] = (10.0 ** ((n - 1) / 6.0)) * ((-1.0) ** (i - 1)) * eta[i - 1]

    return res

def talbot_get_alpha(n: int) -> np.ndarray:
    """
    Get alpha coefficients for Talbot method.

    Args:
        n: Number of terms

    Returns:
        Complex array of alpha values
    """
    arr = np.zeros(n, dtype=complex)

    # For k = 1
    arr[0] = complex(2.0 * n / 5.0, 0.0)

    # For k = 2 onwards
    for i in range(2, n + 1):
        theta = (i - 1) * pi / n
        real_part = 2 * (i - 1) * pi / 5 * (1.0 / tan(theta))
        imag_part = 2 * (i - 1) * pi / 5
        arr[i - 1] = complex(real_part, imag_part)

    return arr

def talbot_get_omega(n: int, alpha: np.ndarray) -> np.ndarray:
    """
    Get omega coefficients for Talbot method.

    Args:
        n: Number of terms
        alpha: Alpha coefficients from talbot_get_alpha

    Returns:
        Complex array of omega values
    """
    arr = np.zeros(n, dtype=complex)

    # For k = 1
    arr[0] = np.exp(alpha[0]) / 5.0

    # For k = 2 onwards
    for i in range(2, n + 1):
        theta = (i - 1) * pi / n
        current_alpha_exp = np.exp(alpha[i - 1])
        tan_theta = tan(theta)
        cot_theta = 1.0 / tan_theta
        multiplier = complex(1.0, theta * (1 + cot_theta**2) - cot_theta)
        arr[i - 1] = 2 * current_alpha_exp / 5.0 * multiplier

    return arr

def gaver_stehfest_get_omega(n: int) -> np.ndarray:
    """
    Get omega coefficients for Gaver-Stehfest method.

    Args:
        n: Number of terms (will be rounded down to even)

    Returns:
        Array of omega values
    """
    if n % 2 == 1:
        n = n - 1  # Gaver-Stehfest only supports even n

    res = np.zeros(n)

    for k in range(1, n + 1):
        val = ((-1.0) ** (n // 2 + k)) * log(2.0)

        # Summation
        sum_val = 0.0
        j = int((k + 1) // 2)
        while j <= min(k, n // 2):
            val2 = (j ** (n // 2 + 1))
            val2 /= factorial(n // 2)
            val2 *= comb(n // 2, j, exact=True)
            val2 *= comb(2 * j, j, exact=True)
            val2 *= comb(j, k - j, exact=True)
            sum_val += val2
            j += 1

        val *= sum_val
        res[k - 1] = val

    return res

def gaver_stehfest_get_alpha(n: int) -> np.ndarray:
    """
    Get alpha coefficients for Gaver-Stehfest method.

    Args:
        n: Number of terms (will be rounded down to even)

    Returns:
        Array of alpha values (real)
    """
    if n % 2 == 1:
        n = n - 1

    res = np.zeros(n)
    for k in range(1, n + 1):
        res[k - 1] = k * log(2.0)

    return res

def laplace_invert_euler(F: Callable[[complex], complex], t: float,
                         n: int = 41) -> float:
    """
    Invert Laplace transform using Euler method.

    The default is 41, odd so the rounding rule below leaves it alone. THE OLD
    DEFAULT OF 99 WAS UNUSABLE IN DOUBLE PRECISION: the weights carry a factor
    10**((n-1)/6) against an ALTERNATING sum, so accuracy is a race between the
    series converging and the cancellation eating the mantissa. Worst relative
    error on F(s) = 2/(s+2), whose inverse is 2*exp(-2t), over
    t in {0.1, 0.5, 1, 2}:

        n   =    11      21      31      41      51      71      99
        err = 4.4e-3  2.1e-6  1.6e-9  1.6e-10 4.7e-8  1.7e-4  1.4e+0

    At 99 the scale is 2.2e16, past what a double resolves, and the result is
    140 per cent wrong and negative at some t. Pass a larger n only in higher
    precision; in double, 41 is the optimum.

    Args:
        F: Laplace transform function F(s)
        t: Time point to evaluate at
        n: Number of terms (default: 99, must be odd)

    Returns:
        Approximate value of f(t)
    """
    if n % 2 == 0:
        n = n + 1  # Ensure odd

    alpha = euler_get_alpha(n)
    omega = euler_get_omega(n)

    # Evaluate F at all s points
    result = 0.0
    for i in range(n):
        s = alpha[i] / t
        result += (omega[i] * F(s)).real

    return result / t

def laplace_invert_talbot(F: Callable[[complex], complex], t: float,
                          n: int = 32) -> float:
    """
    Invert Laplace transform using Talbot method.

    Args:
        F: Laplace transform function F(s)
        t: Time point to evaluate at
        n: Number of terms (default: 32)

    Returns:
        Approximate value of f(t)
    """
    alpha = talbot_get_alpha(n)
    omega = talbot_get_omega(n, alpha)

    result = 0.0
    for i in range(n):
        s = alpha[i] / t
        result += (omega[i] * F(s)).real

    return result / t

def laplace_invert_gaver_stehfest(F: Callable[[float], float], t: float,
                                   n: int = 12) -> float:
    """
    Invert Laplace transform using Gaver-Stehfest method.

    This method only works for real-valued Laplace transforms.

    Args:
        F: Laplace transform function F(s) (real)
        t: Time point to evaluate at
        n: Number of terms (default: 12, will be rounded to even)

    Returns:
        Approximate value of f(t)
    """
    if n % 2 == 1:
        n = n - 1

    alpha = gaver_stehfest_get_alpha(n)
    omega = gaver_stehfest_get_omega(n)

    result = 0.0
    for i in range(n):
        s = alpha[i] / t
        result += omega[i] * F(s)

    return result / t

def laplace_weeks_coeffs(F: Callable[[complex], complex], sigma: float = 0.0,
                         b: float = 1.0, p0: int = 200) -> np.ndarray:
    """
    Laguerre coefficients q_n, n = 0..2*p0-1, of the damped and scaled function
    f_{sigma,b}(t) = exp(-sigma t) f(t/b), whose Laguerre generating function is

        Q_{sigma,b}(z) = b/(1-z) * L( b(1+z)/(2(1-z)) + b*sigma )

    (Harrison and Knottenbelt 2002, Sec. 4.2) and q_n = the n-th Laurent
    coefficient of Q on the circle |z| = r.

    NOTE ON THE PAPER. Eq. 10 as printed carries the factor (1-z) rather than
    1/(1-z). The scaled form above, printed later in the same section, carries
    1/(1-z) and is the correct one: with l_n(t) = exp(-t/2) L_n(t) the
    transform of l_n is (s-1/2)^n/(s+1/2)^{n+1}, so L(s) = Q(z)/(s+1/2) with
    z = (s-1/2)/(s+1/2) and s+1/2 = 1/(1-z). Implementing the printed (1-z)
    is wrong at every t (163 per cent at t = 0.1 on Exp(2)).

    Sec. 4.3 fixes the trapezoid count at 2*p0 and the radius at
    r = 0.1^(4/p0) for every n instead of letting them grow with n. The
    quadrature is then a DFT of Q sampled on the circle, so one FFT yields all
    the coefficients and the transform is evaluated 2*p0 times in total.
    """
    if b <= 0:
        raise ValueError("The scaling parameter b must be positive.")
    N = 2 * p0
    r = 0.1 ** (4.0 / p0)
    j = np.arange(N)
    z = r * np.exp(2j * np.pi * j / N)
    s = b * (1.0 + z) / (2.0 * (1.0 - z)) + b * sigma
    Qz = np.array([b / (1.0 - z[k]) * F(s[k]) for k in range(N)])
    q = np.fft.fft(Qz).real / N
    return q / (r ** j)

def laplace_weeks_scaling(F: Callable[[complex], complex], p0: int = 200,
                          tol: float = 1e-10) -> Tuple[float, float, np.ndarray]:
    """
    Automatic damping sigma and scaling b for the Laguerre inversion, following
    the algorithm of Fig. 1 of Harrison and Knottenbelt 2002: accept the first
    (sigma, b) at which the coefficients have decayed by term p0, doubling
    sigma from 0.001, and stepping b by 4 whenever sigma passes 0.2.

    Raising b too far is counterproductive and excessive damping is unstable in
    finite precision, so the search is bounded. When it exhausts the box the
    failure is REFUSED BY NAME: a density with a discontinuity in itself or its
    derivatives has no usable Laguerre representation (Sec. 4.2), and returning
    the last iterate would report noise as an answer. Use 'euler' for those.
    """
    sigma, b = 0.0, 1.0
    while True:
        q = laplace_weeks_coeffs(F, sigma, b, p0)
        if abs(q[p0]) <= tol and abs(q[p0 + 1]) <= tol:
            return sigma, b, q
        sigma = 0.001 if sigma == 0.0 else 2.0 * sigma
        if sigma > 0.2:
            b = b + 4.0
            if b > 10.0:
                raise ValueError(
                    "no suitable scaling parameters were found for the Laguerre "
                    "inversion: the transform's density is not smooth enough for "
                    "a Laguerre series. Use laplace_invert(F, t, 'euler') instead.")
            sigma = 0.0

def _weeks_nterms(q: np.ndarray) -> int:
    """
    Truncate at the FIRST index where the coefficients have decayed, never the
    last. The quadrature divides by r^n with r < 1, so past the genuine decay
    the entries are rounding noise amplified by r^-n: at n = 2*p0 that factor is
    1e8, and scanning for the last entry above a threshold sums 1e-8 of pure
    noise (worst error on Exp(2) 2.7e-09 instead of 1.9e-14).
    """
    p0 = len(q) // 2
    for n in range(1, p0 - 1):
        if abs(q[n]) <= 1e-13 and abs(q[n + 1]) <= 1e-13:
            return n
    return p0

def _laguerre_functions(t: float, N: int) -> np.ndarray:
    """l_n(t) = exp(-t/2) L_n(t) by the stable recursion of Sec. 4.1."""
    l = np.zeros(N)
    l[0] = np.exp(-t / 2.0)
    if N > 1:
        l[1] = (1.0 - t) * l[0]
    for n in range(2, N):
        l[n] = ((2 * n - 1 - t) / n) * l[n - 1] - ((n - 1.0) / n) * l[n - 2]
    return l

def laplace_invert_weeks(F: Callable[[complex], complex], t,
                         q: Optional[np.ndarray] = None,
                         sigma: Optional[float] = None,
                         b: Optional[float] = None):
    """
    Invert a Laplace transform by the Laguerre (Weeks) series
    f(t) = sum_n q_n l_n(t), recovered as exp(sigma*b*t) f_{sigma,b}(b*t).

    Unlike Euler and Talbot the coefficients do not depend on t, so ONE
    coefficient set serves an arbitrary number of time points: the transform is
    evaluated 2*p0 times in total, not 2*p0 times per t. That is the property
    this method is here for, so pass a q from laplace_weeks_scaling when
    inverting on a grid rather than letting each call rebuild it.

    Reference: Weeks, JACM 13, 1966; Abate, Choudhury and Whitt, INFORMS J.
    Computing 8(4), 1996; Harrison and Knottenbelt 2002, Sec. 4.1-4.3.
    """
    if q is None:
        sigma, b, q = laplace_weeks_scaling(F)
    elif sigma is None or b is None:
        raise ValueError(
            "When q is supplied, sigma and b must be supplied with it: they are "
            "the parameters q was built at, and evaluating q at other values "
            "silently returns a different function.")
    nterms = _weeks_nterms(q)
    scalar = np.isscalar(t)
    tv = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.zeros(len(tv))
    for i, tval in enumerate(tv):
        if tval <= 0:
            continue
        l = _laguerre_functions(b * tval, nterms)
        out[i] = np.exp(sigma * b * tval) * float(np.dot(q[:nterms], l))
    return float(out[0]) if scalar else out

def laplace_invert(F: Callable, t: float, method: str = 'euler',
                   n: Optional[int] = None) -> float:
    """
    Invert Laplace transform using specified method.

    Args:
        F: Laplace transform function F(s)
        t: Time point to evaluate at
        method: 'euler', 'talbot', 'gaver-stehfest', 'cme' or 'weeks'
        n: Number of terms (default depends on method: 41 euler, 32 talbot,
           12 gaver-stehfest, 25 cme, 200 weeks -- the p0 of its search)

    Returns:
        Approximate value of f(t)
    """
    if method == 'weeks':
        n = n or 200
        sigma, b, q = laplace_weeks_scaling(F, n)
        return laplace_invert_weeks(F, t, q, sigma, b)
    if method == 'euler':
        n = n or 41
        return laplace_invert_euler(F, t, n)
    elif method == 'talbot':
        n = n or 32
        return laplace_invert_talbot(F, t, n)
    elif method == 'gaver-stehfest' or method == 'gaver_stehfest':
        n = n or 12
        return laplace_invert_gaver_stehfest(F, t, n)
    elif method == 'cme':
        n = n or 25
        return laplace_invert_cme(F, t, n)
    else:
        raise ValueError(f"Unknown method: {method}")

def laplace_invert_cdf(F: Callable[[complex], complex], t_values: np.ndarray,
                       method: str = 'euler', n: Optional[int] = None
                      ) -> np.ndarray:
    """
    Invert Laplace transform of a CDF at multiple time points.

    Args:
        F: Laplace transform function F(s) of PDF (not CDF)
        t_values: Array of time points
        method: Inversion method
        n: Number of terms

    Returns:
        Array of CDF values
    """
    t_values = np.asarray(t_values)
    result = np.zeros_like(t_values, dtype=float)

    # For CDF, use F(s)/s where F(s) is the Laplace transform of PDF
    def F_cdf(s):
        if abs(s) < 1e-15:
            return 1.0
        return F(s) / s

    if method == 'weeks':
        # One coefficient set serves the whole grid; this is the point of Weeks.
        sigma, b, q = laplace_weeks_scaling(F_cdf, n or 200)
        result = laplace_invert_weeks(F_cdf, t_values, q, sigma, b)
        return np.maximum.accumulate(np.clip(result, 0, 1))

    for i, t in enumerate(t_values):
        if t <= 0:
            result[i] = 0.0
        else:
            result[i] = laplace_invert(F_cdf, t, method, n)

    # Ensure CDF properties
    result = np.clip(result, 0, 1)
    result = np.maximum.accumulate(result)  # Ensure monotonicity

    return result

def laplace_invert_pdf(F: Callable[[complex], complex], t_values: np.ndarray,
                       method: str = 'euler', n: Optional[int] = None
                      ) -> np.ndarray:
    """
    Invert Laplace transform of a PDF at multiple time points.

    Args:
        F: Laplace transform function F(s)
        t_values: Array of time points
        method: Inversion method
        n: Number of terms

    Returns:
        Array of PDF values
    """
    t_values = np.asarray(t_values)
    result = np.zeros_like(t_values, dtype=float)

    if method == 'weeks':
        sigma, b, q = laplace_weeks_scaling(F, n or 200)
        return np.maximum(laplace_invert_weeks(F, t_values, q, sigma, b), 0)

    for i, t in enumerate(t_values):
        if t <= 0:
            result[i] = 0.0
        else:
            result[i] = laplace_invert(F, t, method, n)

    # Ensure non-negative
    result = np.maximum(result, 0)

    return result

# CME parameter cache (lazy-loaded singleton)
_cme_params_cache = None

def iltcme_load_params():
    """
    Load CME parameters from iltcme.json.

    Returns a list of dicts with keys: n, a, b, c, omega, mu1, mu2, cv2, etc.
    Results are cached after first load.
    """
    global _cme_params_cache
    if _cme_params_cache is None:
        json_path = os.path.join(os.path.dirname(__file__), 'iltcme.json')
        with open(json_path) as f:
            _cme_params_cache = json.load(f)
    return _cme_params_cache

def laplace_invert_cme(F: Callable[[complex], complex], t: float,
                        maxFnEvals: int = 25) -> float:
    """
    Invert Laplace transform using CME (Concentrated Matrix Exponential) method.

    Uses pre-computed CME parameters from iltcme.json for the Abate-Whitt framework.

    Args:
        F: Laplace transform function F(s)
        t: Time point to evaluate at (must be positive)
        maxFnEvals: Maximum number of function evaluations allowed

    Returns:
        Approximate value of f(t)
    """
    params_list = iltcme_load_params()

    # Find the most steep CME satisfying maxFnEvals
    best = params_list[0]
    for p in params_list:
        if p['cv2'] < best['cv2'] and p['n'] + 1 <= maxFnEvals:
            best = p

    n = best['n']
    mu1 = best['mu1']
    a = best['a']
    b = best['b']
    c = best['c']
    omega = best['omega']

    # eta = [c*mu1, (a + i*b)*mu1]
    eta = np.concatenate(([c], np.array(a) + 1j * np.array(b))) * mu1
    # beta = [1, 1 + i*(1:n)*omega] * mu1
    beta = np.concatenate(([1.0], 1.0 + 1j * np.arange(1, n + 1) * omega)) * mu1

    # Abate-Whitt evaluation
    result = 0.0
    for j in range(len(eta)):
        s = beta[j] / t
        result += (eta[j] * F(s)).real
    return result / t

def ilt(F: Callable, T, maxFnEvals: int, method: str = 'cme'):
    """
    Unified inverse Laplace transform using Abate-Whitt framework.

    Matches the MATLAB matlab_ilt.m function signature and behavior exactly.
    Supports CME, Euler, and Gaver-Stehfest methods.

    Args:
        F: Laplace transform function F(s) (must accept complex arguments)
        T: Time point(s) to evaluate at (scalar or array-like, must be positive)
        maxFnEvals: Maximum number of function evaluations allowed
        method: 'cme' (default), 'euler', or 'gaver'

    Returns:
        numpy array of f(t) values
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))

    if method == 'cme':
        params_list = iltcme_load_params()
        best = params_list[0]
        for p in params_list:
            if p['cv2'] < best['cv2'] and p['n'] + 1 <= maxFnEvals:
                best = p
        n = best['n']
        mu1 = best['mu1']
        eta = np.concatenate(([best['c']], np.array(best['a']) + 1j * np.array(best['b']))) * mu1
        beta = np.concatenate(([1.0], 1.0 + 1j * np.arange(1, n + 1) * best['omega'])) * mu1

    elif method == 'euler':
        n_euler = int(floor((maxFnEvals - 1) / 2))
        eta = np.concatenate(([0.5], np.ones(n_euler), np.zeros(n_euler - 1), [2 ** -n_euler]))
        for k in range(1, n_euler):
            log_binom = sum(log(i) for i in range(1, n_euler + 1)) \
                        - n_euler * log(2) \
                        - sum(log(i) for i in range(1, k + 1)) \
                        - sum(log(i) for i in range(1, n_euler - k + 1))
            eta[2 * n_euler - k] = eta[2 * n_euler - k + 1] + exp(log_binom)
        k = np.arange(2 * n_euler + 1)
        beta = n_euler * log(10) / 3 + 1j * pi * k
        eta = (10 ** (n_euler / 3)) * (1 - (k % 2) * 2) * eta

    elif method == 'gaver':
        mfe = maxFnEvals
        if mfe % 2 == 1:
            mfe -= 1
        ndiv2 = mfe // 2
        eta = np.zeros(mfe)
        beta = np.zeros(mfe)
        logsum = np.concatenate(([0], np.cumsum(np.log(np.arange(1, mfe + 1)))))
        for k in range(1, mfe + 1):
            inside_sum = 0.0
            for j in range(int(floor((k + 1) / 2)), min(k, ndiv2) + 1):
                inside_sum += exp(
                    (ndiv2 + 1) * log(j)
                    - logsum[ndiv2 - j]
                    + logsum[2 * j]
                    - 2 * logsum[j]
                    - logsum[k - j]
                    - logsum[2 * j - k]
                )
            eta[k - 1] = log(2) * ((-1) ** (k + ndiv2)) * inside_sum
            beta[k - 1] = k * log(2)
    else:
        raise ValueError(
            f"Unknown inverse Laplace transform method: {method}. Supported: cme, euler, gaver"
        )

    # Common Abate-Whitt evaluation: f(t) = (1/t) * sum(real(eta * F(beta/t)))
    # Use np.dot for numerical stability with large alternating-sign terms (Gaver)
    result = np.zeros(len(T))
    for i, t_val in enumerate(T):
        f_vals = np.array([F(b / t_val) for b in beta])
        result[i] = np.dot(eta, f_vals).real / t_val

    return result

__all__ = [
    'euler_get_alpha',
    'euler_get_eta',
    'euler_get_omega',
    'talbot_get_alpha',
    'talbot_get_omega',
    'gaver_stehfest_get_alpha',
    'gaver_stehfest_get_omega',
    'laplace_invert_euler',
    'laplace_invert_talbot',
    'laplace_invert_gaver_stehfest',
    'laplace_weeks_coeffs',
    'laplace_weeks_scaling',
    'laplace_invert_weeks',
    'laplace_invert',
    'laplace_invert_cdf',
    'laplace_invert_pdf',
    'iltcme_load_params',
    'laplace_invert_cme',
    'ilt',
]
