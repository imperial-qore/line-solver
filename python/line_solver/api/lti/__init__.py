"""
LTI: Laplace Transform Inversion algorithms.

Native Python implementations of numerical Laplace transform inversion methods:
- Euler method (based on Abate-Whitt acceleration)
- Talbot method (contour integration)
- Gaver-Stehfest method

These methods are used for computing distributions and performance metrics
from their Laplace transforms in queueing theory applications.
"""

import numpy as np
from numpy.typing import NDArray
from typing import Callable, Optional, Tuple, List
from scipy.special import comb
from math import factorial, log, tan, pi, exp


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
                         n: int = 99) -> float:
    """
    Invert Laplace transform using Euler method.

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


def laplace_invert(F: Callable, t: float, method: str = 'euler',
                   n: Optional[int] = None) -> float:
    """
    Invert Laplace transform using specified method.

    Args:
        F: Laplace transform function F(s)
        t: Time point to evaluate at
        method: 'euler', 'talbot', or 'gaver-stehfest'
        n: Number of terms (default depends on method)

    Returns:
        Approximate value of f(t)
    """
    if method == 'euler':
        n = n or 99
        return laplace_invert_euler(F, t, n)
    elif method == 'talbot':
        n = n or 32
        return laplace_invert_talbot(F, t, n)
    elif method == 'gaver-stehfest' or method == 'gaver_stehfest':
        n = n or 12
        return laplace_invert_gaver_stehfest(F, t, n)
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

    for i, t in enumerate(t_values):
        if t <= 0:
            result[i] = 0.0
        else:
            result[i] = laplace_invert(F, t, method, n)

    # Ensure non-negative
    result = np.maximum(result, 0)

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
    'laplace_invert',
    'laplace_invert_cdf',
    'laplace_invert_pdf',
]
