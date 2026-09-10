"""
Hypothesis tests used by the output-analysis procedures.

Both are general purpose and usable on their own. They are implemented here
rather than taken from ``scipy.stats`` so that the three LINE codebases run the
same algorithm; ``shapirowilk`` agrees with ``scipy.stats.shapiro`` to 5e-10 in
the statistic and 1.5e-7 in the p-value, which is the check that validates the
port.

References:
    Original MATLAB: matlab/src/api/sim/sim_vonneumann.m, sim_shapirowilk.m
    J. von Neumann, "Distribution of the Ratio of the Mean Square Successive
    Difference to the Variance", Ann. Math. Statist. 12(4), 1941.
    L. C. Young, "Randomness in Ordered Sequences", Ann. Math. Statist. 12, 1941.
    J. P. Royston, "Approximating the Shapiro-Wilk W-test for Non-normality",
    Statistics and Computing 2, 1992; "Remark AS R94", Applied Statistics 44(4),
    1995.
"""

from math import asin, exp, log, pi, sqrt
from typing import Dict, Sequence

import numpy as np

from .dist import normcdf, norminv

__all__ = ['vonneumann', 'shapirowilk', 'shapirowilk_weights']

# Royston's correction polynomials for the two extreme weights
_C1 = (0.0, 0.221157, -0.147981, -2.071190, 4.434685, -2.706056)
_C2 = (0.0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633)


def _poly(coefficients: Sequence[float], x: float) -> float:
    value = 0.0
    power = 1.0
    for c in coefficients:
        value += c * power
        power *= x
    return value


def vonneumann(x: Sequence[float], alpha: float = 0.05) -> Dict[str, object]:
    """
    Von Neumann ratio test for randomness of a sequence.

    The statistic is the ratio of the mean square successive difference to the
    variance,

        ratio = sum_{i=1}^{b-1} (x_{i+1}-x_i)^2 / sum_{i=1}^{b} (x_i - xbar)^2.

    Under the null hypothesis that the sequence is i.i.d. normal the ratio has
    mean 2 and variance ``4(b-2)/((b-1)(b+1))``, and the standardized statistic
    is asymptotically normal, giving a two-sided p-value. Serial correlation of
    either sign moves the ratio away from 2: positive correlation shrinks the
    successive differences and pushes the ratio below 2, negative correlation
    pushes it above. Those null moments were confirmed by Monte Carlo over
    b = 10, 16, 24, 32, 50 to within 0.3%.

    Args:
        x: The sequence, at least 3 finite and not all equal values
        alpha: Significance level in (0,1)

    Returns:
        Dict with keys ``ratio``, ``zscore``, ``pvalue``, ``reject``, ``nobs``
        and ``analyzer``.

    Raises:
        ValueError: If the input is too short, not finite, or constant.
    """
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0,1), got %r" % (alpha,))
    v = np.asarray(x, dtype=float).ravel()
    b = v.size
    if b < 3:
        raise ValueError("At least 3 observations are required, got %d" % b)
    if not np.all(np.isfinite(v)):
        raise ValueError("The sequence must be finite")

    den = float(np.sum((v - v.mean()) ** 2))
    if den <= 0.0:
        raise ValueError("The sequence is constant, the ratio is undefined")

    ratio = float(np.sum(np.diff(v) ** 2)) / den
    sd = sqrt(4.0 * (b - 2) / ((b - 1) * (b + 1)))
    zscore = (ratio - 2.0) / sd
    pvalue = 2.0 * (1.0 - normcdf(abs(zscore)))

    return {'ratio': ratio, 'zscore': zscore, 'pvalue': pvalue,
            'reject': bool(pvalue < alpha), 'nobs': int(b),
            'analyzer': 'vonneumann'}


def shapirowilk_weights(n: int) -> np.ndarray:
    """
    Royston AS R94 antisymmetric weight vector, ``a[n-1-i] == -a[i]``.

    Args:
        n: Sample size, at least 3

    Returns:
        The weight vector, ascending with the order statistics.
    """
    if n < 3:
        raise ValueError("At least 3 observations are required, got %d" % n)
    if n == 3:
        return np.array([-sqrt(0.5), 0.0, sqrt(0.5)])

    m = np.array([norminv((i - 0.375) / (n + 0.25)) for i in range(1, n + 1)])
    mm = float(m @ m)
    c = m / sqrt(mm)
    u = 1.0 / sqrt(n)

    a = m.copy()
    an = c[n - 1] + _poly(_C1, u)
    if n > 5:
        anm1 = c[n - 2] + _poly(_C2, u)
        phi = ((mm - 2.0 * m[n - 1] ** 2 - 2.0 * m[n - 2] ** 2)
               / (1.0 - 2.0 * an ** 2 - 2.0 * anm1 ** 2))
        a[2:n - 2] = m[2:n - 2] / sqrt(phi)
        a[n - 1] = an
        a[n - 2] = anm1
        a[0] = -an
        a[1] = -anm1
    else:
        phi = (mm - 2.0 * m[n - 1] ** 2) / (1.0 - 2.0 * an ** 2)
        a[1:n - 1] = m[1:n - 1] / sqrt(phi)
        a[n - 1] = an
        a[0] = -an
    return a


def shapirowilk(x: Sequence[float], alpha: float = 0.05) -> Dict[str, object]:
    """
    Shapiro-Wilk test for univariate normality, Royston's AS R94 algorithm.

    The statistic is

        W = (sum_i a_i x_(i))^2 / sum_i (x_i - xbar)^2,

    with ``x_(i)`` the order statistics and ``a`` the antisymmetric weight
    vector obtained by correcting the normalized expected normal order
    statistics ``m_i = Phi^{-1}((i-3/8)/(n+1/4))`` in their two extreme
    components. Small W means departure from normality, so the test is one-sided
    in W and the p-value is an upper normal tail after Royston's normalizing
    transform, which has three branches: n = 3 exact, 4 <= n <= 11, and n >= 12.
    Valid for 3 <= n <= 5000.

    Args:
        x: The sample, 3 to 5000 finite and not all equal values
        alpha: Significance level in (0,1)

    Returns:
        Dict with keys ``W``, ``pvalue``, ``zscore``, ``reject``, ``nobs`` and
        ``analyzer``. ``zscore`` is NaN for n = 3, where the exact null law is
        used instead of the transform.

    Raises:
        ValueError: If the input is out of range, not finite, or constant.
    """
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0,1), got %r" % (alpha,))
    v = np.sort(np.asarray(x, dtype=float).ravel())
    n = v.size
    if n < 3:
        raise ValueError("At least 3 observations are required, got %d" % n)
    if n > 5000:
        raise ValueError("The AS R94 approximation is valid up to n = 5000, got %d" % n)
    if not np.all(np.isfinite(v)):
        raise ValueError("The sample must be finite")

    ssd = float(np.sum((v - v.mean()) ** 2))
    if ssd <= 0.0:
        raise ValueError("The sample is constant, W is undefined")

    a = shapirowilk_weights(n)
    W = min(float(a @ v) ** 2 / ssd, 1.0)

    if n == 3:
        # exact null law, W is supported on [3/4, 1]
        pvalue = 6.0 / pi * (asin(sqrt(W)) - asin(sqrt(0.75)))
        pvalue = min(max(pvalue, 0.0), 1.0)
        zscore = float('nan')
    else:
        if n <= 11:
            g = -2.273 + 0.459 * n
            w = -log(g - log(1.0 - W))
            mu = 0.5440 - 0.39978 * n + 0.025054 * n ** 2 - 0.0006714 * n ** 3
            sigma = exp(1.3822 - 0.77857 * n + 0.062767 * n ** 2 - 0.0020322 * n ** 3)
        else:
            ln = log(n)
            w = log(1.0 - W)
            mu = -1.5861 - 0.31082 * ln - 0.083751 * ln ** 2 + 0.0038915 * ln ** 3
            sigma = exp(-0.4803 - 0.082676 * ln + 0.0030302 * ln ** 2)
        zscore = (w - mu) / sigma
        pvalue = 1.0 - normcdf(zscore)

    return {'W': W, 'pvalue': pvalue, 'zscore': zscore,
            'reject': bool(pvalue < alpha), 'nobs': int(n),
            'analyzer': 'shapirowilk'}
