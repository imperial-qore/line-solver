"""
Normal and Student t quantiles used by the output-analysis routines.

Thin wrappers over SciPy, kept in one place so the MATLAB twins (which reach the
same values through ``erfc``, ``erfinv`` and ``betaincinv`` to avoid a toolbox
dependency) and the JAR twins (which use commons-math3) have a single point of
comparison. Agreement across the three is to within 1e-10.

References:
    Original MATLAB: matlab/src/api/sim/sim_normcdf.m, sim_norminv.m, sim_tinv.m
"""

from math import sqrt

from scipy.special import erfc, erfinv
from scipy.stats import t as _student_t

__all__ = ['normcdf', 'norminv', 'tinv']


def normcdf(z: float) -> float:
    """
    Standard normal cumulative distribution function.

    Args:
        z: The argument

    Returns:
        Phi(z)
    """
    return float(0.5 * erfc(-z / sqrt(2.0)))


def norminv(p: float) -> float:
    """
    Standard normal quantile function.

    Args:
        p: Probability in [0,1]

    Returns:
        The p-quantile, infinite at the endpoints.

    Raises:
        ValueError: If p lies outside [0,1].
    """
    if not 0.0 <= p <= 1.0:
        raise ValueError("Probability must lie in [0,1], got %r" % (p,))
    if p == 0.0:
        return float('-inf')
    if p == 1.0:
        return float('inf')
    return float(sqrt(2.0) * erfinv(2.0 * p - 1.0))


def tinv(p: float, nu: float) -> float:
    """
    Quantile function of Student's t distribution.

    Args:
        p: Probability in [0,1]
        nu: Degrees of freedom, positive

    Returns:
        The p-quantile of t with nu degrees of freedom.

    Raises:
        ValueError: If p lies outside [0,1] or nu is not positive.
    """
    if not nu > 0.0:
        raise ValueError("Degrees of freedom must be positive, got %r" % (nu,))
    if not 0.0 <= p <= 1.0:
        raise ValueError("Probability must lie in [0,1], got %r" % (p,))
    return float(_student_t.ppf(p, nu))
