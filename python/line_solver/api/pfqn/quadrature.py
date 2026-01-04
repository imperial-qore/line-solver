"""
Quadrature Methods for Normalizing Constant Computation.

Native Python implementations of numerical integration methods for
computing normalizing constants in product-form queueing networks.

Key functions:
    pfqn_mmint2: McKenna-Mitra integral using scipy.integrate
    pfqn_mmint2_gausslegendre: Gauss-Legendre quadrature
    pfqn_mmint2_gausslaguerre: Gauss-Laguerre quadrature
    pfqn_mmsample2: Monte Carlo sampling approximation

References:
    McKenna, J., and Mitra, D. "Integral representations and asymptotic
    expansions for closed Markovian queueing networks: Normal usage."
    Bell System Technical Journal 61.5 (1982): 661-683.
    Original MATLAB: matlab/src/api/pfqn/pfqn_mmint2*.m
"""

import numpy as np
from scipy import integrate
from scipy.special import gammaln
from typing import Tuple, Optional
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.laguerre import laggauss


def factln(n: float) -> float:
    """Compute log(n!) using log-gamma function."""
    if n <= 0:
        return 0.0
    return gammaln(n + 1)


def factln_vec(arr: np.ndarray) -> np.ndarray:
    """Compute log(n!) element-wise for an array."""
    arr = np.asarray(arr, dtype=float)
    result = np.zeros_like(arr)
    mask = arr > 0
    result[mask] = gammaln(arr[mask] + 1)
    return result


def logsumexp(x: np.ndarray) -> float:
    """
    Compute log(sum(exp(x))) in a numerically stable way.

    Args:
        x: Array of values

    Returns:
        log(sum(exp(x)))
    """
    x = np.asarray(x).flatten()
    x_max = np.max(x)
    if not np.isfinite(x_max):
        return float('-inf')
    return x_max + np.log(np.sum(np.exp(x - x_max)))


def pfqn_mmint2(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                m: int = 1) -> Tuple[float, float]:
    """
    McKenna-Mitra integral form using scipy.integrate.

    Computes the normalizing constant using numerical integration of the
    McKenna-Mitra integral representation.

    Args:
        L: Service demand vector (R,)
        N: Population vector (R,)
        Z: Think time vector (R,)
        m: Replication factor (default: 1)

    Returns:
        Tuple of (G, lG):
            - G: Normalizing constant
            - lG: Logarithm of normalizing constant

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_mmint2.m
    """
    L = np.asarray(L, dtype=float).flatten()
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()

    # Find non-zero classes
    nnzClasses = np.where(N > 0)[0]

    if len(nnzClasses) == 0:
        return 1.0, 0.0

    order = 12

    def f(u):
        """Integrand function."""
        terms = Z[nnzClasses] + L[nnzClasses] * u
        terms = np.maximum(terms, 1e-300)  # Avoid log(0)
        return np.exp(-u) * np.prod(np.power(terms, N[nnzClasses]))

    # Cutoff for exponential term (99.999...% percentile)
    p = 1 - 10**(-order)
    exp1prctile = -np.log(1 - p)

    # Perform integration
    try:
        result, _ = integrate.quad(f, 0, exp1prctile, limit=500)
        lG = np.log(max(result, 1e-300)) - np.sum(factln_vec(N))
    except Exception:
        # Fallback to simple quadrature
        result = 0.0
        n_points = 1000
        x = np.linspace(0, exp1prctile, n_points)
        for xi in x:
            result += f(xi)
        result *= exp1prctile / n_points
        lG = np.log(max(result, 1e-300)) - np.sum(factln_vec(N))

    G = np.exp(lG) if lG > -700 else 0.0

    return G, lG


def pfqn_mmint2_gausslegendre(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                               m: int = 1) -> Tuple[float, float]:
    """
    McKenna-Mitra integral with Gauss-Legendre quadrature.

    Uses Gauss-Legendre quadrature for improved accuracy in computing
    the McKenna-Mitra integral.

    Args:
        L: Service demand vector (R,)
        N: Population vector (R,)
        Z: Think time vector (R,)
        m: Replication factor (default: 1)

    Returns:
        Tuple of (G, lG):
            - G: Normalizing constant
            - lG: Logarithm of normalizing constant

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_mmint2_gausslegendre.m
    """
    L = np.asarray(L, dtype=float).flatten()
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()

    # Number of quadrature points (use at least 300)
    n = max(300, min(500, 2 * int(np.sum(N) + m - 1) - 1))

    # Get Gauss-Legendre nodes and weights on [-1, 1]
    nodes, weights = leggauss(n)

    # Transform to [0, upper_limit]
    # For repairman integrals, we use a finite upper limit
    # The integrand decays exponentially, so we use a large but finite bound
    upper_limit = 20 + np.sum(N)  # Adaptive upper limit

    # Transform nodes and weights: x = (u + 1) * upper_limit / 2
    x = (nodes + 1) * upper_limit / 2
    w = weights * upper_limit / 2

    # Compute integrand at each node
    y = np.zeros(n)
    for i in range(n):
        terms = Z + L * x[i]
        terms = np.maximum(terms, 1e-300)
        y[i] = np.dot(N, np.log(terms))

    # Compute log of integrand
    g = np.log(w) - x + y + (m - 1) * np.log(np.maximum(x, 1e-300))

    # Compute coefficient
    coeff = -np.sum(factln_vec(N)) - factln(m - 1)

    # Use logsumexp for numerical stability
    lG = logsumexp(g) + coeff

    G = np.exp(lG) if np.isfinite(lG) and lG > -700 else 0.0

    return G, lG


def pfqn_mmint2_gausslaguerre(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                               m: int = 1) -> Tuple[float, float]:
    """
    McKenna-Mitra integral with Gauss-Laguerre quadrature.

    Uses Gauss-Laguerre quadrature which is naturally suited for
    integrals with exponential decay.

    Args:
        L: Service demand vector (R,)
        N: Population vector (R,)
        Z: Think time vector (R,)
        m: Replication factor (default: 1)

    Returns:
        Tuple of (G, lG):
            - G: Normalizing constant
            - lG: Logarithm of normalizing constant

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_mmint2_gausslaguerre.m
    """
    L = np.asarray(L, dtype=float).flatten()
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()

    # Find non-zero classes
    nonzeroClasses = np.where(N > 0)[0]

    if len(nonzeroClasses) == 0:
        return 1.0, 0.0

    # Number of quadrature points
    n = min(300, 2 * int(np.sum(N)) + 1)
    n = max(n, 50)  # At least 50 points

    # Get Gauss-Laguerre nodes and weights
    # For Laguerre polynomials, the integral is exp(-x) * f(x) from 0 to inf
    x, w = laggauss(n)

    # Compute integrand values
    F = np.zeros(n)
    for i in range(n):
        terms = Z[nonzeroClasses] + L[nonzeroClasses] * x[i]
        terms = np.maximum(terms, 1e-300)
        F[i] = (m - 1) * np.log(max(x[i], 1e-300)) + np.dot(N[nonzeroClasses], np.log(terms))

    # Log of weighted integrand
    # Note: Gauss-Laguerre weights include the exp(-x) factor
    g = np.log(np.maximum(w, 1e-300)) + F - np.sum(factln_vec(N)) - factln(m - 1)

    # Use logsumexp for numerical stability
    lG = logsumexp(g)

    G = np.exp(lG) if np.isfinite(lG) and lG > -700 else 0.0

    return G, lG


def pfqn_mmsample2(L: np.ndarray, N: np.ndarray, Z: np.ndarray,
                   m: int = 1, n_samples: int = 10000) -> Tuple[float, float]:
    """
    Monte Carlo sampling approximation for normalizing constant.

    Uses Monte Carlo sampling to approximate the McKenna-Mitra integral.
    Useful for very large populations where quadrature becomes expensive.

    Args:
        L: Service demand vector (R,)
        N: Population vector (R,)
        Z: Think time vector (R,)
        m: Replication factor (default: 1)
        n_samples: Number of Monte Carlo samples (default: 10000)

    Returns:
        Tuple of (G, lG):
            - G: Normalizing constant (approximate)
            - lG: Logarithm of normalizing constant

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_mmsample2.m
    """
    L = np.asarray(L, dtype=float).flatten()
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()

    # Find non-zero classes
    nonzeroClasses = np.where(N > 0)[0]

    if len(nonzeroClasses) == 0:
        return 1.0, 0.0

    # Sample from exponential distribution (for importance sampling)
    # The integral is of the form: int exp(-u) * f(u) du from 0 to inf
    samples = np.random.exponential(1.0, n_samples)

    # Evaluate log of integrand at each sample
    log_values = np.zeros(n_samples)
    for i, u in enumerate(samples):
        terms = Z[nonzeroClasses] + L[nonzeroClasses] * u
        terms = np.maximum(terms, 1e-300)
        log_values[i] = (m - 1) * np.log(max(u, 1e-300)) + np.dot(N[nonzeroClasses], np.log(terms))

    # For importance sampling with exp(-u), the weights are uniform
    # (since we're already sampling from exp(-u))

    # Compute coefficient
    coeff = -np.sum(factln_vec(N)) - factln(m - 1)

    # Use logsumexp and subtract log(n_samples) for averaging
    lG = logsumexp(log_values) - np.log(n_samples) + coeff

    G = np.exp(lG) if np.isfinite(lG) and lG > -700 else 0.0

    return G, lG


__all__ = [
    'pfqn_mmint2',
    'pfqn_mmint2_gausslegendre',
    'pfqn_mmint2_gausslaguerre',
    'pfqn_mmsample2',
    'logsumexp',
]
