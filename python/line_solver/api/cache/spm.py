"""
Singular Perturbation Method (SPM) for Cache Analysis.

Native Python implementations of SPM techniques for cache system analysis
using partial differential equation methods. Provides numerical solutions
for cache performance in high-dimensional parameter spaces.

References:
    Gast, N. and Van Houdt, B., "Transient and Steady-state Regime of a
    Family of List-based Cache Replacement Algorithms." SIGMETRICS 2015.
"""

import numpy as np
from scipy import special
from typing import Tuple


def cache_xi_iter(gamma: np.ndarray, m: np.ndarray,
                  tol: float = 1e-12, max_iter: int = 1000) -> np.ndarray:
    """
    Compute cache xi terms using iterative method from Gast-van Houdt.

    This method calculates the xi values which are important for understanding
    the distribution of items in the cache. The algorithm assumes that the
    access factors are monotone with the list index.

    Args:
        gamma: Cache access factors matrix (n x h).
        m: Cache capacity vector (h,).
        tol: Convergence tolerance (default: 1e-12).
        max_iter: Maximum iterations (default: 1000).

    Returns:
        Vector of xi values (h,).
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    n = gamma.shape[0]  # Number of items
    h = len(m)  # Number of cache levels

    # Normalize capacity
    f = m / n

    # Build pp matrix (h+1 x n)
    # pp[0, :] = ones, pp[j+1, :] = gamma[:, j]
    pp = np.zeros((h + 1, n))
    pp[0, :] = 1.0
    for i in range(h):
        pp[i + 1, :] = gamma[:, i]

    z_old = np.zeros(h + 1)
    z = np.ones(h + 1)

    iteration = 0
    while np.max(np.abs(z - z_old)) > tol * np.max(np.abs(z_old)) and iteration < max_iter:
        z_old = z.copy()

        temp = n * (z @ pp)  # Element-wise: n * sum_j z[j] * pp[j, :]

        for i in range(h):
            # a = temp - n * z[i+1] * pp[i+1, :]
            a = temp - n * z[i + 1] * pp[i + 1, :]

            # Fi = sum_k (pp[i+1, k] / (n * pp[i+1, k] + a[k]))
            denom = n * pp[i + 1, :] + a
            denom = np.where(np.abs(denom) < 1e-14, 1e-14, denom)
            Fi = np.sum(pp[i + 1, :] / denom)

            # Bisection to find z[i+1]
            if Fi > f[i]:
                zi_min = 0.0
                zi_max = 1.0
            else:
                zi_min = 1.0
                zi_max = 2.0
                # Expand search range with safety limit to prevent infinite loop
                expand_iter = 0
                max_expand = 100  # Prevent infinite loop if condition never met
                while expand_iter < max_expand:
                    denom_test = n * zi_max * pp[i + 1, :] + a
                    denom_test = np.where(np.abs(denom_test) < 1e-14, 1e-14, denom_test)
                    Fi_test = np.sum(zi_max * pp[i + 1, :] / denom_test)
                    if Fi_test >= f[i] or not np.isfinite(Fi_test) or not np.isfinite(f[i]):
                        break
                    zi_min = zi_max
                    zi_max *= 2
                    expand_iter += 1

            # Bisection search
            for _ in range(50):
                zi = (zi_min + zi_max) / 2
                z[i + 1] = zi
                denom_zi = n * zi * pp[i + 1, :] + a
                denom_zi = np.where(np.abs(denom_zi) < 1e-14, 1e-14, denom_zi)
                Fi_zi = np.sum(zi * pp[i + 1, :] / denom_zi)
                if Fi_zi < f[i]:
                    zi_min = zi
                else:
                    zi_max = zi

        iteration += 1

    # Return xi values (without the first element)
    return z[1:]


def cache_spm(gamma: np.ndarray, m: np.ndarray
              ) -> Tuple[float, float, np.ndarray]:
    """
    Approximate the normalizing constant using SPM method.

    Computes the normalizing constant of the cache steady state distribution
    using the singular perturbation method (also known as ray integration).

    Args:
        gamma: Cache access factors matrix (n x h).
        m: Cache capacity vector (h,).

    Returns:
        Tuple of (Z, lZ, xi) where:
            - Z: Approximated normalizing constant
            - lZ: Log of normalizing constant
            - xi: Xi terms vector
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    # Filter out items with zero access
    row_sums = np.sum(gamma, axis=1)
    keep_rows = row_sums > 0
    gamma = gamma[keep_rows, :]

    h = len(m)  # Number of cache levels
    n = gamma.shape[0]  # Number of items
    mt = np.sum(m)  # Total capacity

    if n == mt:
        print("Warning: Number of items equals the cache capacity.")

    # Compute xi terms
    xi = cache_xi_iter(gamma, m)

    # Compute S_k = sum_l gamma[k, l] * xi[l]
    S = gamma @ xi  # Shape: (n,)

    # Compute phi = sum_k log(1 + S_k) - sum_l m_l * log(xi_l)
    phi = np.sum(np.log(1 + S)) - np.sum(m * np.log(xi))

    # Compute matrix C for determinant
    delta = np.eye(h)
    C = np.zeros((h, h))

    for j in range(h):
        for l in range(h):
            # C1 = sum_k gamma[k, j] / (1 + S_k)
            C1 = np.sum(gamma[:, j] / (1 + S))

            # C2 = sum_k gamma[k, j] * gamma[k, l] / (1 + S_k)^2
            C2 = np.sum(gamma[:, j] * gamma[:, l] / (1 + S)**2)

            C[j, l] = delta[j, l] * C1 - xi[j] * C2

    # Compute Z using the formula (matching MATLAB: lZ=real(lZ) to handle roundoff)
    det_C = np.linalg.det(C)

    # m! (factorial of each element)
    m_fact = special.factorial(m)
    m_fact_prod = np.prod(m_fact)

    # Log of Z
    lZ = (-h * np.log(np.sqrt(2 * np.pi)) + phi +
          np.sum(special.gammaln(m + 1)) -  # log of m!
          np.sum(0.5 * np.log(xi)) -
          0.5 * np.log(det_C + 0j))  # +0j to handle negative det via complex log
    lZ = np.real(lZ)  # Remove small imaginary parts from roundoff (MATLAB: lZ=real(lZ))

    # Z
    Z = np.exp(lZ)

    return Z, lZ, xi


def cache_prob_spm(gamma: np.ndarray, m: np.ndarray) -> np.ndarray:
    """
    Compute cache miss probabilities using singular perturbation method.

    Uses the SPM (ray integration) method to compute miss probabilities
    for cache systems.

    Args:
        gamma: Cache access factors matrix (n x h).
        m: Cache capacity vector (h,).

    Returns:
        Vector of miss probabilities for each item.
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    n = gamma.shape[0]
    h = gamma.shape[1]

    # Get xi values
    _, _, xi = cache_spm(gamma, m)

    # Compute S_k
    S = gamma @ xi

    # Miss probability for item k: 1 / (1 + S_k)
    miss_prob = 1.0 / (1 + S)

    return miss_prob


def cache_prob_fpi(gamma: np.ndarray, m: np.ndarray,
                   tol: float = 1e-8, max_iter: int = 1000) -> np.ndarray:
    """
    Compute cache hit probabilities using fixed point iteration (FPI method).

    Matches MATLAB cache_prob_fpi: Uses fixed point iteration to compute
    cache hit probability distribution.

    Args:
        gamma: Cache access factors matrix (n x h).
        m: Cache capacity vector (h,).
        tol: Convergence tolerance (unused, kept for API compatibility).
        max_iter: Maximum iterations (unused, kept for API compatibility).

    Returns:
        Matrix (n x h+1) where:
            - Column 0: miss probabilities
            - Columns 1:h: hit probabilities at each level
    """
    from .miss import cache_xi_fp

    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    n = gamma.shape[0]
    h = gamma.shape[1]

    # Get xi using fixed point method (matches MATLAB cache_xi_fp)
    xi, _, _, _ = cache_xi_fp(gamma, m)

    # Build probability matrix (n x h+1)
    # Column 0: miss probability, columns 1:h: hit at each level
    prob = np.zeros((n, h + 1))

    for i in range(n):
        # S_i = gamma(i,:) * xi
        S_i = np.dot(gamma[i, :], xi)

        # Miss probability: 1 / (1 + S_i)
        prob[i, 0] = 1.0 / (1.0 + S_i)

        # Hit probabilities at each level: gamma(i,j)*xi(j) / (1 + S_i)
        for j in range(h):
            prob[i, j + 1] = gamma[i, j] * xi[j] / (1.0 + S_i)

    return prob


__all__ = [
    'cache_xi_iter',
    'cache_spm',
    'cache_prob_spm',
    'cache_prob_fpi',
]
