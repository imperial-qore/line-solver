"""
Laplace Approximation for Normalizing Constant Computation.

Implements the Laplace approximation method for computing normalizing
constants in closed queueing networks, including the NRP (Normal Radius-Probit)
and NRL (Normal Radius-Logistic) methods.
"""

import numpy as np
from typing import Tuple, Callable, Optional
from math import log, exp, pi


def num_hess(f: Callable, x0: np.ndarray, tol: float = 1e-5) -> np.ndarray:
    """
    Compute numerical Hessian matrix using finite differences.

    Uses central differences for improved accuracy.

    Args:
        f: Function to differentiate (maps x -> scalar)
        x0: Point at which to compute Hessian
        tol: Step size for finite differences

    Returns:
        Hessian matrix (d x d) where d = len(x0)
    """
    x0 = np.atleast_1d(np.asarray(x0, dtype=float))
    d = len(x0)
    H = np.zeros((d, d))

    h = tol

    for i in range(d):
        for j in range(d):
            if i == j:
                # Diagonal: second derivative
                x_plus = x0.copy()
                x_minus = x0.copy()
                x_plus[i] += h
                x_minus[i] -= h

                f_plus = f(x_plus)
                f_minus = f(x_minus)
                f_0 = f(x0)

                H[i, i] = (f_plus - 2 * f_0 + f_minus) / (h * h)
            else:
                # Off-diagonal: mixed partial derivative
                x_pp = x0.copy()
                x_pm = x0.copy()
                x_mp = x0.copy()
                x_mm = x0.copy()

                x_pp[i] += h
                x_pp[j] += h
                x_pm[i] += h
                x_pm[j] -= h
                x_mp[i] -= h
                x_mp[j] += h
                x_mm[i] -= h
                x_mm[j] -= h

                H[i, j] = (f(x_pp) - f(x_pm) - f(x_mp) + f(x_mm)) / (4 * h * h)

    # Ensure symmetry
    H = (H + H.T) / 2

    return H


def laplaceapprox(h: Callable, x0: np.ndarray,
                  tol: float = 1e-5) -> Tuple[float, np.ndarray, float]:
    """
    Laplace approximation for multidimensional integrals.

    Approximates I = ∫ h(x) dx using the Laplace method:
        I ≈ h(x0) * sqrt((2π)^d / det(-H))

    where H is the Hessian of log(h) at x0.

    Args:
        h: Function to integrate (must be positive)
        x0: Point for Laplace approximation (typically the mode)
        tol: Tolerance for numerical Hessian computation

    Returns:
        Tuple (I, H, logI) where:
            I: Approximate integral value
            H: Hessian matrix at x0
            logI: Logarithm of integral value (more stable)
    """
    x0 = np.atleast_1d(np.asarray(x0, dtype=float))
    d = len(x0)

    # Compute Hessian of log(h) at x0
    def log_h(x):
        val = h(x)
        if isinstance(val, np.ndarray):
            val = val[0] if len(val) > 0 else 0.0
        return log(max(val, 1e-300))

    H = num_hess(log_h, x0, tol)

    # Compute determinant of negative Hessian
    detNegH = np.linalg.det(-H)

    # Handle numerical issues
    if detNegH < 0:
        # Try with larger tolerance
        H = num_hess(log_h, x0, tol * 10)
        detNegH = np.linalg.det(-H)

    if detNegH < 0:
        # Try even larger tolerance
        H = num_hess(log_h, x0, tol * 100)
        detNegH = np.linalg.det(-H)

    if detNegH < 0:
        # Last resort: use absolute value
        detNegH = abs(detNegH)

    # Evaluate h at x0
    h_x0 = h(x0)
    if isinstance(h_x0, np.ndarray):
        h_x0 = h_x0[0] if len(h_x0) > 0 else 0.0

    # Compute approximation
    if detNegH > 0 and h_x0 > 0:
        I = h_x0 * np.sqrt((2 * pi) ** d / detNegH)
        logI = log(h_x0) + (d / 2) * log(2 * pi) - 0.5 * log(detNegH)
    else:
        I = 0.0
        logI = -np.inf

    return I, H, logI


def pfqn_nrl(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
             alpha: np.ndarray = None) -> float:
    """
    Normal Radius-Logistic (NRL) approximation for normalizing constant.

    Computes the logarithm of the normalizing constant using Laplace
    approximation with logistic transformation.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,) - optional
        alpha: Load-dependent rate matrix (M x Ntot) - optional

    Returns:
        lG: Logarithm of normalizing constant

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_nrl.m
    """
    from .ljd import infradius_h
    from .ncld import pfqn_gld

    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()

    # Handle edge cases
    if np.sum(N) < 0:
        return -np.inf
    if np.sum(N) == 0:
        return 0.0

    M, R = L.shape
    Nt = int(np.sum(N))

    # Build alpha if not provided (default: identity service rates)
    if alpha is None:
        alpha = np.tile(np.arange(1, Nt + 1), (M, 1))

    # Add think time as an extra station
    if Z is not None:
        Z = np.asarray(Z, dtype=float).flatten()
        if np.sum(Z) > 0:
            L = np.vstack([L, Z])
            alpha_row = np.arange(1, Nt + 1)
            alpha = np.vstack([alpha, alpha_row])

    # Handle single-station case
    if M == 1 and (Z is None or np.sum(Z) == 0):
        _, lG = pfqn_gld(L, N, alpha)
        return lG

    # Scale demands
    Lmax = np.max(L, axis=0)
    Lmax[Lmax == 0] = 1.0
    L_scaled = L / Lmax

    # Initial point for Laplace approximation
    x0 = np.zeros(R)

    # Laplace approximation
    def h_func(x):
        return infradius_h(x.reshape(1, -1), L_scaled, N, alpha)

    _, _, lG = laplaceapprox(h_func, x0)

    # Adjust for scaling
    lG = np.real(lG + np.dot(N, np.log(Lmax)))

    return lG


def pfqn_nrp(L: np.ndarray, N: np.ndarray, Z: np.ndarray = None,
             alpha: np.ndarray = None) -> float:
    """
    Normal Radius-Probit (NRP) approximation for normalizing constant.

    Computes the logarithm of the normalizing constant using Laplace
    approximation with probit (normalized) transformation.

    Args:
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        Z: Think time vector (R,) - optional
        alpha: Load-dependent rate matrix (M x Ntot) - optional

    Returns:
        lG: Logarithm of normalizing constant

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_nrp.m
    """
    from .ljd import infradius_hnorm
    from .ncld import pfqn_gld

    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()

    # Handle edge cases
    if np.sum(N) < 0:
        return -np.inf
    if np.sum(N) == 0:
        return 0.0

    M, R = L.shape
    Nt = int(np.sum(N))

    # Build alpha if not provided
    if alpha is None:
        alpha = np.tile(np.arange(1, Nt + 1), (M, 1))

    # Add think time as an extra station
    if Z is not None:
        Z = np.asarray(Z, dtype=float).flatten()
        if np.sum(Z) > 0:
            L = np.vstack([L, Z])
            alpha_row = np.arange(1, Nt + 1)
            alpha = np.vstack([alpha, alpha_row])

    # Handle single-station case
    if M == 1 and (Z is None or np.sum(Z) == 0):
        _, lG = pfqn_gld(L, N, alpha)
        return lG

    # Scale demands
    Lmax = np.max(L, axis=0)
    Lmax[Lmax == 0] = 1.0
    L_scaled = L / Lmax

    # Initial point for Laplace approximation
    x0 = np.zeros(R)

    # Get a normalization scale from initial evaluation
    from .ljd import infradius_h
    h0 = infradius_h(x0.reshape(1, -1), L_scaled, N, alpha)
    logNormConstScale = np.log(max(h0[0], 1e-300)) if len(h0) > 0 else 0.0

    # Laplace approximation with normalized function
    def h_func(x):
        return infradius_hnorm(x.reshape(1, -1), L_scaled, N, alpha, logNormConstScale)

    _, _, lG_normalized = laplaceapprox(h_func, x0)

    # Adjust for normalization and scaling
    lG = np.real(lG_normalized + logNormConstScale + np.dot(N, np.log(Lmax)))

    return lG


def pfqn_lap(L: np.ndarray, N: np.ndarray, Z: np.ndarray) -> float:
    """
    Laplace approximation for normalizing constant.

    Computes the logarithm of the normalizing constant using the
    classical Laplace approximation with root finding.

    This method uses a saddle-point approximation to estimate the
    normalizing constant integral representation.

    Args:
        L: Service demand vector (R,) - single station
        N: Population vector (R,)
        Z: Think time vector (R,)

    Returns:
        lG: Logarithm of normalizing constant approximation

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_lap.m
    """
    from scipy.optimize import brentq
    from scipy.special import gammaln

    L = np.asarray(L, dtype=float).flatten()
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()

    R = len(N)

    # Compute total population
    Ntot = np.sum(N)

    if Ntot <= 0:
        return 0.0

    # Log factorial function
    def factln(n):
        if n <= 0:
            return 0.0
        return float(gammaln(n + 1))

    # Define the equation to solve for the saddle point
    # f(x) = 1 - sum(N*L / (Z + Ntot*L*x))
    def f(x):
        denom = Z + Ntot * L * x
        # Avoid division by zero
        denom = np.maximum(denom, 1e-300)
        return 1 - np.sum(N * L / denom)

    # Find root u0
    try:
        # Try brentq first for robust root finding
        # Look for sign change in a reasonable interval
        x_low = 1e-10
        x_high = 10.0

        f_low = f(x_low)
        f_high = f(x_high)

        if f_low * f_high < 0:
            u0 = brentq(f, x_low, x_high)
        else:
            # Try to find a bracket
            found = False
            for x_test in [0.1, 0.5, 1.0, 2.0, 5.0]:
                f_test = f(x_test)
                if f_low * f_test < 0:
                    u0 = brentq(f, x_low, x_test)
                    found = True
                    break
                elif f_test * f_high < 0:
                    u0 = brentq(f, x_test, x_high)
                    found = True
                    break

            if not found:
                # Fall back to grid search
                u0 = 1.0
                init_sign = np.sign(f(0.001)) if f(0.001) != 0 else 1.0
                for x_test in np.arange(1e-4, 10, 1e-4):
                    f_val = f(x_test)
                    if f_val != 0 and np.sign(f_val) != init_sign:
                        u0 = x_test
                        break
    except Exception:
        # Grid search fallback
        u0 = 1.0
        try:
            init_sign = np.sign(f(0.001)) if f(0.001) != 0 else 1.0
            for x_test in np.arange(1e-4, 10, 1e-4):
                f_val = f(x_test)
                if f_val != 0 and np.sign(f_val) != init_sign:
                    u0 = x_test
                    break
        except Exception:
            pass

    # Check if u0 is valid
    if u0 < 0 or not np.isfinite(u0):
        return np.nan

    # Compute log of normalizing constant using Laplace approximation
    # logI = log(Ntot) - sum(factln(N)) - Ntot*u0 + sum(N*log(Z+L*u0*Ntot))
    #      + 0.5*log(2*pi) - 0.5*log(f''(u0)) - 0.5*log(Ntot)

    # Coefficient term
    logI = np.log(Ntot) - sum(factln(n) for n in N)

    # f(u0) term
    terms = Z + L * u0 * Ntot
    terms = np.maximum(terms, 1e-300)
    logI = logI - Ntot * u0 + np.sum(N * np.log(terms))

    # Second derivative term: f''(u0) = sum((N/Ntot) / (Z/(Ntot*L) + u0)^2)
    denom = Z / (Ntot * L + 1e-300) + u0
    denom = np.maximum(denom, 1e-300)
    f2 = np.sum((N / Ntot) / (denom ** 2))

    if f2 <= 0:
        f2 = 1e-300

    logI = logI + 0.5 * np.log(2 * pi) - 0.5 * np.log(f2)

    # Ntot correction term
    logI = logI - 0.5 * np.log(Ntot)

    return float(logI)


__all__ = [
    'num_hess',
    'laplaceapprox',
    'pfqn_nrl',
    'pfqn_nrp',
    'pfqn_lap',
]
