"""
PFQN Utility Functions.

Native Python implementations of utility functions for Product-Form
Queueing Network analysis, including load-dependent scaling, multiserver
scaling, and normalizing constant sanitization.

Key functions:
    pfqn_lldfun: Load and queue-dependent scaling function
    pfqn_mu_ms: Multi-server load-dependent service rate
    pfqn_nc_sanitize: Sanitize and preprocess NC parameters
    pfqn_cdfun: Class-dependence function
    multichoose: Combinations with repetition
    matchrow: Find row index in matrix
    oner: Increment vector at position

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_*.m
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.special import gammaln
from typing import Tuple, Optional, List
from itertools import combinations_with_replacement


def factln(n: float) -> float:
    """
    Compute log(n!) using log-gamma function.

    Args:
        n: Non-negative number

    Returns:
        log(n!) = log(Gamma(n+1))
    """
    if n <= 0:
        return 0.0
    return gammaln(n + 1)


def factln_vec(arr: np.ndarray) -> np.ndarray:
    """
    Compute log(n!) element-wise for an array.

    Args:
        arr: Array of non-negative numbers

    Returns:
        Array of log(n!) values
    """
    arr = np.asarray(arr, dtype=float)
    result = np.zeros_like(arr)
    mask = arr > 0
    result[mask] = gammaln(arr[mask] + 1)
    return result


def softmin(a: float, b: float, alpha: float = 20.0) -> float:
    """
    Compute a smooth approximation to min(a, b).

    Uses the log-sum-exp trick for numerical stability.

    Args:
        a: First value
        b: Second value
        alpha: Smoothing parameter (larger = sharper approximation)

    Returns:
        Smooth approximation of min(a, b)
    """
    if alpha <= 0:
        return min(a, b)

    # Use log-sum-exp for numerical stability
    max_val = max(-alpha * a, -alpha * b)
    return -1.0 / alpha * (max_val + np.log(np.exp(-alpha * a - max_val) +
                                              np.exp(-alpha * b - max_val)))


def oner(n: np.ndarray, s: int) -> np.ndarray:
    """
    Return a copy of n with position s reduced by 1.

    Args:
        n: Population vector
        s: Index to decrement (0-based)

    Returns:
        Copy of n with n[s] -= 1
    """
    n_copy = np.asarray(n, dtype=float).copy()
    if 0 <= s < len(n_copy):
        n_copy[s] = max(0, n_copy[s] - 1)
    return n_copy


def multichoose(r: int, n: int) -> np.ndarray:
    """
    Generate all combinations with repetition.

    Returns all ways to choose n items from r categories with repetition,
    where the result is a matrix with each row being a combination.

    Args:
        r: Number of categories
        n: Number of items to choose

    Returns:
        Matrix (C x r) where C = C(n+r-1, r-1) is the number of combinations
    """
    if r <= 0 or n < 0:
        return np.array([[]])

    if n == 0:
        return np.zeros((1, r), dtype=int)

    # Generate all combinations with repetition
    # Each combination represents how many of each category
    result = []

    def generate_combinations(remaining: int, categories_left: int, current: List[int]):
        if categories_left == 1:
            result.append(current + [remaining])
            return
        for i in range(remaining + 1):
            generate_combinations(remaining - i, categories_left - 1, current + [i])

    generate_combinations(n, r, [])

    return np.array(result, dtype=int)


def matchrow(matrix: np.ndarray, row: np.ndarray) -> int:
    """
    Find the index of a row in a matrix.

    Args:
        matrix: 2D array to search in
        row: 1D array to find

    Returns:
        1-based index of the matching row, or 0 if not found
    """
    matrix = np.atleast_2d(np.asarray(matrix))
    row = np.asarray(row).flatten()

    for i in range(matrix.shape[0]):
        if np.array_equal(matrix[i, :len(row)], row):
            return i + 1  # 1-based indexing like MATLAB
    return 0


def pfqn_lldfun(n: np.ndarray, lldscaling: Optional[np.ndarray] = None,
                nservers: Optional[np.ndarray] = None) -> np.ndarray:
    """
    AMVA-QD load and queue-dependent scaling function.

    Computes the scaling factor for load-dependent queueing stations,
    accounting for multi-server stations and general load-dependent
    service rate scaling.

    Args:
        n: Queue population vector (M,)
        lldscaling: Load-dependent scaling matrix (M x Nmax), optional
        nservers: Number of servers per station (M,), optional

    Returns:
        Scaling factor vector (M,)

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_lldfun.m
    """
    n = np.asarray(n, dtype=float).flatten()
    M = len(n)
    r = np.ones(M)
    alpha = 20.0  # softmin parameter

    if lldscaling is not None:
        lldscaling = np.atleast_2d(np.asarray(lldscaling, dtype=float))
        smax = lldscaling.shape[1]
    else:
        smax = 0

    for i in range(M):
        # Handle servers
        if nservers is not None:
            nservers_arr = np.asarray(nservers, dtype=float).flatten()
            if i < len(nservers_arr):
                if np.isinf(nservers_arr[i]):
                    # Delay server
                    r[i] = 1.0
                else:
                    # Multi-server station
                    sm = softmin(n[i], nservers_arr[i], alpha)
                    if np.isnan(sm) or sm <= 0:
                        # Fallback if numerical problems in soft-min
                        sm = min(n[i], nservers_arr[i])
                        if sm <= 0:
                            sm = 1.0
                    r[i] = r[i] / sm

        # Handle generic load-dependent scaling
        if lldscaling is not None and smax > 0:
            if i < lldscaling.shape[0]:
                row = lldscaling[i, :smax]
                # Check if there's variation in the scaling
                if np.ptp(row) > 0:  # ptp = max - min = range
                    # Interpolate the scaling factor
                    x_vals = np.arange(1, smax + 1)
                    try:
                        interp_func = interp1d(x_vals, row, kind='cubic',
                                              fill_value='extrapolate')
                        scale_val = float(interp_func(n[i]))
                        if scale_val > 0 and np.isfinite(scale_val):
                            r[i] = r[i] / scale_val
                    except Exception:
                        # Fallback to linear interpolation
                        scale_val = np.interp(n[i], x_vals, row)
                        if scale_val > 0 and np.isfinite(scale_val):
                            r[i] = r[i] / scale_val

    return r


def pfqn_mu_ms(N: int, m: int, c: int) -> np.ndarray:
    """
    Compute load-dependent rates for m identical c-server FCFS stations.

    Calculates the effective service rate as a function of the number
    of jobs in the system for a network of m identical stations, each
    with c parallel servers.

    Args:
        N: Maximum population
        m: Number of identical stations
        c: Number of servers per station

    Returns:
        Load-dependent service rate vector (1 x N)

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_mu_ms.m
    """
    if N <= 0:
        return np.array([])

    mu = np.zeros(N)
    g = np.zeros((m, N + 1))  # g[station, population]

    # Auxiliary function to compute g values
    def gnaux(n: int, station: int, c: int, g: np.ndarray) -> float:
        if n == 0:
            return 1.0
        elif station == 0:  # First station (0-indexed)
            # Product of min(1:n, c)
            vals = np.minimum(np.arange(1, n + 1), c)
            return 1.0 / np.prod(vals)
        else:
            result = 0.0
            for k in range(n + 1):
                if k == 0:
                    a = 1.0
                else:
                    vals = np.minimum(np.arange(1, k + 1), c)
                    a = 1.0 / np.prod(vals)

                b = g[station - 1, n - k]
                result += a * b
            return result

    # Fill in the g table
    for n in range(N + 1):
        for station in range(m):
            g[station, n] = gnaux(n, station, c, g)

    # Compute mu from g
    for n in range(1, N + 1):
        if g[m - 1, n] > 0:
            mu[n - 1] = g[m - 1, n - 1] / g[m - 1, n]
        else:
            mu[n - 1] = 0.0

    return mu


def pfqn_nc_sanitize(lam: np.ndarray, L: np.ndarray, N: np.ndarray,
                     Z: np.ndarray, atol: float = 1e-8
                     ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                np.ndarray, float]:
    """
    Sanitize and preprocess network parameters for NC solvers.

    Removes empty/ill-defined classes, rescales demands for numerical
    stability, and reorders classes by think time.

    Args:
        lam: Arrival rate vector
        L: Service demand matrix (M x R)
        N: Population vector
        Z: Think time vector
        atol: Absolute tolerance for numerical comparisons

    Returns:
        Tuple of (lambda, L, N, Z, lGremaind) where:
            - lambda: Sanitized arrival rates
            - L: Sanitized service demands (rescaled)
            - N: Sanitized populations
            - Z: Sanitized think times (rescaled)
            - lGremaind: Log normalization factor from removed classes

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_nc_sanitize.m
    """
    lam = np.asarray(lam, dtype=float).flatten()
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    Z = np.asarray(Z, dtype=float).flatten()

    # Replace NaN with 0
    L = np.nan_to_num(L, nan=0.0)
    Z = np.nan_to_num(Z, nan=0.0)

    # Erase empty classes (N == 0)
    nnzclasses = np.where(N > 0)[0]
    if len(nnzclasses) == 0:
        return lam, L, N, Z, 0.0

    L = L[:, nnzclasses]
    N = N[nnzclasses]
    Z = Z[nnzclasses] if len(Z) > 0 else np.zeros(len(nnzclasses))
    lam = lam[nnzclasses] if len(lam) > 0 else np.zeros(len(nnzclasses))

    # Erase ill-defined classes (sum of L and Z < atol)
    class_totals = np.sum(L, axis=0) + Z
    valid_classes = np.where(class_totals >= atol)[0]
    if len(valid_classes) < L.shape[1]:
        L = L[:, valid_classes]
        N = N[valid_classes]
        Z = Z[valid_classes]
        lam = lam[valid_classes]

    lGremaind = 0.0

    # Find zero demand classes (all stations have demand < atol)
    R = L.shape[1]
    zerodemands = []
    for r in range(R):
        if np.max(L[:, r]) < atol:
            zerodemands.append(r)

    if len(zerodemands) > 0:
        zerodemands = np.array(zerodemands)
        # Contribution from zero-demand classes
        for r in zerodemands:
            if Z[r] > 0:
                lGremaind += N[r] * np.log(Z[r]) - factln(N[r])

        # Remove zero-demand classes
        keep = np.setdiff1d(np.arange(R), zerodemands)
        if len(keep) > 0:
            L = L[:, keep]
            Z = Z[keep]
            N = N[keep]
            lam = lam[keep] if len(lam) > 0 else np.array([])
        else:
            return lam, np.zeros((L.shape[0], 0)), np.array([]), np.array([]), lGremaind

    if L.shape[1] == 0:
        return lam, L, N, Z, lGremaind

    # Rescale demands for numerical stability
    Lmax = np.max(L, axis=0)
    Lmax[Lmax < atol] = 1.0  # Avoid division by zero

    L = L / Lmax
    Z = Z / Lmax
    lGremaind += np.dot(N, np.log(Lmax))

    # Sort from smallest to largest think time
    if len(Z) > 0 and np.sum(Z) > 0:
        rsort = np.argsort(Z)
        L = L[:, rsort]
        Z = Z[rsort]
        N = N[rsort]
        if len(lam) > 0:
            lam = lam[rsort]

    # Ensure zero think time classes are first
    zerothinktimes = np.where(Z < atol)[0]
    nonzerothinktimes = np.setdiff1d(np.arange(L.shape[1]), zerothinktimes)
    order = np.concatenate([zerothinktimes, nonzerothinktimes])

    L = L[:, order]
    N = N[order]
    Z = Z[order]
    if len(lam) > 0:
        lam = lam[order]

    return lam, L, N, Z, lGremaind


def pfqn_cdfun(nvec: np.ndarray, cdscaling: Optional[List] = None) -> np.ndarray:
    """
    AMVA-QD class-dependence function for queue-dependent scaling.

    Computes the scaling factor at each station based on the class-dependent
    population state.

    Args:
        nvec: Population state matrix (M x R) or vector (M,)
        cdscaling: List of class-dependent scaling functions, one per station.
                  Each function takes a class population vector and returns a scalar.

    Returns:
        Scaling factor vector (M,)

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_cdfun.m
    """
    nvec = np.atleast_2d(np.asarray(nvec, dtype=float))
    M = nvec.shape[0]
    r = np.ones(M)

    if cdscaling is not None and len(cdscaling) > 0:
        for i in range(M):
            if i < len(cdscaling) and cdscaling[i] is not None:
                try:
                    scale_val = cdscaling[i](nvec[i, :])
                    if scale_val > 0 and np.isfinite(scale_val):
                        r[i] = 1.0 / scale_val
                except Exception:
                    pass

    return r


def pfqn_ljdfun(n: np.ndarray, ljdscaling: Optional[np.ndarray] = None) -> np.ndarray:
    """
    AMVA-QD limited joint-dependence function.

    Computes scaling factors for stations with limited joint-dependence,
    where the service rate depends on the total population across stations.

    Args:
        n: Queue population vector (M,)
        ljdscaling: Limited joint-dependence scaling matrix (M x Nmax), optional

    Returns:
        Scaling factor vector (M,)

    References:
        Original MATLAB: matlab/src/api/pfqn/pfqn_ljdfun.m
    """
    n = np.asarray(n, dtype=float).flatten()
    M = len(n)
    r = np.ones(M)

    if ljdscaling is None:
        return r

    ljdscaling = np.atleast_2d(np.asarray(ljdscaling, dtype=float))
    smax = ljdscaling.shape[1]

    # Total population
    nt = int(np.sum(n))

    for i in range(M):
        if i < ljdscaling.shape[0] and nt > 0:
            # Get the scaling value at total population
            idx = min(nt, smax) - 1
            if idx >= 0 and ljdscaling[i, idx] > 0:
                r[i] = 1.0 / ljdscaling[i, idx]

    return r


__all__ = [
    'factln',
    'factln_vec',
    'softmin',
    'oner',
    'multichoose',
    'matchrow',
    'pfqn_lldfun',
    'pfqn_mu_ms',
    'pfqn_nc_sanitize',
    'pfqn_cdfun',
    'pfqn_ljdfun',
]
