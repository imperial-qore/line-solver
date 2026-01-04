"""
Limited Joint-Dependence (LJD) Indexing Functions.

Native Python implementations of linearization/delinearization functions
for LJD scaling tables and multi-linear interpolation.

Key functions:
    ljd_linearize: Convert population vector to linear index
    ljd_delinearize: Convert linear index to population vector
    ljcd_interpolate: Multi-linear interpolation for scaling tables

References:
    Original MATLAB: matlab/src/api/pfqn/ljd_*.m, ljcd_interpolate.m
"""

import numpy as np
from typing import Union


def ljd_linearize(nvec: np.ndarray, cutoffs: np.ndarray) -> int:
    """
    Convert per-class population vector to linearized index.

    Maps a multi-dimensional population vector to a single linear index
    for efficient table lookups in LJD scaling tables.

    Index formula: idx = 1 + n1 + n2*(N1+1) + n3*(N1+1)*(N2+1) + ...

    Args:
        nvec: Per-class populations [n1, n2, ..., nK]
        cutoffs: Per-class cutoffs [N1, N2, ..., NK]

    Returns:
        1-based linearized index

    References:
        Original MATLAB: matlab/src/api/pfqn/ljd_linearize.m
    """
    nvec = np.asarray(nvec, dtype=int).flatten()
    cutoffs = np.asarray(cutoffs, dtype=int).flatten()

    K = len(nvec)
    idx = 1  # 1-indexed for MATLAB compatibility
    multiplier = 1

    for k in range(K):
        nk = min(nvec[k], cutoffs[k])  # Clamp to cutoff
        idx += nk * multiplier
        multiplier *= (cutoffs[k] + 1)

    return idx


def ljd_delinearize(idx: int, cutoffs: np.ndarray) -> np.ndarray:
    """
    Convert linearized index back to per-class population vector.

    Inverse of ljd_linearize: recovers the multi-dimensional population
    vector from a linear index.

    Args:
        idx: 1-based linearized index
        cutoffs: Per-class cutoffs [N1, N2, ..., NK]

    Returns:
        Per-class populations [n1, n2, ..., nK]

    References:
        Original MATLAB: matlab/src/api/pfqn/ljd_delinearize.m
    """
    cutoffs = np.asarray(cutoffs, dtype=int).flatten()

    K = len(cutoffs)
    nvec = np.zeros(K, dtype=int)
    idx = idx - 1  # Convert to 0-based for computation

    for k in range(K):
        base = cutoffs[k] + 1
        nvec[k] = idx % base
        idx = idx // base

    return nvec


def ljcd_interpolate(nvec: np.ndarray, cutoffs: np.ndarray,
                     table: np.ndarray, K: int) -> float:
    """
    Multi-linear interpolation for LJCD scaling tables.

    Performs multi-linear interpolation of a throughput value from an LJCD
    scaling table for non-integer population vectors.

    For K classes, interpolates between 2^K corner points of the hypercube
    containing the population vector using multi-linear interpolation.

    Args:
        nvec: Continuous population vector [n1, n2, ..., nK] (clamped to cutoffs)
        cutoffs: Per-class cutoffs [N1, N2, ..., NK]
        table: Linearized throughput table (1-D array indexed by ljd_linearize)
        K: Number of classes

    Returns:
        Interpolated throughput value

    References:
        Original MATLAB: matlab/src/api/pfqn/ljcd_interpolate.m
    """
    nvec = np.asarray(nvec, dtype=float).flatten()
    cutoffs = np.asarray(cutoffs, dtype=int).flatten()
    table = np.asarray(table, dtype=float).flatten()

    # Get floor and ceiling for each dimension
    nFloor = np.floor(nvec).astype(int)
    nCeil = np.ceil(nvec).astype(int)

    # Clamp to valid range
    nFloor = np.maximum(0, np.minimum(nFloor, cutoffs))
    nCeil = np.maximum(0, np.minimum(nCeil, cutoffs))

    # Compute fractional parts (weights)
    frac = nvec - nFloor

    # Handle edge case: all integer values
    if np.all(np.abs(frac) < 1e-10):
        idx = ljd_linearize(nFloor, cutoffs)
        if idx <= len(table):
            return table[idx - 1]  # Convert to 0-based
        else:
            return 0.0

    # Multi-linear interpolation over 2^K corners
    Xval = 0.0
    numCorners = 2 ** K

    for corner in range(numCorners):
        # Build corner point: bit i determines floor (0) or ceil (1) for class i
        cornerPoint = nFloor.copy()
        weight = 1.0

        for i in range(K):
            if (corner >> i) & 1:
                # Use ceiling for this dimension
                cornerPoint[i] = nCeil[i]
                weight *= frac[i]
            else:
                # Use floor for this dimension
                weight *= (1 - frac[i])

        # Look up table value at corner point
        idx = ljd_linearize(cornerPoint, cutoffs)
        if idx <= len(table):
            Xval += weight * table[idx - 1]  # Convert to 0-based

    return Xval


def infradius_h(x: np.ndarray, L: np.ndarray, N: np.ndarray,
                alpha: np.ndarray) -> np.ndarray:
    """
    Helper function for infinite radius computation with logistic transformation.

    Used in normalizing constant computation via integration methods.

    Args:
        x: Logistic transformation parameters
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        alpha: Load-dependent rate matrix

    Returns:
        Evaluated function value for integration

    References:
        Original MATLAB: matlab/src/api/pfqn/infradius_h.m
    """
    # Import here to avoid circular dependency
    from .ncld import pfqn_gld

    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=float).flatten()
    x = np.atleast_2d(np.asarray(x, dtype=float))

    M = L.shape[0]
    Nt = int(np.sum(N))
    beta = N / Nt if Nt > 0 else np.zeros_like(N)

    y = np.zeros(x.shape[0])

    for i in range(x.shape[0]):
        xi = x[i, :]

        # Logistic transformation
        t = np.exp(xi) / (1 + np.exp(xi))
        tb = np.sum(beta * t)

        # Evaluate h function
        z = np.sum(L * np.tile(np.exp(2 * np.pi * 1j * (t - tb)), (M, 1)), axis=1)
        gld_result = pfqn_gld(z, Nt, alpha)
        # Extract scalar G value from PfqnNcResult
        gld_value = gld_result.G if hasattr(gld_result, 'G') else gld_result

        # Jacobian of transformation
        jacobian = np.prod(np.exp(xi) / (1 + np.exp(xi)) ** 2)

        y[i] = np.real(gld_value * jacobian)

    return y


def infradius_hnorm(x: np.ndarray, L: np.ndarray, N: np.ndarray,
                    alpha: np.ndarray, logNormConstScale: float) -> np.ndarray:
    """
    Normalized helper function for infinite radius computation.

    Like infradius_h but with normalization for numerical stability.

    Args:
        x: Logistic transformation parameters
        L: Service demand matrix (M x R)
        N: Population vector (R,)
        alpha: Load-dependent rate matrix
        logNormConstScale: Logarithm of normalizing scale factor

    Returns:
        Normalized evaluated function value

    References:
        Original MATLAB: matlab/src/api/pfqn/infradius_hnorm.m
    """
    y = infradius_h(x, L, N, alpha)

    # Normalize
    if logNormConstScale != 0:
        y = y / np.exp(logNormConstScale)

    return y


__all__ = [
    'ljd_linearize',
    'ljd_delinearize',
    'ljcd_interpolate',
    'infradius_h',
    'infradius_hnorm',
]
