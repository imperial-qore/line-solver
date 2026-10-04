"""
Laplace-Stieltjes Transform (LST) Functions for AoI Analysis.

Native Python implementations of LST functions for common distributions
used in Age of Information analysis.

Key functions:
    aoi_lst_exp: LST for exponential distribution
    aoi_lst_erlang: LST for Erlang distribution
    aoi_lst_det: LST for deterministic (constant) distribution
    aoi_lst_ph: LST for phase-type distribution

References:
    Original MATLAB: matlab/src/api/aoi/aoi_lst_*.m
"""

import numpy as np
from typing import Callable


def aoi_lst_exp(mu: float) -> Callable[[np.ndarray], np.ndarray]:
    """
    Laplace-Stieltjes transform for exponential distribution.

    Returns a function for the LST of an exponential distribution
    with rate mu.

    Args:
        mu: Rate parameter (mean = 1/mu)

    Returns:
        LST function: H*(s) = mu / (mu + s)

    Raises:
        ValueError: If mu <= 0

    References:
        Original MATLAB: matlab/src/api/aoi/aoi_lst_exp.m
    """
    if mu <= 0:
        raise ValueError("Rate mu must be positive")

    def lst(s):
        s = np.asarray(s, dtype=complex)
        return mu / (mu + s)

    return lst


def aoi_lst_erlang(k: int, mu: float) -> Callable[[np.ndarray], np.ndarray]:
    """
    Laplace-Stieltjes transform for Erlang distribution.

    Returns a function for the LST of an Erlang-k distribution
    with rate parameter mu.

    Args:
        k: Shape parameter (number of phases, must be positive integer)
        mu: Rate parameter per phase (mean = k/mu)

    Returns:
        LST function: H*(s) = (mu / (mu + s))^k

    Raises:
        ValueError: If k < 1 or not an integer, or mu <= 0

    References:
        Original MATLAB: matlab/src/api/aoi/aoi_lst_erlang.m
    """
    if k < 1 or k != int(k):
        raise ValueError("Shape k must be a positive integer")
    if mu <= 0:
        raise ValueError("Rate mu must be positive")

    k = int(k)

    def lst(s):
        s = np.asarray(s, dtype=complex)
        return (mu / (mu + s)) ** k

    return lst


def aoi_lst_det(d: float) -> Callable[[np.ndarray], np.ndarray]:
    """
    Laplace-Stieltjes transform for deterministic (constant) distribution.

    Returns a function for the LST of a deterministic distribution
    with constant value d.

    Args:
        d: Constant value (must be positive)

    Returns:
        LST function: H*(s) = exp(-s * d)

    Raises:
        ValueError: If d <= 0

    References:
        Original MATLAB: matlab/src/api/aoi/aoi_lst_det.m
    """
    if d <= 0:
        raise ValueError("Constant d must be positive")

    def lst(s):
        s = np.asarray(s, dtype=complex)
        return np.exp(-s * d)

    return lst


def aoi_lst_ph(alpha: np.ndarray, T: np.ndarray) -> Callable[[np.ndarray], np.ndarray]:
    """
    Laplace-Stieltjes transform for phase-type distribution.

    Returns a function for the LST of a phase-type distribution
    with initial probability vector alpha and sub-generator T.

    Args:
        alpha: Initial probability vector (m,)
        T: Sub-generator matrix (m x m)

    Returns:
        LST function: H*(s) = alpha @ inv(s*I - T) @ (-T @ ones)

    Raises:
        ValueError: If dimensions mismatch

    References:
        Original MATLAB: matlab/src/api/aoi/aoi_lst_ph.m
    """
    alpha = np.asarray(alpha, dtype=float).flatten()
    T = np.asarray(T, dtype=float)

    m = len(alpha)
    if T.shape != (m, m):
        raise ValueError(f"Dimension mismatch: T shape {T.shape}, expected ({m}, {m})")

    # Exit rate vector
    t0 = -T @ np.ones(m)

    def lst(s):
        s_scalar = np.isscalar(s)
        s = np.atleast_1d(np.asarray(s, dtype=complex))

        result = np.zeros(len(s), dtype=complex)
        I = np.eye(m)

        for i, si in enumerate(s):
            try:
                inv_mat = np.linalg.solve(si * I - T, t0)
                result[i] = np.dot(alpha, inv_mat)
            except np.linalg.LinAlgError:
                result[i] = 0.0

        if s_scalar:
            return result[0]
        return result

    return lst


__all__ = [
    'aoi_lst_exp',
    'aoi_lst_erlang',
    'aoi_lst_det',
    'aoi_lst_ph',
]
