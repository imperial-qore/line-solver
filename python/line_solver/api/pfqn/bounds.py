"""
Asymptotic Bounds for Product-Form Queueing Networks.

Implements various bounds on throughput and queue lengths for closed
queueing networks, including:
- Balanced Job Bounds (Zahorjan)
- Asymptotic Bounds (Zahorjan-Gittelsohn-Bryant)
- ZGSB Bounds (Zahorjan-Gittelsohn-Schweitzer-Bryant)
"""

import numpy as np
from typing import Union


def pfqn_xzabalow(
    L: np.ndarray,
    N: Union[int, float],
    Z: float
) -> float:
    """
    Lower asymptotic bound on throughput (Zahorjan-Balanced).

    Provides a simple lower bound on system throughput for single-class
    closed queueing networks.

    Args:
        L: Service demand vector (M,).
        N: Population (scalar).
        Z: Think time.

    Returns:
        Lower bound on throughput.
    """
    L = np.asarray(L, dtype=float).ravel()
    Ltot = np.sum(L)
    return float(N / (Z + Ltot * N))


def pfqn_xzabaup(
    L: np.ndarray,
    N: Union[int, float],
    Z: float
) -> float:
    """
    Upper asymptotic bound on throughput (Zahorjan-Balanced).

    Provides a simple upper bound on system throughput for single-class
    closed queueing networks based on bottleneck analysis.

    Args:
        L: Service demand vector (M,).
        N: Population (scalar).
        Z: Think time.

    Returns:
        Upper bound on throughput.
    """
    L = np.asarray(L, dtype=float).ravel()
    return float(min(1.0 / np.max(L), N / (np.sum(L) + Z)))


def pfqn_qzgblow(
    L: np.ndarray,
    N: Union[int, float],
    Z: float,
    i: int
) -> float:
    """
    Lower asymptotic bound on queue length (Zahorjan-Gittelsohn-Bryant).

    Args:
        L: Service demand vector (M,).
        N: Population (scalar).
        Z: Think time.
        i: Station index (0-based).

    Returns:
        Lower bound on mean queue length at station i.
    """
    L = np.asarray(L, dtype=float).ravel()
    yi = N * L[i] / (Z + np.sum(L) + np.max(L) * N)

    if yi >= 1:
        return float(N)

    Qgb = yi / (1 - yi) - (yi ** (N + 1)) / (1 - yi)
    return float(max(0, Qgb))


def pfqn_qzgbup(
    L: np.ndarray,
    N: Union[int, float],
    Z: float,
    i: int
) -> float:
    """
    Upper asymptotic bound on queue length (Zahorjan-Gittelsohn-Bryant).

    Args:
        L: Service demand vector (M,).
        N: Population (scalar).
        Z: Think time.
        i: Station index (0-based).

    Returns:
        Upper bound on mean queue length at station i.
    """
    L = np.asarray(L, dtype=float).ravel()
    sigma = np.sum(L ** 2) / np.sum(L)

    # Compute upper bound on throughput at N-1
    if N > 1:
        X_N_minus_1 = pfqn_xzabaup(L, N - 1, Z)
    else:
        X_N_minus_1 = 0

    Yi = L[i] * min(
        1.0 / np.max(L),
        N / (Z + np.sum(L) + sigma * (N - 1 - Z * X_N_minus_1))
    )

    if Yi < 1:
        Qgb = Yi / (1 - Yi) - (Yi ** (N + 1)) / (1 - Yi)
        return float(max(0, Qgb))
    else:
        return float(N)


def pfqn_xzgsblow(
    L: np.ndarray,
    N: Union[int, float],
    Z: float
) -> float:
    """
    Lower asymptotic bound on throughput (Zahorjan-Gittelsohn-Schweitzer-Bryant).

    Provides a tighter lower bound than pfqn_xzabalow by accounting for
    queue length bounds.

    Args:
        L: Service demand vector (M,).
        N: Population (scalar).
        Z: Think time.

    Returns:
        Lower bound on throughput.
    """
    L = np.asarray(L, dtype=float).ravel()
    M = len(L)
    max_L = np.max(L)

    R = Z + np.sum(L) + max_L * (N - 1)

    for i in range(M):
        if L[i] < max_L:
            R = R + (L[i] - max_L) * pfqn_qzgblow(L, N - 1, Z, i)

    discriminant = R ** 2 - 4 * Z * max_L * (N - 1)
    if discriminant < 0:
        # Fall back to simple bound
        return pfqn_xzabalow(L, N, Z)

    X = 2 * N / (R + np.sqrt(discriminant))
    return float(X)


def pfqn_xzgsbup(
    L: np.ndarray,
    N: Union[int, float],
    Z: float
) -> float:
    """
    Upper asymptotic bound on throughput (Zahorjan-Gittelsohn-Schweitzer-Bryant).

    Provides a tighter upper bound than pfqn_xzabaup by accounting for
    queue length bounds.

    Args:
        L: Service demand vector (M,).
        N: Population (scalar).
        Z: Think time.

    Returns:
        Upper bound on throughput.
    """
    L = np.asarray(L, dtype=float).ravel()
    M = len(L)
    max_L = np.max(L)

    R = Z + np.sum(L) + max_L * (N - 1)

    for i in range(M):
        if L[i] < max_L:
            R = R + (L[i] - max_L) * pfqn_qzgbup(L, N - 1, Z, i)

    discriminant = R ** 2 - 4 * Z * max_L * N
    if discriminant < 0:
        # Fall back to simple bound
        return pfqn_xzabaup(L, N, Z)

    X = 2 * N / (R + np.sqrt(discriminant))
    return float(X)


__all__ = [
    'pfqn_xzabalow',
    'pfqn_xzabaup',
    'pfqn_qzgblow',
    'pfqn_qzgbup',
    'pfqn_xzgsblow',
    'pfqn_xzgsbup',
]
