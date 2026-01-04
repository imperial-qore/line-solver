"""
MAP Derivative Computations.

Native Python implementations for computing derivatives of MAP
probability distributions used in joint moments analysis.

Key algorithms:
    map_ccdf_derivative: Derivative at zero of MAP CCDF
    map_jointpdf_derivative: Partial derivative at zero of MAP joint PDF

References:
    Original MATLAB: matlab/src/api/mam/map_*.m
    A. Horvath et al., "A Joint Moments Based Analysis of Networks
    of MAP/MAP/1 Queues", QEST 2008
"""

import numpy as np
from typing import List, Union
from .map_analysis import map_pie


def map_ccdf_derivative(MAP: List[np.ndarray], i: int) -> float:
    """
    Compute derivative at zero of a MAP's complementary CDF.

    Calculates the i-th derivative at t=0 of the complementary
    cumulative distribution function (CCDF) of the inter-arrival
    time distribution.

    Formula: ν_i = π_e * D_0^i * e

    where π_e is the embedded stationary vector and e is the
    column vector of ones.

    Args:
        MAP: Markovian Arrival Process as [D0, D1]
        i: Order of the derivative

    Returns:
        Value of the i-th derivative at zero

    References:
        Original MATLAB: matlab/src/api/mam/map_ccdf_derivative.m
        A. Horvath et al., "A Joint Moments Based Analysis of
        Networks of MAP/MAP/1 Queues"
    """
    D0 = np.asarray(MAP[0], dtype=np.float64)
    n = D0.shape[0]

    pie = map_pie(MAP)

    # Compute π_e * D_0^i * e
    D0_power = np.linalg.matrix_power(D0, i)
    nu = pie @ D0_power @ np.ones(n)

    return float(nu)


def map_jointpdf_derivative(MAP: List[np.ndarray],
                            iset: Union[List[int], np.ndarray]) -> float:
    """
    Compute partial derivative at zero of a MAP's joint PDF.

    Calculates the partial derivative at t=0 of the joint probability
    density function of consecutive inter-arrival times.

    For index set {i_1, i_2, ..., i_k}:
    γ = π_e * D_0^{i_1} * D_1 * D_0^{i_2} * D_1 * ... * D_0^{i_k} * D_1 * e

    Args:
        MAP: Markovian Arrival Process as [D0, D1]
        iset: Vector of indices defining the partial derivative orders

    Returns:
        Value of the partial derivative at zero

    References:
        Original MATLAB: matlab/src/api/mam/map_jointpdf_derivative.m
        A. Horvath et al., "A Joint Moments Based Analysis of
        Networks of MAP/MAP/1 Queues"
    """
    D0 = np.asarray(MAP[0], dtype=np.float64)
    D1 = np.asarray(MAP[1], dtype=np.float64)
    n = D0.shape[0]

    iset = np.asarray(iset, dtype=int).flatten()

    gamma = map_pie(MAP)

    for j in iset:
        D0_power = np.linalg.matrix_power(D0, j)
        gamma = gamma @ D0_power @ D1

    gamma_val = gamma @ np.ones(n)

    return float(gamma_val)


def map_factorial_moment(MAP: List[np.ndarray], k: int) -> float:
    """
    Compute the k-th factorial moment of a MAP.

    The k-th factorial moment is computed using derivatives of
    the inter-arrival time distribution.

    Args:
        MAP: Markovian Arrival Process as [D0, D1]
        k: Order of the factorial moment

    Returns:
        k-th factorial moment

    References:
        Based on MAP moment formulas from matrix-analytic methods
    """
    D0 = np.asarray(MAP[0], dtype=np.float64)
    D1 = np.asarray(MAP[1], dtype=np.float64)
    n = D0.shape[0]

    pie = map_pie(MAP)

    # Compute (-D0)^{-k} iteratively
    try:
        neg_D0_inv = np.linalg.inv(-D0)
    except np.linalg.LinAlgError:
        neg_D0_inv = np.linalg.pinv(-D0)

    result = pie.copy()
    for _ in range(k):
        result = result @ neg_D0_inv

    moment_k = np.math.factorial(k) * result @ np.ones(n)

    return float(moment_k)


def map_joint_moment(MAP: List[np.ndarray], k: int, l: int) -> float:
    """
    Compute the (k,l)-th joint moment of consecutive inter-arrival times.

    E[X_n^k * X_{n+1}^l] for a MAP with inter-arrival times X_n.

    Args:
        MAP: Markovian Arrival Process as [D0, D1]
        k: Power for first inter-arrival time
        l: Power for second inter-arrival time

    Returns:
        Joint moment E[X_n^k * X_{n+1}^l]

    References:
        Based on joint moment formulas for MAPs
    """
    D0 = np.asarray(MAP[0], dtype=np.float64)
    D1 = np.asarray(MAP[1], dtype=np.float64)
    n = D0.shape[0]

    pie = map_pie(MAP)

    try:
        neg_D0_inv = np.linalg.inv(-D0)
    except np.linalg.LinAlgError:
        neg_D0_inv = np.linalg.pinv(-D0)

    # E[X_n^k * X_{n+1}^l] = k! * l! * π_e * (-D0)^{-k} * D1 * (-D0)^{-l} * e
    result = pie.copy()
    for _ in range(k):
        result = result @ neg_D0_inv

    result = result @ D1

    for _ in range(l):
        result = result @ neg_D0_inv

    joint_moment = (np.math.factorial(k) * np.math.factorial(l) *
                    result @ np.ones(n))

    return float(joint_moment)


__all__ = [
    'map_ccdf_derivative',
    'map_jointpdf_derivative',
    'map_factorial_moment',
    'map_joint_moment',
]
