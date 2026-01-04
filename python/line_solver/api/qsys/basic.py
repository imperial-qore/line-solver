"""
Basic single-queue system analysis algorithms.

Native Python implementations for M/M/1, M/M/k, M/G/1, G/M/1, and M/G/inf queues.
These implementations match the API of the JPype wrapper for compatibility.
"""

import numpy as np
from math import factorial
from typing import Dict, Optional, Any

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    def njit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator


@njit(cache=True)
def _erlang_c(k: int, rho: float) -> float:
    """
    Erlang-C formula (probability all servers busy).

    Args:
        k: Number of servers
        rho: Utilization per server (lambda / (k * mu))

    Returns:
        Probability that an arriving customer must wait
    """
    if rho >= 1.0:
        return 1.0

    # Sum from j=0 to k-1 of (k*rho)^j / j!
    kr = k * rho
    S = 0.0
    term = 1.0  # (k*rho)^0 / 0! = 1
    S += term
    for j in range(1, k):
        term *= kr / j
        S += term

    # (k*rho)^k / k!
    numerator_term = term * kr / k

    # C = numerator_term / (numerator_term + (1-rho) * S)
    C = numerator_term / (numerator_term + (1 - rho) * S)
    return C


@njit(cache=True)
def _erlang_b(k: int, a: float) -> float:
    """
    Erlang-B formula (blocking probability for M/M/k/k).

    Args:
        k: Number of servers (and capacity)
        a: Offered load (lambda / mu)

    Returns:
        Blocking probability
    """
    # Recursive formula: B(k,a) = a*B(k-1,a) / (k + a*B(k-1,a))
    B = 1.0
    for i in range(1, k + 1):
        B = a * B / (i + a * B)
    return B


def qsys_mm1(lambda_val: float, mu: float) -> Dict[str, float]:
    """
    Analyze M/M/1 queue (Poisson arrivals, exponential service).

    Args:
        lambda_val: Arrival rate (lambda)
        mu: Service rate

    Returns:
        dict: Performance measures including:
            - L: Mean number in system
            - Lq: Mean number in queue
            - W: Mean response time (time in system)
            - Wq: Mean waiting time (time in queue)
            - rho: Utilization (lambda/mu)

    Example:
        >>> result = qsys_mm1(0.5, 1.0)
        >>> print(f"Utilization: {result['rho']:.2f}")
        Utilization: 0.50
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return {
            'L': np.inf,
            'Lq': np.inf,
            'W': np.inf,
            'Wq': np.inf,
            'rho': rho
        }

    # Little's law and M/M/1 formulas
    L = rho / (1 - rho)
    Lq = rho**2 / (1 - rho)
    W = 1 / (mu - lambda_val)  # = L / lambda
    Wq = rho / (mu - lambda_val)  # = Lq / lambda

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_mmk(lambda_val: float, mu: float, k: int) -> Dict[str, float]:
    """
    Analyze M/M/k queue (Poisson arrivals, k exponential servers).

    Args:
        lambda_val: Arrival rate (lambda)
        mu: Service rate per server
        k: Number of parallel servers

    Returns:
        dict: Performance measures including:
            - L: Mean number in system
            - Lq: Mean number in queue
            - W: Mean response time
            - Wq: Mean waiting time
            - rho: Utilization per server (lambda/(k*mu))
            - P0: Probability of empty system

    Example:
        >>> result = qsys_mmk(2.0, 1.0, 3)
        >>> print(f"Utilization: {result['rho']:.2f}")
        Utilization: 0.67
    """
    rho = lambda_val / (mu * k)
    a = lambda_val / mu  # Offered load

    if rho >= 1.0:
        return {
            'L': np.inf,
            'Lq': np.inf,
            'W': np.inf,
            'Wq': np.inf,
            'rho': rho,
            'P0': 0.0
        }

    # Erlang-C formula
    C = _erlang_c(k, rho)

    # Queue length in queue
    Lq = C * rho / (1 - rho)

    # Total in system
    L = Lq + a

    # Waiting times via Little's law
    Wq = Lq / lambda_val
    W = L / lambda_val

    # P0: probability of empty system
    # Sum of (k*rho)^j/j! for j=0 to k-1, plus (k*rho)^k/(k!(1-rho))
    kr = k * rho
    sum_term = 0.0
    term = 1.0
    sum_term += term
    for j in range(1, k):
        term *= kr / j
        sum_term += term
    term *= kr / k  # Now term = (k*rho)^k / k!
    sum_term += term / (1 - rho)
    P0 = 1.0 / sum_term

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho,
        'P0': P0
    }


def qsys_mg1(lambda_val: float, mu: float, cs: float) -> Dict[str, float]:
    """
    Analyze M/G/1 queue using Pollaczek-Khinchine formula.

    Args:
        lambda_val: Arrival rate (lambda)
        mu: Service rate (mean service time = 1/mu)
        cs: Coefficient of variation of service time (std/mean)

    Returns:
        dict: Performance measures including:
            - L: Mean number in system
            - Lq: Mean number in queue
            - W: Mean response time
            - Wq: Mean waiting time
            - rho: Utilization (lambda/mu)

    Example:
        >>> result = qsys_mg1(0.5, 1.0, 1.0)  # cs=1 is exponential (M/M/1)
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return {
            'L': np.inf,
            'Lq': np.inf,
            'W': np.inf,
            'Wq': np.inf,
            'rho': rho
        }

    # Pollaczek-Khinchine formula for Lq
    # Lq = (rho^2 + lambda^2 * Var[S]) / (2 * (1 - rho))
    # where Var[S] = cs^2 / mu^2
    cs2 = cs ** 2
    var_s = cs2 / (mu ** 2)

    Lq = (rho**2 + lambda_val**2 * var_s) / (2 * (1 - rho))
    L = Lq + rho

    # Waiting times via Little's law
    Wq = Lq / lambda_val
    W = L / lambda_val

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_gm1(lambda_val: float, mu: float, ca: float) -> Dict[str, float]:
    """
    Analyze G/M/1 queue (general arrivals, exponential service).

    Args:
        lambda_val: Arrival rate (lambda)
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time

    Returns:
        dict: Performance measures including:
            - L: Mean number in system
            - Lq: Mean number in queue
            - W: Mean response time
            - Wq: Mean waiting time
            - rho: Utilization

    Note:
        Uses approximation based on the coefficient of variation.
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return {
            'L': np.inf,
            'Lq': np.inf,
            'W': np.inf,
            'Wq': np.inf,
            'rho': rho
        }

    # For G/M/1, we need to solve for sigma where sigma = A*(mu*(1-sigma))
    # where A*(s) is the LST of inter-arrival times
    # For general G, we use approximation

    ca2 = ca ** 2

    # Kramer-Langenbach-Belz approximation
    if ca2 <= 1:
        g = np.exp(-2 * (1 - rho) * (1 - ca2) / (3 * rho * (ca2 + 1)))
    else:
        g = np.exp(-(1 - rho) * (ca2 - 1) / (ca2 + 4 * rho))

    # M/M/1 values
    L_mm1 = rho / (1 - rho)
    Lq_mm1 = rho**2 / (1 - rho)

    # Apply correction factor
    Lq = g * Lq_mm1
    L = Lq + rho

    # Waiting times via Little's law
    Wq = Lq / lambda_val
    W = L / lambda_val

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_mminf(lambda_val: float, mu: float) -> Dict[str, float]:
    """
    Analyze M/M/inf queue (infinite servers / delay station).

    Args:
        lambda_val: Arrival rate (lambda)
        mu: Service rate

    Returns:
        dict: Performance measures including:
            - L: Mean number in system (= lambda/mu)
            - Lq: Mean number in queue (= 0)
            - W: Mean time in system (= 1/mu)
            - Wq: Mean waiting time (= 0)
            - P0: Probability of empty system
    """
    rho = lambda_val / mu

    return {
        'L': rho,
        'Lq': 0.0,
        'W': 1.0 / mu,
        'Wq': 0.0,
        'P0': np.exp(-rho)
    }


def qsys_mginf(lambda_val: float, mu: float, k: Optional[int] = None) -> Dict[str, Any]:
    """
    Analyze M/G/inf queue (infinite servers, general service).

    Performance is independent of service time distribution shape.
    Number of customers follows Poisson distribution.

    Args:
        lambda_val: Arrival rate (lambda)
        mu: Service rate (mean service time = 1/mu)
        k: Optional state for probability computation

    Returns:
        dict: Performance measures including:
            - L: Mean number in system
            - Lq: Mean number in queue (= 0)
            - W: Mean time in system (= 1/mu)
            - Wq: Mean waiting time (= 0)
            - P0: Probability of empty system
            - Pk: Probability of k customers (if k provided)
    """
    rho = lambda_val / mu

    result = {
        'L': rho,
        'Lq': 0.0,
        'W': 1.0 / mu,
        'Wq': 0.0,
        'P0': np.exp(-rho)
    }

    if k is not None:
        # Poisson probability P(X=k) = exp(-rho) * rho^k / k!
        result['Pk'] = np.exp(-rho) * (rho ** k) / factorial(k)

    return result
