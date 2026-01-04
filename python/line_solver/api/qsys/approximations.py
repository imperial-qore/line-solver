"""
G/G/1 and G/G/k approximation algorithms.

Native Python implementations for various approximations of general
queueing systems.
"""

import numpy as np
from typing import Dict, Tuple

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    def njit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator


def qsys_gig1_approx_allencunneen(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Dict[str, float]:
    """
    Allen-Cunneen approximation for G/G/1 queue.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        dict: Performance measures including L, Lq, W, Wq, rho
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

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Allen-Cunneen formula
    Lq = rho**2 * (ca2 + cs2) / (2 * (1 - rho))
    L = Lq + rho

    Wq = Lq / lambda_val
    W = L / lambda_val

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_gig1_approx_kingman(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Dict[str, float]:
    """
    Kingman's approximation for G/G/1 queue.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        dict: Performance measures including L, Lq, W, Wq, rho
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

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Kingman's heavy-traffic approximation
    Wq = (rho / (1 - rho)) * (1 / mu) * (ca2 + cs2) / 2
    Lq = lambda_val * Wq

    L = Lq + rho
    W = L / lambda_val

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_gig1_approx_marchal(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Dict[str, float]:
    """
    Marchal's approximation for G/G/1 queue.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        dict: Performance measures including L, Lq, W, Wq, rho
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

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Marchal's approximation with correction factor
    if cs2 >= 1:
        eta = 1.0
    else:
        eta = np.exp(-2 * (1 - rho) * (1 - cs2)**2 / (3 * rho * (ca2 + cs2)))

    Lq = eta * rho**2 * (ca2 + cs2) / (2 * (1 - rho))
    L = Lq + rho

    Wq = Lq / lambda_val
    W = L / lambda_val

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_gig1_approx_whitt(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Dict[str, float]:
    """
    Whitt's approximation for G/G/1 queue.

    Uses QNA (Queueing Network Analyzer) approximation.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        dict: Performance measures including L, Lq, W, Wq, rho
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

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Whitt's correction factor
    if ca2 <= 1 and cs2 <= 1:
        phi = np.exp(-2 * (1 - rho) * (1 - ca2)**2 / (3 * rho * (ca2 + cs2)))
    elif ca2 > 1 and cs2 <= 1:
        phi = np.exp(-(1 - rho) * (ca2 - 1) / (ca2 + 4 * cs2))
    elif ca2 <= 1 and cs2 > 1:
        phi = 1.0
    else:  # ca2 > 1 and cs2 > 1
        phi = 1.0

    Lq = phi * rho**2 * (ca2 + cs2) / (2 * (1 - rho))
    L = Lq + rho

    Wq = Lq / lambda_val
    W = L / lambda_val

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_gig1_approx_heyman(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Dict[str, float]:
    """
    Heyman's approximation for G/G/1 queue.

    A simple two-moment approximation.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        dict: Performance measures including W, Wq, rho, rhohat
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return {
            'W': np.inf,
            'Wq': np.inf,
            'rho': rho,
            'rhohat': 1.0
        }

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Heyman's formula
    W = rho / (1 - rho) / mu * (ca2 + cs2) / 2 + 1.0 / mu
    Wq = W - 1.0 / mu
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return {
        'W': W,
        'Wq': Wq,
        'rho': rho,
        'rhohat': rhohat
    }


def qsys_gig1_approx_kobayashi(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Dict[str, float]:
    """
    Kobayashi's approximation for G/G/1 queue.

    Uses an exponential decay formulation for the modified utilization.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        dict: Performance measures including W, Wq, rho, rhohat
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return {
            'W': np.inf,
            'Wq': np.inf,
            'rho': rho,
            'rhohat': 1.0
        }

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Kobayashi's formula
    rhohat = np.exp(-2 * (1 - rho) / (rho * (ca2 + cs2 / rho)))
    W = rhohat / (1 - rhohat) / lambda_val if rhohat < 1 else np.inf
    Wq = W - 1.0 / mu if np.isfinite(W) else np.inf

    return {
        'W': W,
        'Wq': Wq,
        'rho': rho,
        'rhohat': rhohat
    }


def qsys_gig1_approx_gelenbe(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Dict[str, float]:
    """
    Gelenbe's approximation for G/G/1 queue.

    A simple two-moment approximation.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        dict: Performance measures including W, Wq, rho
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return {
            'W': np.inf,
            'Wq': np.inf,
            'rho': rho
        }

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Gelenbe's formula
    W = (rho * ca2 + cs2) / 2.0 / (1.0 - rho) / lambda_val
    Wq = W - 1.0 / mu

    return {
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_gig1_approx_kimura(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Dict[str, float]:
    """
    Kimura's approximation for G/G/1 queue.

    Modified utilization-based approximation.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        dict: Performance measures including W, Wq, rho
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return {
            'W': np.inf,
            'Wq': np.inf,
            'rho': rho
        }

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Kimura's formula (using sigma = rho)
    W = rho * (ca2 + cs2) / mu / (1.0 - rho) / (1.0 + ca2) if ca2 >= 0 else np.inf
    Wq = W - 1.0 / mu if np.isfinite(W) else np.inf

    return {
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_gigk_approx(
    lambda_val: float, mu: float, k: int, ca: float, cs: float
) -> Dict[str, float]:
    """
    Approximation for G/G/k queue.

    Uses Allen-Cunneen approximation extended to multi-server.

    Args:
        lambda_val: Arrival rate
        mu: Service rate per server
        k: Number of servers
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        dict: Performance measures including L, Lq, W, Wq, rho
    """
    rho = lambda_val / (k * mu)

    if rho >= 1.0:
        return {
            'L': np.inf,
            'Lq': np.inf,
            'W': np.inf,
            'Wq': np.inf,
            'rho': rho
        }

    from .basic import qsys_mmk

    # Get M/M/k values
    mmk = qsys_mmk(lambda_val, mu, k)

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Allen-Cunneen correction
    Lq_mmk = mmk['Lq']
    Lq = Lq_mmk * (ca2 + cs2) / 2

    L = Lq + lambda_val / mu
    Wq = Lq / lambda_val
    W = L / lambda_val

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': rho
    }


def qsys_gig1_approx_klb(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Dict[str, float]:
    """
    Kraemer-Langenbach-Belz (KLB) approximation for G/G/1 queue.

    Reference: Kraemer, W., & Langenbach-Belz, M. (1976).
    "Approximate Formulas for the Delay in the Queuing System G/G/1".
    In: Proceedings of the 8th International Teletraffic Congress.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        dict: Performance measures including W (waiting time) and rhohat (effective utilization)
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return {
            'W': np.inf,
            'rhohat': 1.0,
            'rho': rho
        }

    ca2 = ca ** 2
    cs2 = cs ** 2

    # KLB formula with approximation function g
    if ca <= 1:
        g = np.exp(-2 * (1 - rho) * (1 - ca2) ** 2 / (3 * rho * (ca2 + cs2)))
    else:
        g = np.exp(-(1 - rho) * (ca2 - 1) / (ca2 + 4 * cs2))

    W = (1 / mu) * ((rho / (1 - rho)) * ((cs2 + ca2) / 2) * g + 1)
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return {
        'W': W,
        'rhohat': rhohat,
        'rho': rho
    }


def qsys_gig1_approx_myskja(
    lambda_val: float, mu: float, ca: float, cs: float,
    q0: float, qa: float
) -> float:
    """
    Myskja's approximation for G/G/1 queue.

    Uses third moments of arrival and service distributions for improved accuracy.

    Reference: Myskja, A. (1991).
    "An Experimental Study of a H₂/H₂/1 Queue".
    Stochastic Models, 7(4), 571-595.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time
        q0: Lowest relative third moment for given mean and SCV
        qa: Third relative moment E[X^3]/6/E[X]^3 of inter-arrival time

    Returns:
        float: Waiting time W
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Myskja formula incorporating third moments
    W = (rho / (2 * mu * (1 - rho))) * (
        (1 + cs2) + (q0 / qa) ** (1 / rho - rho) * (1 / rho) * (ca2 - 1)
    )

    return W


def qsys_gig1_approx_myskja2(
    lambda_val: float, mu: float, ca: float, cs: float,
    q0: float, qa: float
) -> float:
    """
    Modified Myskja (Myskja2) approximation for G/G/1 queue.

    An improved version of Myskja's approximation with better accuracy
    for a wider range of parameter combinations.

    Reference: Myskja, A. (1991).
    "An Experimental Study of a H₂/H₂/1 Queue".
    Stochastic Models, 7(4), 571-595.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time
        q0: Lowest relative third moment for given mean and SCV
        qa: Third relative moment E[X^3]/6/E[X]^3 of inter-arrival time

    Returns:
        float: Waiting time W
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Intermediate calculations
    ra = (1 + ca2) / 2
    rs = (1 + cs2) / 2
    theta = (rho * (qa - ra) - (qa - ra ** 2)) / (2 * rho * (ra - 1))
    d = (1 + 1 / ra) * (1 - rs) * (1 - (q0 / qa) ** 3) * (1 - rho ** 3)
    D = (rs - theta) ** 2 + (2 * rs - 1 + d) * (ra - 1)

    # Myskja2 formula
    W = (rho / (1 - rho)) / lambda_val * (rs + (1 / rho) * (np.sqrt(D) - (rs - theta)))

    return W


def qsys_gig1_ubnd_kingman(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Tuple[float, float]:
    """
    Kingman's upper bound for G/G/1 queue waiting time.

    This is an upper bound approximation that provides conservative
    estimates for the mean waiting time.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Upper bound on mean response time
            rhohat: Effective utilization (so M/M/1 formulas still hold)

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_gig1_ubnd_kingman.m
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf, 1.0

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Kingman's upper bound formula
    W = rho / (1 - rho) * (ca2 + cs2) / 2 * (1 / mu) + (1 / mu)
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat


def qsys_gigk_approx_kingman(
    lambda_val: float, mu: float, ca: float, cs: float, k: int
) -> Tuple[float, float]:
    """
    Kingman's approximation for G/G/k queue waiting time.

    Extends Kingman's approximation to multi-server queues using
    M/M/k waiting time as a base.

    Args:
        lambda_val: Arrival rate
        mu: Service rate per server
        k: Number of servers
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Approximate mean response time
            rhohat: Effective utilization

    References:
        Original MATLAB: matlab/src/api/qsys/qsys_gigk_approx_kingman.m
    """
    from .basic import qsys_mmk

    rho = lambda_val / (k * mu)

    if rho >= 1.0:
        return np.inf, 1.0

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Get M/M/k waiting time
    mmk_result = qsys_mmk(lambda_val, mu, k)
    W_mmk = mmk_result['W']

    # Kingman's approximation for G/G/k
    W = (ca2 + cs2) / 2 * (W_mmk - 1 / mu) + 1 / mu
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat
