"""
G/G/1 and G/G/k approximation algorithms.

Native Python implementations for various approximations of general
queueing systems.
"""

import numpy as np
from typing import Tuple

def qsys_gig1_approx_allencunneen(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Tuple[float, float]:
    """
    Allen-Cunneen approximation for G/G/1 queue.

    Matches MATLAB qsys_gig1_approx_allencunneen.m exactly.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Mean response time (time in system)
            rhohat: Effective utilization (so that M/M/1 formulas still hold)
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf, 1.0

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Allen-Cunneen formula (matches MATLAB qsys_gig1_approx_allencunneen.m)
    W = (rho / (1 - rho)) / mu * ((cs2 + ca2) / 2) + 1 / mu
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat

def qsys_gig1_approx_kingman(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Tuple[float, float]:
    """
    Kingman's upper bound approximation for G/G/1 queue.

    Note: In MATLAB, 'gig1' and 'gig1.kingman' map to qsys_gig1_ubnd_kingman.
    This function provides the same formula.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Mean response time (time in system)
            rhohat: Effective utilization (so that M/M/1 formulas still hold)
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf, 1.0

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Kingman's upper bound (matches MATLAB qsys_gig1_ubnd_kingman.m)
    W = rho / (1 - rho) * (ca2 + cs2) / 2 * (1 / mu) + (1 / mu)
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat

def qsys_gig1_approx_marchal(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Tuple[float, float]:
    """
    Marchal's approximation for G/G/1 queue.

    Matches MATLAB qsys_gig1_approx_marchal.m exactly.
    Note: MATLAB formula uses ca (not ca^2) in the numerator factor.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Mean response time (time in system)
            rhohat: Effective utilization (so that M/M/1 formulas still hold)
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf, 1.0

    cs2 = cs ** 2

    # Marchal's approximation (matches MATLAB qsys_gig1_approx_marchal.m exactly)
    # MATLAB: W = Wmm1*(1+cs^2)/2/mu*(ca+rho^2*cs^2)/(1+rho^2*cs^2)+1/mu
    Wmm1 = rho / (1 - rho)
    W = Wmm1 * (1 + cs2) / 2 / mu * (ca + rho**2 * cs2) / (1 + rho**2 * cs2) + 1 / mu
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat

def qsys_gig1_approx_whitt(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Tuple[float, float]:
    """
    Whitt's approximation for G/G/1 queue.

    Uses QNA (Queueing Network Analyzer) approximation.
    Note: No direct MATLAB counterpart (qsys_gig1_approx_whitt.m does not exist).

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Mean response time (time in system)
            rhohat: Effective utilization (so that M/M/1 formulas still hold)
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf, 1.0

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
    W = L / lambda_val

    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat

def qsys_gig1_approx_heyman(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Tuple[float, float]:
    """
    Heyman's approximation for G/G/1 queue.

    Matches MATLAB qsys_gig1_approx_heyman.m exactly.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Mean response time (time in system)
            rhohat: Effective utilization (so that M/M/1 formulas still hold)
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf, 1.0

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Heyman's formula (matches MATLAB)
    W = rho / (1 - rho) / mu * (ca2 + cs2) / 2 + 1.0 / mu
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat

def qsys_gig1_approx_kobayashi(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Tuple[float, float]:
    """
    Kobayashi's approximation for G/G/1 queue.

    Matches MATLAB qsys_gig1_approx_kobayashi.m exactly.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Mean response time (time in system)
            rhohat: Effective utilization
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf, 1.0

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Kobayashi's formula (matches MATLAB)
    rhohat = np.exp(-2 * (1 - rho) / (rho * (ca2 + cs2 / rho)))
    W = rhohat / (1 - rhohat) / lambda_val if rhohat < 1 else np.inf

    return W, rhohat

def qsys_gig1_approx_gelenbe(
    lambda_val: float, mu: float, ca: float, cs: float
) -> float:
    """
    Gelenbe's approximation for G/G/1 queue.

    Matches MATLAB qsys_gig1_approx_gelenbe.m exactly.
    Note: MATLAB returns only W (no rhohat).

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        float: Mean response time W
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Gelenbe's formula (matches MATLAB)
    W = (rho * ca2 + cs2) / 2.0 / (1.0 - rho) / lambda_val

    return W

def qsys_gig1_approx_kimura(
    lambda_val: float, mu: float, ca: float, cs: float
) -> float:
    """
    Kimura's approximation for G/G/1 queue.

    Matches MATLAB qsys_gig1_approx_kimura.m exactly.
    Note: MATLAB uses sigma (=rho) as first parameter, returns only W.

    Args:
        lambda_val: Arrival rate (used as sigma = lambda/mu = rho)
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        float: Mean response time W
    """
    # MATLAB: sigma = first arg, here we compute rho = lambda/mu
    sigma = lambda_val / mu

    if sigma >= 1.0:
        return np.inf

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Kimura's formula (matches MATLAB: W=sigma*(ca^2+cs^2)/mu/(1-sigma)/(1+ca^2))
    W = sigma * (ca2 + cs2) / mu / (1.0 - sigma) / (1.0 + ca2)

    return W

def qsys_gigk_approx(
    lambda_val: float, mu: float, ca: float, cs: float, k: int
) -> Tuple[float, float]:
    """
    Approximation for G/G/k queue.

    Matches MATLAB qsys_gigk_approx.m formula using alpha-factor correction.

    Args:
        lambda_val: Arrival rate
        mu: Service rate per server
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time
        k: Number of servers

    Returns:
        Tuple of (W, rhohat) where:
            W: Approximate mean response time
            rhohat: Effective utilization
    """
    rho = lambda_val / (k * mu)

    if rho >= 1.0:
        return np.inf, 1.0

    ca2 = ca ** 2
    cs2 = cs ** 2

    # MATLAB formula: alpha depends on rho
    if rho > 0.7:
        alpha = (rho**k + rho) / 2
    else:
        alpha = rho**((k + 1) / 2)

    W = (alpha / mu) * (1 / (1 - rho)) * (ca2 + cs2) / (2 * k) + 1 / mu
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat

def qsys_gig1_approx_klb(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Tuple[float, float]:
    """
    Kraemer-Langenbach-Belz (KLB) approximation for G/G/1 queue.

    Matches MATLAB qsys_gig1_approx_klb.m exactly.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Mean response time (time in system)
            rhohat: Effective utilization (so that M/M/1 formulas still hold)
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf, 1.0

    ca2 = ca ** 2
    cs2 = cs ** 2

    # KLB formula (matches MATLAB)
    if ca <= 1:
        g = np.exp(-2 * (1 - rho) * (1 - ca2) ** 2 / (3 * rho * (ca2 + cs2)))
    else:
        g = np.exp(-(1 - rho) * (ca2 - 1) / (ca2 + 4 * cs2))

    W = 1 / mu * ((rho / (1 - rho)) * ((cs2 + ca2) / 2) * g + 1)
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat

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

    # Handle edge case: ra = 1 (ca = 1, exponential arrivals)
    # theta = (rho*(qa-ra) - (qa-ra^2)) / (2*rho*(ra-1))
    # When ra -> 1: numerator -> (rho-1)*(qa-1), denominator -> 0
    # For exponential arrivals qa = ra = 1, so 0/0 form.
    # Use L'Hopital: d/d(ra)[num] = -rho + 2*ra, d/d(ra)[den] = 2*rho
    # At ra=1: theta_limit = (-rho + 2) / (2*rho) = (2 - rho) / (2*rho)
    if abs(ra - 1) < 1e-12:
        theta = (2 - rho) / (2 * rho)
    else:
        theta = (rho * (qa - ra) - (qa - ra ** 2)) / (2 * rho * (ra - 1))

    d = (1 + 1 / ra) * (1 - rs) * (1 - (q0 / qa) ** 3) * (1 - rho ** 3)
    D = (rs - theta) ** 2 + (2 * rs - 1 + d) * (ra - 1)

    # Myskja2 formula
    # Ensure D is non-negative (can become slightly negative due to numerics)
    D = max(D, 0.0)
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

def qsys_gg1(
    lambda_val: float, mu: float, ca2: float, cs2: float
) -> Tuple[float, float]:
    """
    G/G/1 queue analysis using exact methods for special cases and
    Allen-Cunneen approximation for the general case.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca2: Squared coefficient of variation of inter-arrival time
        cs2: Squared coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Mean response time (time in system)
            rhohat: Effective utilization

    References:
        Original JAR: jar/src/main/kotlin/jline/api/qsys/Qsys_gg1.kt
    """
    from .basic import qsys_mm1, qsys_mg1

    tol = 1e-8

    if abs(ca2 - 1.0) < tol and abs(cs2 - 1.0) < tol:
        # M/M/1 case
        result = qsys_mm1(lambda_val, mu)
        W = result['W']
        rhohat = result['rho']
    elif abs(ca2 - 1.0) < tol:
        # M/G/1 case (ca2 = 1)
        result = qsys_mg1(lambda_val, mu, np.sqrt(cs2))
        W = result['W']
        rhohat = result.get('rhohat', result.get('rho', lambda_val / mu))
    elif abs(cs2 - 1.0) < tol:
        # G/M/1 case (cs2 = 1) - matches MATLAB: W = qsys_gm1(rho, mu)
        rho = lambda_val / mu
        W = 1.0 / ((1 - rho) * mu)
        rhohat = W * lambda_val / (1 + W * lambda_val)
    else:
        # General case - Allen-Cunneen approximation
        W, rhohat = qsys_gig1_approx_allencunneen(lambda_val, mu, np.sqrt(ca2), np.sqrt(cs2))

    return W, rhohat

def qsys_gig1_lbnd(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Tuple[float, float]:
    """
    Fundamental theoretical lower bounds for G/G/1 queues.

    These are the minimum possible values that performance measures
    cannot fall below for any realization of the arrival and service
    processes.

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Lower bound on mean response time (= 1/mu)
            rhohat: Effective utilization

    References:
        Original JAR: jar/src/main/kotlin/jline/api/qsys/Qsys_gig1_lbnd.kt
    """
    W = 1.0 / mu  # At least the mean service time
    rhohat = W * lambda_val / (1 + W * lambda_val)
    return W, rhohat

def qsys_gigk_approx_cosmetatos(
    lambda_val: float, mu: float, ca: float, cs: float, k: int
) -> Tuple[float, float]:
    """
    G/G/k queue approximation using the Cosmetatos method.

    Adjusts M/M/k results based on the variability of arrival
    and service processes.

    Args:
        lambda_val: Arrival rate
        mu: Service rate per server
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time
        k: Number of servers

    Returns:
        Tuple of (W, rhohat) where:
            W: Approximate mean response time
            rhohat: Effective utilization

    References:
        Original JAR: jar/src/main/kotlin/jline/api/qsys/Qsys_gigk_approx_cosmetatos.kt
    """
    from .basic import qsys_mmk

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Get M/M/k baseline
    mmk_result = qsys_mmk(lambda_val, mu, k)
    W_mmk = mmk_result['W']

    # Cosmetatos variability adjustment
    Wq_mmk = W_mmk - 1.0 / mu
    variability_factor = (ca2 + cs2) / 2.0
    Wq = Wq_mmk * variability_factor
    W = Wq + 1.0 / mu

    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat

def qsys_gigk_approx_whitt(
    lambda_val: float, mu: float, ca: float, cs: float, k: int
) -> Tuple[float, float]:
    """
    G/G/k queue approximation using Whitt's QNA method.

    Uses the Queue Network Analyzer (QNA) methodology to provide
    accurate estimates for multi-server queues with general distributions.

    Args:
        lambda_val: Arrival rate
        mu: Service rate per server
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time
        k: Number of servers

    Returns:
        Tuple of (W, rhohat) where:
            W: Approximate mean response time
            rhohat: Effective utilization

    References:
        Original JAR: jar/src/main/kotlin/jline/api/qsys/Qsys_gigk_approx_whitt.kt
    """
    from .basic import qsys_mmk

    ca2 = ca ** 2
    cs2 = cs ** 2
    rho = lambda_val / (k * mu)

    if rho >= 1.0:
        return np.inf, 1.0

    # Get M/M/k baseline
    mmk_result = qsys_mmk(lambda_val, mu, k)
    W_mmk = mmk_result['W']
    Wq_mmk = W_mmk - 1.0 / mu

    # Heavy traffic factor
    if rho > 0.7:
        ht_factor = 1.0 + (4.0 * (rho - 0.7)) ** 2
    else:
        ht_factor = 1.0

    # Variability correction factor (Whitt's g function)
    if ca2 >= 1.0 and cs2 >= 1.0:
        g = 1.0
    elif ca2 <= 1.0 and cs2 <= 1.0:
        phi = (1.0 - ca2) * (1.0 - cs2) / (1.0 + cs2)
        g = 1.0 - phi * np.sqrt(k)
    else:
        g = (ca2 + cs2) / 2.0

    # QNA formula
    Wq = Wq_mmk * ((ca2 + cs2) / 2.0) * g * ht_factor
    Wq = max(Wq, 0.0)

    W = Wq + 1.0 / mu
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat
