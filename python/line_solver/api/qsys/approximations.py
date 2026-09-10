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

    Note: alias of qsys_gig1_ubnd_kingman ('gig1.kingman' method).

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Upper bound on mean response time (time in system)
            rhohat: Effective utilization (so that M/M/1 formulas still hold)
    """
    return qsys_gig1_ubnd_kingman(lambda_val, mu, ca, cs)

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
) -> Tuple[float, float]:
    """
    Gelenbe's diffusion approximation for G/G/1 with instantaneous-return
    boundary::

        p(0) = 1-rho,  p(n) = rho*(1-rhat)*rhat^(n-1), n>=1
        rhat = exp(-2*(1-rho)/(rho*ca^2+cs^2))

    hence E[N] = rho/(1-rhat) and the mean response time (time in system)
    is W = E[N]/lambda = 1/(mu*(1-rhat)).

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Mean response time (time in system)
            rhohat: Effective utilization (so that M/M/1 formulas still hold)

    References:
        Gelenbe, E. (1975). "On approximate computer system models".
        Journal of the ACM 22(2), 261-269.
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf, 1.0

    ca2 = ca ** 2
    cs2 = cs ** 2

    rhat = np.exp(-2 * (1 - rho) / (rho * ca2 + cs2))
    W = 1.0 / (mu * (1.0 - rhat))
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat

def qsys_gig1_approx_kimura(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Tuple[float, float]:
    """
    Kimura's diffusion-interpolation approximation for G/G/1::

        Wq = rho*(ca^2+cs^2)/(mu*(1-rho)*(1+ca^2))

    exact for M/M/1 and M/G/1. The returned W adds the mean service time
    (response time, time in system).

    Args:
        lambda_val: Arrival rate
        mu: Service rate
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time

    Returns:
        Tuple of (W, rhohat) where:
            W: Mean response time (time in system)
            rhohat: Effective utilization (so that M/M/1 formulas still hold)

    References:
        Kimura, T. (1986). "A two-moment approximation for the mean waiting
        time in the GI/G/s queue". Management Science 32(6), 751-763.
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf, 1.0

    ca2 = ca ** 2
    cs2 = cs ** 2

    Wq = rho * (ca2 + cs2) / mu / (1.0 - rho) / (1.0 + ca2)
    W = Wq + 1.0 / mu
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat

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

    # rho >= 1 is NOT short-circuited to Inf here: the reference (and the JAR and
    # C++ ports) evaluate the formula unconditionally, so an overloaded queue
    # returns the negative W the expression yields, which the caller's saturation
    # rule then reports as Inf QLen/RespT and zero residence time. Substituting
    # Inf makes ResidT = Inf * V instead, which is a different reported model.
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
) -> Tuple[float, float]:
    """
    Myskja's third-moment approximation for G/G/1::

        Wq = rho/(2*mu*(1-rho))*((1+cs^2)+(q0/qa)^(1/rho-rho)*(1/rho)*(ca^2-1))

    exact for M/G/1 (ca=1). The returned W adds the mean service time
    (response time, time in system).

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
        Tuple of (W, rhohat) where:
            W: Mean response time (time in system)
            rhohat: Effective utilization (so that M/M/1 formulas still hold)
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf, 1.0

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Myskja formula incorporating third moments
    Wq = (rho / (2 * mu * (1 - rho))) * (
        (1 + cs2) + (q0 / qa) ** (1 / rho - rho) * (1 / rho) * (ca2 - 1)
    )
    W = Wq + 1.0 / mu
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat

def qsys_gig1_approx_myskja2(
    lambda_val: float, mu: float, ca: float, cs: float,
    q0: float, qa: float
) -> Tuple[float, float]:
    """
    Modified Myskja (Myskja2) approximation for G/G/1, returning the mean
    response time (time in system). For ca=1 the interpolation parameter
    theta is a 0/0 form, so the exact M/G/1 result is returned instead
    (also the interpolation anchor of the method).

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
        Tuple of (W, rhohat) where:
            W: Mean response time (time in system)
            rhohat: Effective utilization (so that M/M/1 formulas still hold)
    """
    from .basic import qsys_mg1

    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf, 1.0

    ca2 = ca ** 2
    cs2 = cs ** 2

    if abs(ca2 - 1) < 1e-8:
        # M/G/1 case: exact (also the interpolation anchor of the method)
        result = qsys_mg1(lambda_val, mu, cs)
        W = result['W']
        rhohat = W * lambda_val / (1 + W * lambda_val)
        return W, rhohat

    # Intermediate calculations
    ra = (1 + ca2) / 2
    rs = (1 + cs2) / 2

    theta = (rho * (qa - ra) - (qa - ra ** 2)) / (2 * rho * (ra - 1))
    d = (1 + 1 / ra) * (1 - rs) * (1 - (q0 / qa) ** 3) * (1 - rho ** 3)
    D = (rs - theta) ** 2 + (2 * rs - 1 + d) * (ra - 1)

    # Myskja2 formula
    # Ensure D is non-negative (can become slightly negative due to numerics)
    D = max(D, 0.0)
    W = (rho / (1 - rho)) / lambda_val * (rs + (1 / rho) * (np.sqrt(D) - (rs - theta)))
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat

def qsys_gig1_ubnd_kingman(
    lambda_val: float, mu: float, ca: float, cs: float
) -> Tuple[float, float]:
    """
    Kingman's upper bound on the mean waiting time of a G/G/1 queue::

        Wq <= lambda*(sa^2+ss^2)/(2*(1-rho)),

    with sa^2=ca^2/lambda^2 and ss^2=cs^2/mu^2. The returned W adds the
    mean service time, so it upper-bounds the mean response time
    (time in system).

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
        Kingman, J.F.C. (1962). "Some inequalities for the queue GI/G/1".
        Biometrika 49(3/4), 315-324.
        Original MATLAB: matlab/src/api/qsys/qsys_gig1_ubnd_kingman.m
    """
    rho = lambda_val / mu

    if rho >= 1.0:
        return np.inf, 1.0

    ca2 = ca ** 2
    cs2 = cs ** 2

    # Kingman's upper bound formula
    Wq = lambda_val * (ca2 / lambda_val ** 2 + cs2 / mu ** 2) / (2 * (1 - rho))
    W = Wq + 1 / mu
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
    Allen-Cunneen approximation for the general case. In the G/M/1 case,
    the interarrival-time distribution is fitted from (lambda, ca2) by a
    two-moment renewal process (H2 with balanced means for ca2>1, mixed
    Erlang for ca2<1) and sigma is the root of sigma = A*(mu*(1-sigma)),
    with A* the interarrival-time LST.

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
        Original MATLAB: matlab/src/api/qsys/qsys_gg1.m
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
        rhohat = W * lambda_val / (1 + W * lambda_val)
    elif abs(cs2 - 1.0) < tol:
        # G/M/1 case (cs2 = 1)
        sigma = _qsys_gm1_sigma(lambda_val, mu, ca2)
        W = 1.0 / ((1 - sigma) * mu)
        rhohat = W * lambda_val / (1 + W * lambda_val)
    else:
        # General case - Allen-Cunneen approximation
        W, rhohat = qsys_gig1_approx_allencunneen(lambda_val, mu, np.sqrt(ca2), np.sqrt(cs2))

    return W, rhohat

def _qsys_gm1_sigma(lambda_val: float, mu: float, ca2: float) -> float:
    """
    Root in (0,1) of sigma = A*(mu*(1-sigma)) for a two-moment fit of the
    interarrival-time LST A*. The map T(x)=A*(mu*(1-x)) is increasing with
    the queue root as its smallest fixed point, so fixed-point iterates
    converge monotonically.
    """
    if ca2 < 1e-6:
        # deterministic interarrival times
        def lst(s):
            return np.exp(-s / lambda_val)
    elif ca2 < 1.0:
        # mixed Erlang(j-1,j) with common rate (Tijms, 1994)
        j = int(np.ceil(1.0 / ca2))
        p = (j * ca2 - np.sqrt(j * (1 + ca2) - j ** 2 * ca2)) / (1 + ca2)
        nu = (j - p) * lambda_val

        def lst(s):
            return p * (nu / (s + nu)) ** (j - 1) + (1 - p) * (nu / (s + nu)) ** j
    else:
        # hyperexponential H2 with balanced means
        p1 = (1 + np.sqrt((ca2 - 1) / (ca2 + 1))) / 2
        l1 = 2 * p1 * lambda_val
        l2 = 2 * (1 - p1) * lambda_val

        def lst(s):
            return p1 * l1 / (s + l1) + (1 - p1) * l2 / (s + l2)

    sigma = lambda_val / mu
    for _ in range(100000):
        signew = lst(mu * (1.0 - sigma))
        if abs(signew - sigma) < 1e-13:
            return signew
        sigma = signew
    return sigma

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
        Original JAR: jar/src/main/java/jline/api/qsys/Qsys_gig1_lbnd.java
    """
    W = 1.0 / mu  # At least the mean service time
    rhohat = W * lambda_val / (1 + W * lambda_val)
    return W, rhohat

def qsys_gigk_approx_cosmetatos(
    lambda_val: float, mu: float, ca: float, cs: float, k: int
) -> Tuple[float, float]:
    """
    GI/G/k approximation by interpolation of the M/M/k, M/D/k and D/M/k
    queues (Cosmetatos 1982; Page 1982)::

        Wq = [ca^2*cs^2 + ca^2*(1-cs^2)*phi1/2
              + (1-ca^2)*cs^2*phi3/2] * Wq(M/M/k)

    where phi1 and phi3 are the Cosmetatos (1975) correction factors for
    M/D/k and D/M/k, with the safeguards of Whitt (1993). The D/D/k corner
    has Wq=0. The interpolation requires ca^2<=1 and cs^2<=1; outside this
    region the Lee-Longton scaling Wq = ((ca^2+cs^2)/2)*Wq(M/M/k) is used.

    Args:
        lambda_val: Arrival rate
        mu: Service rate per server
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time
        k: Number of servers

    Returns:
        Tuple of (W, rhohat) where:
            W: Approximate mean response time (time in system)
            rhohat: Effective utilization

    References:
        Cosmetatos, G.P. (1975). "Approximate explicit formulae for the
        average queueing time in the processes (M/D/r) and (D/M/r)".
        INFOR 13, 328-331.
        Page, E. (1982). "Tables of waiting times for M/M/n, M/D/n and D/M/n
        and their use to give approximate waiting times in more general
        queues". J. Opl. Res. Soc. 33, 453-473.
    """
    from .basic import qsys_mmk

    ca2 = ca ** 2
    cs2 = cs ** 2
    rho = lambda_val / (k * mu)

    if rho >= 1.0:
        return np.inf, 1.0

    # Exact M/M/k baseline waiting time (Erlang-C based)
    mmk_result = qsys_mmk(lambda_val, mu, k)
    W_mmk = mmk_result['W']
    Wq_mmk = W_mmk - 1.0 / mu

    if ca2 <= 1 and cs2 <= 1:
        # Cosmetatos correction, as modified by Whitt (1993), eq. (2.17)
        gamma = min(0.24, (1 - rho) * (k - 1) * (np.sqrt(4 + 5 * k) - 2) / (16 * k * rho))
        phi1 = 1 + gamma  # M/D/k factor
        phi3 = (1 - 4 * gamma) * np.exp(-2 * (1 - rho) / (3 * rho))  # D/M/k factor
        Wq = (ca2 * cs2 + ca2 * (1 - cs2) * phi1 / 2 + (1 - ca2) * cs2 * phi3 / 2) * Wq_mmk
    else:
        # Interpolation weights are invalid outside the unit box
        Wq = ((ca2 + cs2) / 2) * Wq_mmk
    W = Wq + 1.0 / mu

    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat

def qsys_gigk_approx_whitt(
    lambda_val: float, mu: float, ca: float, cs: float, k: int
) -> Tuple[float, float]:
    """
    GI/G/k approximation of Whitt (1993), eqs. (2.16)-(2.25)::

        Wq = phi(rho,ca^2,cs^2,k) * ((ca^2+cs^2)/2) * Wq(M/M/k)

    where phi interpolates the Cosmetatos M/D/k (phi1) and D/M/k (phi3)
    correction factors. Exact for M/M/k; reduces to the Cosmetatos M/D/k
    approximation for cs=0. Implements eq. (2.25) as printed, which was
    validated against the paper's Tables 5-7 (New column).

    Args:
        lambda_val: Arrival rate
        mu: Service rate per server
        ca: Coefficient of variation of inter-arrival time
        cs: Coefficient of variation of service time
        k: Number of servers

    Returns:
        Tuple of (W, rhohat) where:
            W: Approximate mean response time (time in system)
            rhohat: Effective utilization

    References:
        Whitt, W. (1993). "Approximations for the GI/G/m queue".
        Production and Operations Management 2(2), 114-161.
    """
    from .basic import qsys_mmk

    ca2 = ca ** 2
    cs2 = cs ** 2
    rho = lambda_val / (k * mu)

    if rho >= 1.0:
        return np.inf, 1.0

    # Exact M/M/k baseline (Erlang-C based)
    mmk_result = qsys_mmk(lambda_val, mu, k)
    W_mmk = mmk_result['W']
    Wq_mmk = W_mmk - 1.0 / mu

    # Cosmetatos correction, as modified by Whitt (1993), eq. (2.17)
    gamma = min(0.24, (1 - rho) * (k - 1) * (np.sqrt(4 + 5 * k) - 2) / (16 * k * rho))
    phi1 = 1 + gamma                             # M/D/k factor, eq. (2.16)
    phi2 = 1 - 4 * gamma                         # eq. (2.18)
    phi3 = phi2 * np.exp(-2 * (1 - rho) / (3 * rho))  # D/M/k factor, eq. (2.20)
    phi4 = min(1.0, (phi1 + phi3) / 2)           # eq. (2.21)

    c2 = (ca2 + cs2) / 2
    if c2 >= 1:
        psi = 1.0                                # eq. (2.22)
    else:
        psi = phi4 ** (2 * (1 - c2))

    if abs(ca2 - cs2) < 1e-12:
        phi = psi                                # eq. (2.25) reduces to psi
    elif ca2 > cs2:
        phi = (4 * (ca2 - cs2) / (4 * ca2 - 3 * cs2)) * phi1 + (cs2 / (4 * ca2 - 3 * cs2)) * psi
    else:
        phi = ((cs2 - ca2) / (2 * (ca2 + cs2))) * phi3 + ((cs2 + 3 * ca2) / (2 * (ca2 + cs2))) * psi

    Wq = phi * c2 * Wq_mmk                       # eq. (2.24)
    W = Wq + 1.0 / mu
    rhohat = W * lambda_val / (1 + W * lambda_val)

    return W, rhohat
