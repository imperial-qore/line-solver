"""
Basic single-queue system analysis algorithms.

Native Python implementations for M/M/1, M/M/k, M/G/1, G/M/1, and M/G/inf queues.
"""

import numpy as np
from math import factorial
from typing import Dict, Optional, Any


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


def qsys_mmck(lambda_val: float, mu: float, c: int, K: int) -> Dict[str, float]:
    """
    Exact closed-form analysis of an M/M/c/K queue (finite capacity K, c servers).

    Port of MATLAB qsys_mmck.m. Stationary distribution (truncated Erlang form):
        a = lambda/mu, rho = a/c
        p_n = a^n/n! * p0                 for 0 <= n <= c
        p_n = a^c/c! * rho^(n-c) * p0     for c <= n <= K
    with p0 normalizing the (K+1)-point distribution.

    Args:
        lambda_val: Poisson arrival rate (> 0)
        mu: Per-server exponential service rate (> 0)
        c: Number of servers (>= 1)
        K: System capacity, total jobs allowed (K >= c)

    Returns:
        dict with L, Lq, W, Wq, rho, P0 plus MATLAB-style aliases
        (meanQueueLength, meanQueueLengthQ, meanWaitingTime, meanSojournTime,
        utilization, throughput, lossProbability, queueLengthDist).
    """
    import math
    c = int(c)
    K = int(K)
    a = lambda_val / mu
    rho = a / c

    p = np.zeros(K + 1)
    for n in range(0, c):
        p[n] = a ** n / math.factorial(n)
    ac_over_cfact = a ** c / math.factorial(c)
    for n in range(c, K + 1):
        p[n] = ac_over_cfact * rho ** (n - c)
    S = float(np.sum(p))
    if not np.isfinite(S) or S <= 0:
        raise ValueError("qsys_mmck: stationary distribution failed to normalize")
    p = p / S

    levels = np.arange(K + 1)
    L = float(levels @ p)
    p_K = float(p[K])
    lambda_eff = lambda_val * (1.0 - p_K)
    n_waiting = np.maximum(0, levels - c)
    Lq = float(n_waiting @ p)
    util = lambda_eff / (c * mu)
    if lambda_eff > 0:
        Wq = Lq / lambda_eff
        W = L / lambda_eff
    else:
        Wq = 0.0
        W = 0.0

    return {
        'L': L,
        'Lq': Lq,
        'W': W,
        'Wq': Wq,
        'rho': util,
        'P0': float(p[0]),
        'meanQueueLength': L,
        'meanQueueLengthQ': Lq,
        'meanWaitingTime': Wq,
        'meanSojournTime': W,
        'utilization': util,
        'throughput': lambda_eff,
        'lossProbability': p_K,
        'queueLengthDist': p,
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

    # see _kb/03-api-layer.md for rationale
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

def qsys_gig1_rq(rho: float, mu: float, cs2: float, IaFun):
    """
    Robust Queueing (RQ) approximation for a single G/GI/1 queue partially
    characterized by its arrival rate, index of dispersion for counts (IDC) and
    the first two moments of the service time. Implements the mean steady-state
    workload

        Z* = sup_{x>=0} { -(1-rho) x + sqrt( 2 rho x (I_a(x) + c2_s) / mu ) }

    and the derived steady-state performance measures.

    Reference:
        W. Whitt and W. You (2018), "A Robust Queueing Network Analyzer Based on
        Indices of Dispersion", eqs. (13),(16)-(18).

    Args:
        rho: Traffic intensity lambda/mu (0<rho<1)
        mu: Service rate
        cs2: Service SCV c2_s
        IaFun: callable, IaFun(x) -> arrival IDC I_a(x) at time argument x>0

    Returns:
        Tuple (Z, W, Q, X) with mean workload E[Z], waiting time E[W], queue
        length E[Q] (waiting + in service), and number in system E[X].
    """
    if rho <= 0:
        return 0.0, 0.0, 0.0, 0.0
    if rho >= 1:
        return np.inf, np.inf, np.inf, np.inf
    lambda_val = rho * mu

    def negf(x):
        if x <= 0:
            return 0.0
        ia = float(IaFun(x))
        return -(-(1 - rho) * x + np.sqrt(max(0.0, 2 * rho * x * (ia + cs2) / mu)))

    # see _kb/03-api-layer.md for rationale
    xs = np.logspace(-6, 8, 200)
    fv = np.array([-negf(x) for x in xs])
    imax = int(np.argmax(fv))
    lo = xs[max(0, imax - 1)]
    hi = xs[min(len(xs) - 1, imax + 1)]
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(negf, bounds=(lo, hi), method='bounded',
                          options={'xatol': 1e-10})
    Z = max(fv[imax], -res.fun)
    Z = max(Z, 0.0)

    # derived measures (eqs. 16-18)
    W = max(0.0, Z / rho - (cs2 + 1) / (2 * mu))
    Q = lambda_val * W          # E[Q] waiting (Little's law on waiting time)
    X = Q + rho                 # E[X] number in system including one in service
    return Z, W, Q, X


def qsys_gm1(sigma: float, mu: float) -> Dict[str, float]:
    """
    Analyze G/M/1 queue (general arrivals, exponential service).

    Matches MATLAB qsys_gm1(sigma, mu) and JAR Qsys_gm1: the number of
    customers found by an arrival is geometric with parameter sigma, so the
    mean response time (time in system) is W = 1/(mu*(1-sigma)).

    Args:
        sigma: Root in (0,1) of sigma = A*(mu*(1-sigma)), where A* is the
            Laplace-Stieltjes transform of the interarrival-time
            distribution.
        mu: Service rate

    Returns:
        dict: {'W': mean response time}

    Note:
        To obtain sigma from the first two moments of the interarrival time,
        use qsys_gg1(lambda_val, mu, ca2, 1.0), which fits a two-moment
        renewal process and solves the fixed point.
    """
    W = 1.0 / ((1.0 - sigma) * mu)
    return {'W': W}

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


def qsys_mmcc_retrial_fp(lambda_val: float, mu: float, c: int,
                          tol: float = 1e-10, maxiter: int = 10000) -> Dict[str, float]:
    """
    Fixed-point approximation for M/M/c/c retrial queue.

    Customers arrive at rate lambda to a system with c servers (no waiting
    room), each with service rate mu. Blocked customers join an orbit and
    retry. Under the assumption that the retrial rate is small relative to
    the service rate, the total arrival flow (fresh + retrial) is
    approximated by a Poisson process with rate lambda + r, where r
    satisfies the fixed-point equation:

        r = (lambda + r) * B((lambda + r) / mu, c)

    and B(a, c) is the Erlang-B blocking probability for offered load a
    and c servers.

    Args:
        lambda_val: Arrival rate
        mu: Service rate per server
        c: Number of servers (= capacity, no waiting room)
        tol: Convergence tolerance (default: 1e-10)
        maxiter: Maximum iterations (default: 10000)

    Returns:
        dict: Performance measures including:
            - blocProb: Blocking probability
            - r: Additional arrival rate due to retrials
            - niter: Number of iterations to converge
            - rho: Offered load (lambda / (c * mu))
            - L: Mean number of busy servers

    References:
        Cohen (1957), fixed-point approximation for M/M/c/c retrial queues.
        Phung-Duc, "Retrial Queueing Models: A Survey on Theory and
        Applications", 2019, Eq. (1).

    Example:
        >>> result = qsys_mmcc_retrial_fp(2.0, 1.0, 3)
        >>> print(f"Blocking: {result['blocProb']:.4f}")
    """
    r = 0.0
    niter = 0
    for i in range(1, maxiter + 1):
        niter = i
        a = (lambda_val + r) / mu  # offered load
        b = _erlang_b(c, a)
        r_new = (lambda_val + r) * b
        if abs(r_new - r) < tol:
            r = r_new
            break
        r = r_new

    a_eff = (lambda_val + r) / mu
    blocProb = _erlang_b(c, a_eff)
    rho = lambda_val / (c * mu)
    L = a_eff * (1 - blocProb)

    return {
        'blocProb': blocProb,
        'r': r,
        'niter': niter,
        'rho': rho,
        'L': L,
    }
