"""
Halfin-Whitt QED (quality-and-efficiency-driven) approximations for M/M/s.

Native Python twin of matlab/src/api/qsys/qsys_mmk_qed.m,
qsys_mmk_qed_alpha.m and qsys_mmk_qed_staffing.m: the many-server heavy-traffic
limit of S. Halfin and W. Whitt (1981), Operations Research 29(3), 567-588, and
the square-root staffing rule that inverts it.
"""

from math import ceil, erfc, exp, floor, pi, sqrt
from typing import Any, Callable, Dict, Union

import numpy as np


def qsys_mmk_qed_alpha(beta: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """
    The Halfin-Whitt delay-probability function
    ``alpha(beta) = [1 + beta*Phi(beta)/phi(beta)]^-1`` for ``beta > 0``, with
    ``phi`` and ``Phi`` the standard normal density and cdf.

    It is the limit of the Erlang C delay probability of the M/M/s queue as
    ``s -> inf`` with ``beta = (1-rho)sqrt(s)`` held fixed, decreasing strictly
    from 1 at ``beta = 0`` to 0 as ``beta -> inf``, which is what makes it
    invertible for staffing. Non-positive ``beta`` returns 1: with no server
    slack every arrival is delayed.

    Evaluated as ``phi/(phi + beta*Phi)`` rather than as the reciprocal of
    ``1 + beta*Phi/phi``: the two are the same function, but the quotient
    ``Phi/phi`` overflows once ``phi`` underflows (``beta`` beyond about 38),
    whereas this form degrades to ``0/(0+beta) = 0``, the correct limit.

    Args:
        beta: the QED server-slack parameter, scalar or array

    Returns:
        alpha(beta), of the same shape as the input.

    References:
        S. Halfin, W. Whitt (1981). Heavy-traffic limits for queues with many
        exponential servers. Operations Research 29(3), 567-588.
    """
    b = np.asarray(beta, dtype=float)
    out = np.ones(b.shape) if b.ndim else np.array(1.0)
    pos = b > 0
    if np.any(pos):
        bp = b[pos] if b.ndim else b
        phi = np.exp(-bp ** 2 / 2.0) / np.sqrt(2.0 * np.pi)
        Phi = np.array([erfc(-float(v) / sqrt(2.0)) / 2.0 for v in np.atleast_1d(bp)])
        val = phi / (phi + bp * Phi.reshape(np.shape(phi)))
        if b.ndim:
            out[pos] = val
        else:
            out = np.array(float(val))
    return float(out) if not b.ndim else out


def qsys_mmk_qed(lambda_val: float, mu: float, s: int) -> Dict[str, Any]:
    """
    Halfin-Whitt QED approximation for the M/M/s queue.

    Let ``s`` grow with the offered load ``a = lambda/mu`` so that the server
    slack ``beta = (1-rho)sqrt(s) = (s-a)/sqrt(s)`` stays fixed. The delay
    probability then has the non-degenerate limit ``alpha(beta)``: servers are
    busy a fraction ``1 - beta/sqrt(s)`` of the time, so efficiency tends to 1,
    and yet the delay probability tends to a constant strictly between 0 and 1,
    so quality does not collapse.

    Useful even though M/M/s is exactly solvable, because Erlang C needs a sum of
    ``s`` terms ``a^j/j!`` that overflows in double precision well before the
    thousands of servers a large contact centre or thread pool has.

    Args:
        lambda_val: arrival rate
        mu: service rate of one server
        s: number of servers

    Returns:
        Dict with ``offeredLoad``, ``trafficIntensity``, ``beta``, ``probDelay``,
        ``meanWaitDelayed``, ``meanWait``, ``meanQueueLength``, ``meanNumber``
        and ``utilization``. An overloaded model (``beta <= 0``) has no QED
        limit: ``probDelay`` is 1 and the waiting-time fields are infinite.

    References:
        S. Halfin, W. Whitt (1981). Heavy-traffic limits for queues with many
        exponential servers. Operations Research 29(3), 567-588.
    """
    if lambda_val <= 0:
        raise ValueError('The arrival rate lambda must be positive.')
    if mu <= 0:
        raise ValueError('The service rate mu must be positive.')
    s = int(round(s))
    if s < 1:
        raise ValueError('The number of servers s must be at least 1.')

    a = lambda_val / mu
    rho = a / s
    beta = (s - a) / sqrt(s)
    result: Dict[str, Any] = {'offeredLoad': a, 'trafficIntensity': rho, 'beta': beta,
                              'utilization': rho}
    if beta <= 0:
        result.update({'probDelay': 1.0, 'meanWaitDelayed': float('inf'),
                       'meanWait': float('inf'), 'meanQueueLength': float('inf'),
                       'meanNumber': float('inf')})
        return result

    result['probDelay'] = float(qsys_mmk_qed_alpha(beta))
    result['meanWaitDelayed'] = 1.0 / (s * mu - lambda_val)
    result['meanWait'] = result['probDelay'] * result['meanWaitDelayed']
    result['meanQueueLength'] = lambda_val * result['meanWait']
    result['meanNumber'] = a + result['meanQueueLength']
    return result


def _solve(f: Callable[[float], float]) -> float:
    """
    Bisection for a root of an increasing f on (0, hi]; the bracket grows until
    the sign changes.
    """
    lo, hi = 1e-9, 1.0
    if f(lo) > 0:
        return lo
    while f(hi) < 0:
        hi *= 2.0
        if hi > 1e6:
            raise ValueError('no server slack meets the target; the target is unattainable')
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _erlang_c_stable(s: int, lambda_val: float, mu: float) -> float:
    """
    Erlang C by the recursion B_j = a B_{j-1}/(j + a B_{j-1}) on the Erlang B
    blocking probability, which never forms a^j/j! and so never overflows.
    """
    a = lambda_val / mu
    b = 1.0
    for j in range(1, s + 1):
        b = a * b / (j + a * b)
    rho = a / s
    if rho >= 1.0:
        return 1.0
    return b / (1.0 - rho * (1.0 - b))


def qsys_mmk_qed_staffing(lambda_val: float, mu: float, target: Any,
                          criterion: str = 'delay', exact: bool = False,
                          maxServers: int = 10 ** 7) -> Dict[str, Any]:
    """
    Square-root staffing of the M/M/s queue.

    Invert ``alpha(beta) = target`` for the server slack and staff
    ``s = ceil(a + beta*sqrt(a))`` with ``a = lambda/mu``: the base ``a`` erlangs
    of work plus a cushion that grows only as the square root of the load.
    Doubling the load needs only ``sqrt(2)`` times the cushion, which is why
    large service systems can be both highly utilized and responsive.

    Args:
        lambda_val: arrival rate
        mu: service rate of one server
        target: the target, read according to ``criterion``: a probability for
            ``'delay'``, a time for ``'meanwait'``, or a dict with keys
            ``deadline`` and ``level`` for ``'servicelevel'``
        criterion: ``'delay'`` (P(W>0) <= target), ``'meanwait'`` (E[W] <= target)
            or ``'servicelevel'`` (P(W <= deadline) >= level)
        exact: walk ``s`` until the EXACT Erlang C measure meets the target,
            starting from the square-root answer
        maxServers: cap on that walk

    Returns:
        Dict with ``numServers``, ``beta``, ``betaTarget``, ``offeredLoad``,
        ``probDelay``, ``meanWait``, ``exactUsed`` and, for the service-level
        criterion, ``serviceLevel``.

    References:
        S. Halfin, W. Whitt (1981). Heavy-traffic limits for queues with many
        exponential servers. Operations Research 29(3), 567-588. The staffing
        form is the standard reading of that limit; see also W. Whitt (2007),
        Naval Research Logistics 54(5), 476-484.
    """
    if lambda_val <= 0:
        raise ValueError('The arrival rate lambda must be positive.')
    if mu <= 0:
        raise ValueError('The service rate mu must be positive.')
    a = lambda_val / mu
    crit = criterion.lower()

    if crit == 'delay':
        if not np.isscalar(target) or not (0 < float(target) < 1):
            raise ValueError('For the delay criterion the target must be in (0,1).')
        beta_target = _solve(lambda b: float(target) - float(qsys_mmk_qed_alpha(b)))
    elif crit == 'meanwait':
        if not np.isscalar(target) or float(target) <= 0:
            raise ValueError('For the meanwait criterion the target must be positive.')
        # E[W] = alpha(beta)/(mu beta sqrt(a)) at s ~ a + beta sqrt(a); the
        # residual is written target - E[W] so that it increases in beta.
        beta_target = _solve(
            lambda b: float(target) - float(qsys_mmk_qed_alpha(b)) / (mu * b * sqrt(a)))
    elif crit == 'servicelevel':
        if not isinstance(target, dict) or 'deadline' not in target or 'level' not in target:
            raise ValueError('For the servicelevel criterion the target must be a dict with keys '
                             'deadline and level.')
        if not (0 < target['level'] < 1) or target['deadline'] <= 0:
            raise ValueError('The service level must be in (0,1) and the deadline positive.')
        # P(W > t) = alpha(beta) exp(-(s mu - lambda) t), s mu - lambda = mu beta sqrt(s).
        beta_target = _solve(lambda b: (1.0 - float(qsys_mmk_qed_alpha(b)) *
                                        exp(-mu * b * sqrt(a + b * sqrt(a)) * target['deadline']))
                             - target['level'])
    else:
        raise ValueError('unknown criterion %s' % criterion)

    s = max(1, int(ceil(a + beta_target * sqrt(a))))
    if s * mu <= lambda_val:
        s = int(floor(a)) + 1

    exact_used = False
    if exact:
        exact_used = True
        while not _meets(lambda_val, mu, s, target, crit):
            s += 1
            if s > maxServers:
                raise ValueError('the exact refinement passed maxServers without meeting the target')
        while s > 1 and _meets(lambda_val, mu, s - 1, target, crit):
            s -= 1

    qed = qsys_mmk_qed(lambda_val, mu, s)
    result: Dict[str, Any] = {
        'numServers': s,
        'beta': qed['beta'],
        'betaTarget': beta_target,
        'offeredLoad': a,
        'probDelay': qed['probDelay'],
        'meanWait': qed['meanWait'],
        'exactUsed': exact_used,
    }
    if crit == 'servicelevel':
        result['serviceLevel'] = 1.0 - qed['probDelay'] * exp(
            -(s * mu - lambda_val) * target['deadline'])
    return result


def _meets(lambda_val: float, mu: float, s: int, target: Any, criterion: str) -> bool:
    """The exact M/M/s measure against the target."""
    if s * mu <= lambda_val:
        return False
    c = _erlang_c_stable(s, lambda_val, mu)
    wq = c / (s * mu - lambda_val)
    if criterion == 'delay':
        return c <= float(target)
    if criterion == 'meanwait':
        return wq <= float(target)
    if criterion == 'servicelevel':
        return (1.0 - c * exp(-(s * mu - lambda_val) * target['deadline'])) >= target['level']
    return False
