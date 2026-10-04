"""
State dependent Bernoulli server on a discrete time scale.

Native Python port of matlab/src/api/dqsys/dqsys_bernoulli1.m.

Reference: H. Daduna, Queueing Networks with Discrete Time Scale, LNCS 2046,
Springer 2001, theorem 2.3, corollaries 2.7 and 2.8, example 2.10 and
theorem 2.11.
"""

import math
from typing import Dict, Sequence, Union

import numpy as np

from .discrete import dqsys_geogeo1, LAS_DA

__all__ = ['dqsys_bernoulli1']

Number = Union[float, Sequence[float], np.ndarray]


def dqsys_bernoulli1(b: Number, p: Number, L: float = math.inf) -> Dict[str, object]:
    """Analyze a state dependent Bernoulli server.

    Args:
        b: per-slot arrival probability; a scalar b, or a vector with
            ``b[n]`` the probability in a slot that finds n jobs present
        p: per-slot service completion probability; a scalar p, or a vector
            with ``p[n-1] = p(n)`` for n = 1, 2, ...
        L: buffer capacity in jobs. ``math.inf`` (the default) leaves the
            buffer unbounded, which requires scalar b and p with b < p.

    Returns:
        dict with the fields listed below.

    Time advances in slots. In the slot starting at t with n jobs present the
    job in service departs with probability p(n) and an arrival occurs with
    probability b(n), independently; both are recorded at the end of the slot
    with the departure resolved first (Daduna's LA rule and D/A rule). The
    queue length at slot boundaries is a discrete birth-death chain with::

        pi(n) = [prod_{m=0}^{n-1} b(m) / prod_{m=0}^{n} c(m)]
              * [prod_{m=1}^{n-1} q(m) / prod_{m=1}^{n} p(m)] / H

    c = 1-b and q = 1-p, which is theorem 2.3, and corollary 2.8 once b(n) = 0
    above the capacity. For constant b and p it collapses to the Geo/Geo/1 law
    of :func:`dqsys_geogeo1` under the LAS_DA convention.

    The law seen by an arriving customer, with himself not counted, is
    theorem 2.11::

        pi_1(n) = [prod_{m=0}^{n} b(m) / prod_{m=0}^{n+1} c(m)]
                * [prod_{m=1}^{n} q(m) / prod_{m=1}^{n} p(m)] / H_1

    and is returned in ``arrivalPmf``. It is not the time-stationary law:
    discrete time has no PASTA analogue, and the two differ even when the
    arrival stream is a state independent Bernoulli process. In that state
    independent case pi_1 is exactly the EAS-convention queue length law of
    :func:`dqsys_geogeo1`, geometric with ratio r = b(1-p)/(p(1-b)).

    Fields:
        - capacity: buffer capacity, ``math.inf`` when unbounded
        - arrivalProb: offered per-slot arrival probability b(n)
        - serviceProb: per-slot service completion probability p(n)
        - pmf: time-stationary queue length law, ``pmf[n]`` for n = 0..L, or a
          callable on an unbounded buffer
        - arrivalPmf: arrival queue length law of theorem 2.11, same shape
        - emptyProb: pi(0)
        - utilization: fraction of slots with the server busy, 1 - pi(0)
        - throughput: carried departures per slot
        - lossProb: fraction of offered arrivals lost, 0 when unbounded
        - meanQueueLength, meanWaitingQueue: mean jobs in system / waiting
        - meanSojournTime, meanWaitingTime: in slots, by Little's law
        - normConst: normalizing constant H of theorem 2.3
        - analyzer: 'dqsys_bernoulli1'

    Examples:
        >>> r = dqsys_bernoulli1(0.2, 0.5)
        >>> round(r['meanQueueLength'], 6)
        0.533333
    """
    if L is None:
        L = math.inf
    if not (math.isinf(L) or (float(L).is_integer() and L >= 1)):
        raise ValueError('L must be a positive integer or math.inf')

    barr = np.atleast_1d(np.asarray(b, dtype=float))
    parr = np.atleast_1d(np.asarray(p, dtype=float))
    if barr.size == 0 or np.any(barr < 0) or np.any(barr > 1):
        raise ValueError('arrival probabilities must be real and in [0,1]')
    if parr.size == 0 or np.any(parr <= 0) or np.any(parr > 1):
        raise ValueError('service probabilities must be real and in (0,1]')

    if math.isinf(L):
        if barr.size != 1 or parr.size != 1:
            raise ValueError('an unbounded buffer requires scalar arrival and service '
                             'probabilities; pass a capacity L to analyze a state dependent server')
        return _bernoulli1_infinite(float(barr[0]), float(parr[0]))

    L = int(L)
    boff = _expand_state(barr, L + 1, 'arrival')      # boff[n] = b(n), n = 0..L
    pv = _expand_service(parr, L)                     # pv[n-1] = p(n), n = 1..L
    badm = boff.copy()
    badm[L] = 0.0                                     # an arrival finding L jobs is lost
    if np.any(badm[:L] >= 1):
        # c(n)=0 makes the weight of theorem 2.3 diverge at n; p(n)=1 is fine
        # and truncates the chain instead, which is example 2.9.
        raise ValueError('arrival probabilities below the capacity must be strictly less than one')

    lw = np.zeros(L + 1)                              # log unnormalized pi
    acc = -math.log(1.0 - badm[0])
    lw[0] = acc
    for n in range(1, L + 1):
        acc += _log(badm[n - 1]) - math.log(1.0 - badm[n]) - math.log(pv[n - 1])
        if n >= 2:
            acc += _log(1.0 - pv[n - 2])
        lw[n] = acc
    w = np.exp(lw - lw.max())
    H = float(w.sum())
    pmf = w / H

    la = np.full(L, -math.inf)                        # log unnormalized pi_1
    if L >= 1:
        acc = _log(badm[0]) - math.log(1.0 - badm[0]) - math.log(1.0 - badm[1])
        la[0] = acc
        for n in range(1, L):
            acc += _log(badm[n]) - math.log(1.0 - badm[n + 1]) \
                + _log(1.0 - pv[n - 1]) - math.log(pv[n - 1])
            la[n] = acc
    finite = np.isfinite(la)
    if not finite.any():
        arrivalPmf = np.zeros(0)
    else:
        wa = np.zeros(L)
        wa[finite] = np.exp(la[finite] - la[finite].max())
        arrivalPmf = wa / wa.sum()

    nvec = np.arange(L + 1)
    meanQueueLength = float(np.dot(pmf, nvec))
    utilization = float(1.0 - pmf[0])
    throughput = float(np.dot(pmf[1:], pv))
    offered = float(np.dot(pmf, boff))
    lossProb = float(pmf[L] * boff[L] / offered) if offered > 0 else 0.0
    meanWaitingQueue = meanQueueLength - utilization

    return {
        'capacity': L,
        'arrivalProb': boff,
        'serviceProb': pv,
        'pmf': pmf,
        'arrivalPmf': arrivalPmf,
        'emptyProb': float(pmf[0]),
        'utilization': utilization,
        'throughput': throughput,
        'lossProb': lossProb,
        'meanQueueLength': meanQueueLength,
        'meanWaitingQueue': meanWaitingQueue,
        'meanSojournTime': meanQueueLength / throughput if throughput > 0 else 0.0,
        'meanWaitingTime': meanWaitingQueue / throughput if throughput > 0 else 0.0,
        'normConst': H * math.exp(float(lw.max())),
        'analyzer': 'dqsys_bernoulli1',
    }


def _log(x: float) -> float:
    return math.log(x) if x > 0 else -math.inf


def _bernoulli1_infinite(b: float, p: float) -> Dict[str, object]:
    """Corollary 2.7 in closed form, with the arrival law of theorem 2.11."""
    if b >= p:
        raise ValueError('load b/p must be strictly less than 1 on an unbounded buffer')
    g = dqsys_geogeo1(b, p, LAS_DA)
    r = g['ratio']

    def arrival_pmf(n):
        arr = np.atleast_1d(np.asarray(n))
        if np.any(arr < 0) or np.any(arr != np.floor(arr)):
            raise ValueError('queue length must be a non-negative integer')
        out = (1.0 - r) * r ** arr
        return float(out[0]) if np.isscalar(n) or arr.size == 1 else out

    return {
        'capacity': math.inf,
        'arrivalProb': b,
        'serviceProb': p,
        'pmf': g['pmf'],
        'arrivalPmf': arrival_pmf,
        'emptyProb': g['emptyProb'],
        'utilization': g['utilization'],
        'throughput': g['throughput'],
        'lossProb': 0.0,
        'meanQueueLength': g['meanQueueLength'],
        'meanWaitingQueue': g['meanWaitingQueue'],
        'meanSojournTime': g['meanSojournTime'],
        'meanWaitingTime': g['meanWaitingTime'],
        'normConst': 1.0 / g['emptyProb'],
        'analyzer': 'dqsys_bernoulli1',
    }


def _expand_state(x: np.ndarray, length: int, what: str) -> np.ndarray:
    """Expand a scalar or a vector indexed by n = 0..length-1 to full length."""
    if x.size == 1:
        return np.full(length, float(x[0]))
    if x.size == length:
        return x.astype(float).ravel()
    raise ValueError('the %s probability vector must have %d entries, one per state 0..%d'
                     % (what, length, length - 1))


def _expand_service(x: np.ndarray, L: int) -> np.ndarray:
    """Expand the service probabilities to p(1..L). A vector of length L+1 is
    accepted with its first entry, which would be p(0), ignored."""
    if x.size == 1:
        return np.full(L, float(x[0]))
    if x.size == L:
        return x.astype(float).ravel()
    if x.size == L + 1:
        return x.astype(float).ravel()[1:]
    raise ValueError('the service probability vector must have %d or %d entries' % (L, L + 1))
