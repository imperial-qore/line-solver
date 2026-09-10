"""
Service rate of a station served by heterogeneous server pools with a
class-compatibility graph, and the peak that normalizes its utilization.

Port of ``matlab/src/api/sn/sn_compat_rate.m``.
"""

import numpy as np

__all__ = ['sn_compat_rate', 'sn_compat_peak', 'sn_compat_scaling']


def sn_compat_rate(compat, counts, rates, n):
    """
    Total service rate of a compatibility-structured station.

    A pool t holds ``counts[t]`` identical servers, each running at
    ``rates[t]``, and may serve operand j when ``compat[t, j]`` is nonzero. The
    rate the station clears in state ``n`` is::

        mu(n) = sum_t counts[t]*rates[t]*min(1, sum_{j: compat[t,j] != 0} n[j])

    the ACTIVATED-SERVER law: a pool contributes its full rate as soon as it is
    compatible with at least one operand PRESENT. This is the order-independent
    reading of a compatibility structure -- at an INTEGER state mu depends on n
    only through its SUPPORT, so it is invariant to the arrival order and to any
    permutation of the microstate, which is exactly the condition an OI station
    has to meet (Dorsman & Gardner, Queueing Systems 107:205-256, 2024, Fig. 1).
    It is also what ``pas_compatibility_5class.m`` encodes for a flat Network, so
    the layered and flat readings of one compatibility matrix agree.

    WHY min(1, .) AND NOT AN INDICATOR. At every integer state the two agree
    exactly -- a pool with at least one compatible job present is fully active,
    one with none is idle -- so nothing about the OI law on the real state
    lattice changes. They part company only at a FRACTIONAL argument, which is
    what a mean-value solver hands this function: AMVA evaluates the rate at a
    MEAN population, and under a hard indicator any operand with a mean above
    zero, however small, activates every pool it touches. A compatibility
    structure would then be invisible to AMVA whenever every operand is a little
    bit busy -- which is nearly always. Scaling linearly below one job keeps the
    structure visible at the evaluation point while leaving the integer-state law
    untouched; it is the ordinary continuous relaxation of a step function, and
    the CTMC and simulation paths, which only ever evaluate at integer states,
    cannot tell the difference.

    IT IS NOT A MATCHING. A pool of two servers compatible with a class holding
    ONE job contributes both servers here, which over-counts against a
    non-redundant system where one server serves one job. That is deliberate:
    the matching size depends on the counts and not only on the support, so it
    is NOT order independent and would take the station outside the product form
    the OI closure is built on. A model that means the matching wants a
    different station, not a different reading of this one.

    :param compat: (npools x noperands) array, nonzero where the pool may serve
    :param counts: (npools) servers held by each pool
    :param rates: (npools) per-server rate of each pool
    :param n: (noperands) per-operand population, integer or fractional
    :return: the total service rate mu(n)
    """
    compat = np.atleast_2d(np.asarray(compat, dtype=float))
    counts = np.atleast_1d(np.asarray(counts, dtype=float)).ravel()
    rates = np.atleast_1d(np.asarray(rates, dtype=float)).ravel()
    n = np.atleast_1d(np.asarray(n, dtype=float)).ravel()
    npools = counts.size
    if rates.size != npools:
        raise ValueError('sn_compat_rate: one rate per pool is required')
    if compat.shape[0] != npools:
        raise ValueError('sn_compat_rate: compat must have one row per pool')
    if n.size != compat.shape[1]:
        raise ValueError('sn_compat_rate: n must have one entry per operand')
    # The pool is activated ONCE by the jobs it can reach, not once per operand:
    # its weight is the compatible load, capped at one job.
    load = np.minimum(1.0, (compat != 0) @ np.maximum(n, 0.0))
    return float(np.sum(counts * rates * load))


def sn_compat_peak(counts, rates):
    """
    Rate a compatibility declaration clears with every pool active,
    ``sum_t counts[t]*rates[t]``.

    Utilization at a rate-scaled station is reported as U = T*S/peak, and the
    peak is a property of the DECLARATION rather than of a state, so it is
    computed once and handed to the solver beside the rate handle rather than
    recovered from :func:`sn_compat_rate` at a guessed state.
    """
    counts = np.atleast_1d(np.asarray(counts, dtype=float)).ravel()
    rates = np.atleast_1d(np.asarray(rates, dtype=float)).ravel()
    if rates.size != counts.size:
        raise ValueError('sn_compat_peak: one rate per pool is required')
    return float(np.sum(counts * rates))


def sn_compat_scaling(compat, counts, rates, n):
    """
    Rate scaling eta(n) a compatibility declaration imposes on its station.

    This is what SolverLN carries onto the layer station, and it is NOT
    ``sn_compat_rate / sn_compat_peak``. The denominator is the rate the SAME
    population would obtain under FULL compatibility::

        eta(n) = mu(n) / (peak * min(1, sum_j n_j / S))

    so eta isolates the effect of the compatibility GRAPH and nothing else. The
    denominator DAMPS BY OCCUPANCY RELATIVE TO THE SERVER COUNT, min(1, N/S),
    because that is precisely what the solver's own multiserver term
    contributes: it applies min(N,S) servers at the average server rate peak/S,
    so::

        min(N,S) * (peak/S) * eta(n) = mu(n)

    and the station clears the activated-server rate exactly, at every state.

    DAMPING BY min(1, N) INSTEAD -- which this did until 2026-08-28 -- leaves the
    effective law at min(N,S)/S * mu(n), which cancels the REDUNDANCY SPEED-UP
    the activated-server law exists to express: a pool of S servers facing one
    compatible job clears S, not 1, because every one of them works on it and
    the first to finish cancels the rest. Under the old normalization a fully
    compatible pool reduced to the plain multiserver, so the OI machinery did no
    work in the homogeneous case and LDES, which simulates mu(n) directly,
    disagreed with it by that factor.

    eta is therefore ABOVE ONE at low occupancy, which is not a defect: it is the
    speed-up carried by servers that would otherwise be idle. A FULLY-COMPATIBLE
    POOL IS THEREFORE NOT THE NEUTRAL eta == 1 -- it is min(1, N) / min(1, N/S),
    which is S below one job, S/N between one job and S, and 1 from S jobs up.
    `tests/test_lqn_server_pools.py` pins that law both directly and through the
    layer station.
    """
    counts = np.atleast_1d(np.asarray(counts, dtype=float)).ravel()
    rates = np.atleast_1d(np.asarray(rates, dtype=float)).ravel()
    n = np.atleast_1d(np.asarray(n, dtype=float)).ravel()
    total = float(np.sum(np.maximum(n, 0.0)))
    if total <= 0:
        return 1.0                      # an empty station: nothing to scale
    nservers = float(np.sum(counts))
    if nservers <= 0:
        return 1.0
    ref = sn_compat_peak(counts, rates) * min(1.0, total / nservers)
    return sn_compat_rate(compat, counts, rates, n) / ref
