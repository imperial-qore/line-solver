"""
The arrival rate of a station-class pair as a function of time.

Native Python twin of matlab/src/api/sn/sn_arrival_rate_fun.m. The time-varying
analyses (Mt/G/inf, the modified offered load, the Gt/Mt/st+GI fluid queue)
consume lambda(t) itself, not a mean rate: their whole content is the LAG
between when work arrives and when it is felt, and a time-averaged rate has no
lag.
"""

from typing import Any, Callable, Tuple

import numpy as np


def sn_arrival_rate_fun(sn, ist: int, r: int) -> Tuple[Callable[[Any], Any], bool, float]:
    """
    Build lambda(t) for station ``ist``, class ``r``.

    LINE carries a time-varying arrival as a MAPt or an NHPP, whose ``sn.proc``
    slot holds a piecewise-constant schedule, so lambda(t) is read off the
    segment in force at t. For any other process the rate is constant and the
    handle returns it, which is what lets a caller ask for the time-varying
    analysis of a stationary model and get the stationary answer rather than an
    error.

    Args:
        sn: the NetworkStruct
        ist: station index
        r: class index

    Returns:
        ``(lambdaFun, isTimeVarying, period)``; ``period`` is the cycle length
        when the schedule is cyclic and ``inf`` otherwise.

    See also:
        matlab/src/api/sn/sn_arrival_rate_fun.m
    """
    from ...solvers.solver_fld.utils.phase_type import (is_mapt, is_pht, is_rate_schedule,
                                                        schedule_segments,
                                                        schedule_map_pair_at)
    rate = float(np.asarray(sn.rates)[ist, r])
    def _shape(t, values):
        # Scalar in, scalar out, as every consumer of a rate function assumes:
        # the fluid integrators call lambda(t) at one instant and wrap the
        # result in float(), which an array of length one no longer satisfies
        # under numpy 2.
        return float(values[0]) if np.ndim(t) == 0 else values

    if not is_rate_schedule(sn, ist, r):
        def const(t, _r=rate):
            tt = np.atleast_1d(np.asarray(t, dtype=float))
            return _shape(t, np.full(tt.shape, _r))
        return const, False, float('inf')

    slot = sn.proc[ist][r]
    if is_mapt(sn, ist, r) or is_pht(sn, ist, r):
        kind = 'MAPt' if is_mapt(sn, ist, r) else 'PHt'
        bp, pairs, cyclic = schedule_segments(slot, kind)
        seg_rate = np.zeros(len(pairs))
        for k in range(len(pairs)):
            D0, D1 = schedule_map_pair_at(slot, kind, k)
            if D0.shape[0] == 1:
                seg_rate[k] = float(D1[0, 0])
            else:
                # The arrival rate of a segment is pie_k D1_k e, the stationary
                # throughput of that segment's own MAP.
                from ..mam.map_analysis import map_pie
                pie = np.asarray(map_pie(D0, D1), dtype=float).ravel()
                seg_rate[k] = float(pie @ D1 @ np.ones(D1.shape[1]))
    else:
        # An NHPP slot is {breakpoints, rates, cyclic}: the rates ARE lambda(t),
        # one per interval, so there is no MAP pair to reduce. The layout
        # differs from the MAPt one and reading it as that raises rather than
        # returning a wrong rate.
        bp = np.asarray(slot[0], dtype=float).ravel()
        seg_rate = np.asarray(slot[1], dtype=float).ravel()
        cyclic = bool(slot[2])
    time_varying = bool(np.any(np.abs(seg_rate - seg_rate[0]) > 1e-12))
    period = float(bp[-1] - bp[0]) if cyclic else float('inf')

    def lam(t, _bp=bp, _sr=seg_rate, _cyc=cyclic):
        tt = np.atleast_1d(np.asarray(t, dtype=float))
        u = tt
        if _cyc and _bp[-1] > _bp[0]:
            u = _bp[0] + np.mod(tt - _bp[0], _bp[-1] - _bp[0])
        # Segment k is in force on [bp[k], bp[k+1]). Before the first
        # breakpoint the first segment holds and after the last the last one
        # does, so a caller integrating over an infinite past (the Mt/G/inf
        # convolution) gets a defined rate everywhere rather than a NaN.
        idx = np.searchsorted(_bp[:-1], u, side='right') - 1
        idx = np.clip(idx, 0, len(_sr) - 1)
        return _shape(t, _sr[idx])

    return lam, time_varying, period
