"""Wiring of the p-norm smoothing in the FLD `matrix` method.

Ruuskanen et al., PEVA 151 (2021) smooth the PS `min()` of eq. (12) with the
inverse p-norm of eq. (26). Two things follow that the code has to respect and
that a mean-only comparison hides:

1. An INF station has k = infinity, so `min(k, sum x) = sum x` identically and
   there is nothing to smooth. Sa carries the total population at a delay
   station, so smoothing it would run it as a k = N queue: at p = 1 that scales
   its departure rate by 1/(1 + sum x / N).
2. Eq. (23) reads the utilization off the SAME share the drift integrated,
   `k rho / E[sum X] = ghat`. Recomputing U and T from the hard min() after a
   smoothed solve breaks flow balance.

The M/M/1 case below is the paper's own Example 1: at p = 1 the smoothed model
IS the Tipper/PSFFA model `xdot = -x/(x+1) + lambda`, whose fixed point is the
exact M/M/1 mean rho/(1-rho), where the unsmoothed mean-field model returns
lambda and is arbitrarily bad as rho -> 1.
"""

import numpy as np

from line_solver import (ClosedClass, Delay, Exp, Network, OpenClass, Queue,
                         SchedStrategy, Sink, SolverFLD, Source)


def _cqn():
    """Delay(rate 1) + PS Queue(rate 2), one closed class of 10."""
    model = Network('cqn')
    d = Delay(model, 'D')
    q = Queue(model, 'Q', SchedStrategy.PS)
    c = ClosedClass(model, 'C', 10, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    model.link(Network.serialRouting(d, q))
    return model


def _mm1(rho):
    model = Network('mm1')
    s = Source(model, 'S')
    q = Queue(model, 'Q', SchedStrategy.PS)
    k = Sink(model, 'K')
    c = OpenClass(model, 'C')
    s.setArrival(c, Exp(rho))
    q.setService(c, Exp(1.0))
    model.link(Network.serialRouting(s, q, k))
    return model


def _avg(model, pstar=None):
    opt = SolverFLD.defaultOptions()
    opt.method = 'matrix'
    if pstar is not None:
        opt.pstar = pstar
    res = SolverFLD(model, opt).getAvg()
    return tuple(np.asarray(r) for r in res[:4])


def test_pnorm_leaves_the_delay_station_alone():
    # every job at a delay is in service, so its departure rate is Q * mu with
    # no share to apply. Smoothing it as a k = N queue would return Q * ghat
    # instead, which at pstar = 1 with Q near N is a factor of about 1/2.
    # (Q itself does move with the exponent: the PS queue IS smoothed, and the
    # population it holds comes off the delay.)
    for pstar in [1.0, 4.0, 20.0]:
        Q, _, _, T = _avg(_cqn(), pstar)
        assert abs(T[0, 0] - Q[0, 0]) < 1e-9, f"delay smoothed at pstar={pstar}"
    # a large exponent recovers the hard min(), i.e. the unsmoothed fixed point
    Q20, _, _, _ = _avg(_cqn(), 20.0)
    Q0, _, _, _ = _avg(_cqn())
    assert np.allclose(Q20, Q0, atol=1e-6)


def test_pnorm_metrics_use_the_smoothed_share():
    # T is read off theta, so a metric computed from the hard min() while x came
    # from the smoothed drift shows up as a flow imbalance around the cycle
    for pstar in [1.0, 4.0]:
        _, _, _, T = _avg(_cqn(), pstar)
        assert abs(T[0, 0] - T[1, 0]) < 1e-6, f"flow imbalance at pstar={pstar}: {T.flatten()}"


def test_pnorm_recovers_the_exact_mm1_mean():
    # paper Example 1: p = 1 is the Tipper/PSFFA model, exact for M/M/1
    for rho in [0.3, 0.5, 0.7]:
        Q, U, _, _ = _avg(_mm1(rho), 1.0)
        assert abs(Q[1, 0] - rho / (1.0 - rho)) < 1e-6
        # ODE tol, not a min() artifact: without the fix U would be 1.0
        assert abs(U[1, 0] - rho) < 1e-6
        # the unsmoothed mean-field model returns lambda instead, the failure
        # the smoothing exists to repair
        Q0, _, _, _ = _avg(_mm1(rho))
        assert abs(Q0[1, 0] - rho) < 1e-6
