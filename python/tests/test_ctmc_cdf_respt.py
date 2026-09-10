"""
Cross-codebase parity of the CTMC response-time distribution.

THE GOLDEN IS MATLAB, C++ AND THE JAR, WHICH AGREE TO TEN DECIMALS. On the model
below (Delay Exp(2) + FCFS Queue Exp(1), one closed class of 2 jobs) all three
report, on the same 0.999-truncated grid,

    Think : n = 3455 points, F(end) = 0.9990002447, E[T] = 0.4995002888
    Q     : n = 8840 points, F(end) = 0.9990007864, E[T] = 1.6655708356

E[T] is the trapezoidal integral of the survival function over the grid, so it
carries the truncation and is NOT the exact mean response time (which is 0.5 and
5/3). That is deliberate: it makes the four codebases comparable point for point
rather than only in the limit.

The LDES simulation of the same model at seed 23000 with 200k samples gives
0.501165 and 1.671497 in MATLAB, Python and C++ alike, i.e. within 0.36 per cent
of the exact law -- an independent confirmation by a different method.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, Network, Queue, SchedStrategy,
                         SolverCTMC)


def _model():
    model = Network('cqn2')
    think = Delay(model, 'Think')
    q = Queue(model, 'Q', SchedStrategy.FCFS)
    c = ClosedClass(model, 'C', 2, think)
    think.setService(c, Exp(2.0))
    q.setService(c, Exp(1.0))
    model.link(Network.serialRouting(think, q))
    return model


def _mean(entry):
    return float(np.trapezoid(1.0 - np.asarray(entry['p']), np.asarray(entry['t'])))


def test_getcdfrespt_returns_curves_not_an_exponential_fit():
    """
    The native path used to fit an exponential to the mean response time, which
    is exact only for an M/M/1. A tagged-chain law is not exponential, and this
    is what tells the two apart: an exponential of mean m has F(m) = 1-1/e =
    0.6321 exactly, at every station and every model.
    """
    rd = SolverCTMC(_model()).getCdfRespT()
    assert len(rd) == 2
    for e in rd:
        t = np.asarray(e['t'])
        F = np.asarray(e['p'])
        assert t.size > 100, "a distribution, not a summary"
        assert np.all(np.diff(F) >= -1e-12), "the CDF must be monotone"
        assert np.all(np.diff(t) > 0), "the grid must increase"
        m = _mean(e)
        Fm = float(np.interp(m, t, F))
        assert abs(Fm - 0.6321205588) > 1e-4, (
            "F(mean) landed on the exponential value, so this is still a fit")


def test_curves_match_the_matlab_cpp_and_jar_golden():
    rd = {e['station']: e for e in SolverCTMC(_model()).getCdfRespT()}

    think = rd[1]
    queue = rd[2]
    assert np.asarray(think['t']).size == 3455
    assert np.asarray(queue['t']).size == 8840
    assert np.asarray(think['p'])[-1] == pytest.approx(0.9990002447, abs=1e-9)
    assert np.asarray(queue['p'])[-1] == pytest.approx(0.9990007864, abs=1e-9)
    assert _mean(think) == pytest.approx(0.4995002888, abs=1e-9)
    assert _mean(queue) == pytest.approx(1.6655708356, abs=1e-9)


def test_the_law_converges_to_the_exact_mean_response_time():
    """
    The grid mean is the exact mean minus the mass past the 0.999 truncation.
    Delay: exactly 0.5. Queue at N=2 with Z=0.5 and D=1: MVA gives R = 5/3.
    """
    rd = {e['station']: e for e in SolverCTMC(_model()).getCdfRespT()}
    assert _mean(rd[1]) == pytest.approx(0.5, rel=2e-3)
    assert _mean(rd[2]) == pytest.approx(5.0 / 3.0, rel=2e-3)


def test_getcdfsysrespt_is_the_cycle_time_law_per_chain():
    """
    The system response time is the CYCLE TIME: one full trip of a tagged job
    round the network, from one arrival at its reference station to the next.
    For this model that is Z + R = 0.5 + 5/3 = 13/6.

    MATLAB, the JAR and C++ all report, on the same 10000-interval grid
    truncated at 1 - FineTol: n = 2183, F(end) = 0.9999999901,
    E[T] = 2.1666666563.
    """
    rd = SolverCTMC(_model()).getCdfSysRespT()
    assert len(rd) == 1, "one law per chain"
    e = rd[0]
    assert e['chain'] == 1
    t = np.asarray(e['t'])
    F = np.asarray(e['p'])
    assert t.size == 2183
    assert F[-1] == pytest.approx(0.9999999901, abs=1e-9)
    assert _mean(e) == pytest.approx(2.1666666563, abs=1e-9)
    assert _mean(e) == pytest.approx(13.0 / 6.0, rel=1e-6)
