"""One load-dependent fixture, every solver, ONE utilization convention.

Utilization has two readings that coincide for a load-independent single server
and part company the moment alpha(n) != 1:

  work-based   U = X * E[S] / max(c, max alpha)   -- the fraction of the
                                                     station's PEAK capacity
                                                     actually being delivered
  time-based   U = P(at least one server busy)    -- the fraction of time the
                                                     station is occupied

A server running alpha(n) times faster does the same work in less time, so the
time-based reading calls it no busier than one at its nominal rate. On the
fixture below that is 0.9587 against 0.6612 -- a 45% spread on the same model.

LINE reports the WORK-BASED number everywhere. It used to be split: MAM's
LD-QBD, the NRM SSA engine and the C++ LDES engine measured busy time while
CTMC, MVA, NC, serial SSA and the Java LDES engine measured work, so SolverSSA
disagreed with itself depending on which engine ran. This file is the guard on
that: the exact solvers are pinned to a closed-form number and every other
solver to them.

Author: QORE Lab, Imperial College London
"""
import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, ClosedClass, Exp, SchedStrategy,
                         SolverCTMC, SolverMVA, SolverNC, SolverMAM, SolverSSA,
                         SolverFluid, circul)

QI = 1                       # 0 = Delay, 1 = Queue
ALPHA = [1.0, 1.5, 2.0, 2.5]
N = 4

# Closed form for the fixture. The queue is a birth-death chain on n = 0..4 with
# birth (N-n)*1.0 and death alpha(n)*1.0, so X = 1.652892561983471 and the
# work-based utilization is X * E[S] / max(alpha) = X / 2.5.
X_EXACT = 1.6528925619834711
U_EXACT = X_EXACT / max(ALPHA)          # 0.6611570247933884
Q_EXACT = 2.3471074380165293


def _model():
    model = Network('lld_util')
    delay = Delay(model, 'Delay')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    jobclass = ClosedClass(model, 'Class1', N, delay, 0)
    delay.setService(jobclass, Exp(1.0))
    queue.setService(jobclass, Exp(1.0))
    queue.setLoadDependence(np.asarray(ALPHA, dtype=float))
    P = model.initRoutingMatrix()
    P[0] = circul(2)
    model.link(P)
    return model


def _queue_row(solver):
    avg = solver.getAvg()
    return (float(np.asarray(avg[0])[QI, 0]),      # QLen
            float(np.asarray(avg[1])[QI, 0]),      # Util
            float(np.asarray(avg[3])[QI, 0]))      # Tput


def test_the_closed_form_is_what_ctmc_computes():
    # anchors the other assertions to arithmetic rather than to a solver
    q, u, x = _queue_row(SolverCTMC(_model()))
    assert q == pytest.approx(Q_EXACT, abs=1e-9)
    assert u == pytest.approx(U_EXACT, abs=1e-9)
    assert x == pytest.approx(X_EXACT, abs=1e-9)
    # and it is emphatically NOT the time-based reading
    assert u == pytest.approx(0.6611570248, abs=1e-9)
    assert abs(u - 0.9586776860) > 0.29


@pytest.mark.parametrize("name,make", [
    ('MVA', lambda m: SolverMVA(m)),
    ('MVA exact', lambda m: SolverMVA(m, method='exact')),
    ('NC', lambda m: SolverNC(m)),
    ('MAM ldqbd', lambda m: SolverMAM(m, method='ldqbd')),
])
def test_exact_solvers_agree_to_machine_precision(name, make):
    # MAM's LD-QBD is an exact chain like CTMC's, so the bar is not a percentage
    q, u, x = _queue_row(make(_model()))
    assert q == pytest.approx(Q_EXACT, abs=1e-9), name
    assert u == pytest.approx(U_EXACT, abs=1e-9), name
    assert x == pytest.approx(X_EXACT, abs=1e-9), name


@pytest.mark.parametrize("method", ['nrm', 'serial'])
def test_both_ssa_engines_report_the_same_convention(method):
    # the two engines measure it differently -- NRM integrates busy time and
    # serial applies the utilization law -- and used to differ by 0.30 here
    q, u, x = _queue_row(SolverSSA(_model(), method=method, samples=20000, seed=23000))
    assert u == pytest.approx(U_EXACT, abs=2e-2), method
    assert q == pytest.approx(Q_EXACT, abs=5e-2), method
    assert x == pytest.approx(X_EXACT, abs=5e-2), method


def test_fluid_is_on_the_same_convention():
    # an approximation, so the bar is loose; what matters is which number it is
    # approximating -- 0.66, not 0.96
    _, u, _ = _queue_row(SolverFluid(_model()))
    assert u == pytest.approx(U_EXACT, abs=2e-2)


def test_utilization_stays_below_one_at_saturation():
    # the time-based reading saturates at 1 while the station still has capacity
    # left; the work-based one reaches 1 only when alpha is at its peak
    model = Network('lld_sat')
    delay = Delay(model, 'Delay')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    jobclass = ClosedClass(model, 'Class1', 6, delay, 0)
    delay.setService(jobclass, Exp(100.0))          # a near-zero think time
    queue.setService(jobclass, Exp(1.0))
    queue.setLoadDependence(np.asarray([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], dtype=float))
    P = model.initRoutingMatrix()
    P[0] = circul(2)
    model.link(P)
    _, u_ctmc, _ = _queue_row(SolverCTMC(model))
    _, u_mam, _ = _queue_row(SolverMAM(model, method='ldqbd'))
    assert u_mam == pytest.approx(u_ctmc, abs=1e-9)
    assert 0.0 < u_mam < 1.0
