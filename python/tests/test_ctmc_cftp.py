"""
Native-Python regression test: the perfect-sampling method of SolverCTMC,
``SolverCTMC(model, method='cftp')``, which draws iid states from the exact
stationary distribution by monotone Coupling From The Past instead of
enumerating the state space.

The reference is SolverCTMC with its default enumerating method on the same
model. Assertions cover

  - exact invariants that hold at any sample size (population conservation,
    flow balance across stations, Little's law, C = N/X, U in [0,1]);
  - accuracy against the exact CTMC, bounded by the Monte Carlo error of the
    configured sample size rather than by the method, since the sampler is
    unbiased;
  - the model-class gate, which must refuse every model outside the closed
    single-class product form rather than approximate it.

The MATLAB twin is line-test.git/test/testsCTMC/test_ctmc_cftp.m and the JAR twin is
CtmcCftpTest.java; the three assert the same invariants and bounds.
"""

import os
import sys

_WORKTREE_PY = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if sys.path and sys.path[0] != _WORKTREE_PY:
    sys.path.insert(0, _WORKTREE_PY)

import numpy as np
import pytest

from line_solver import (
    Network, Delay, Queue, Source, Sink, Exp, Erlang, ClosedClass, OpenClass,
    SchedStrategy, SolverCTMC,
)

import line_solver as _ls
assert os.path.abspath(_ls.__file__).startswith(_WORKTREE_PY), (
    f'wrong line_solver imported: {_ls.__file__} (expected under {_WORKTREE_PY})')

SAMPLES = 20000


def _cqn_model(sched_q2=SchedStrategy.FCFS, nservers=2):
    """Delay -> PS queue -> multiserver queue, one closed class of 8 jobs."""
    model = Network('CQN')
    delay = Delay(model, 'Think')
    q1 = Queue(model, 'Q1', SchedStrategy.PS)
    q2 = Queue(model, 'Q2', sched_q2)
    q2.setNumberOfServers(nservers)
    cl = ClosedClass(model, 'C1', 8, delay, 0)
    delay.setService(cl, Exp(1 / 0.5))
    q1.setService(cl, Exp(1 / 1.0))
    q2.setService(cl, Exp(1 / 0.6))
    model.link(Network.serialRouting(delay, q1, q2))
    return model


def test_invariants():
    solver = SolverCTMC(_cqn_model(), method='cftp', samples=SAMPLES, seed=7)
    QN, UN, RN, TN = solver.getAvg()[:4]
    XN = np.ravel(solver.getAvgSysTput())[0]
    CN = np.ravel(solver.getAvgSysRespT())[0]
    QN, UN, RN, TN = (np.ravel(m) for m in (QN, UN, RN, TN))
    N = 8

    # every sampled state holds the whole closed population
    assert np.isclose(QN.sum(), N, atol=1e-9)
    # a single-class series network has V = 1 everywhere, so throughput is common
    assert np.allclose(TN, TN[0], rtol=1e-12)
    # Little's law per station, and at system level
    assert np.allclose(RN, QN / TN, rtol=1e-12)
    assert np.isclose(CN, N / XN, rtol=1e-12)
    # utilization is a busy-server fraction, so it cannot leave [0,1]
    assert UN.min() >= 0.0
    assert UN[1:].max() <= 1.0


def test_accuracy_vs_exact():
    Qe, Ue = (np.ravel(m) for m in SolverCTMC(_cqn_model()).getAvg()[:2])
    Xe = np.ravel(SolverCTMC(_cqn_model()).getAvgSysTput())[0]

    solver = SolverCTMC(_cqn_model(), method='cftp', samples=SAMPLES, seed=7)
    Qc, Uc = (np.ravel(m) for m in solver.getAvg()[:2])
    Xc = np.ravel(solver.getAvgSysTput())[0]

    # 20000 iid draws put the standard error of a queue length below 1e-2 of the
    # population; the bounds below were measured at implementation time.
    assert np.abs(Qc - Qe).max() < 0.1
    assert np.abs(Uc - Ue).max() < 0.05
    assert abs(Xc - Xe) / Xe < 0.05


def test_approx_sampler_agrees():
    """The rapidly-mixing sampler targets the same distribution as the exact one."""
    model = _cqn_model(SchedStrategy.PS, 1)
    Qe = np.ravel(SolverCTMC(model).getAvgQLen())
    Qa = np.ravel(SolverCTMC(_cqn_model(SchedStrategy.PS, 1), method='cftp.approx',
                             samples=SAMPLES, seed=7).getAvgQLen())
    assert np.abs(Qa - Qe).max() < 0.15
    assert np.isclose(Qa.sum(), 8, atol=1e-9)


def test_multiserver_and_delay():
    """The balance function must carry min(n,c) at the multiserver station and
    n! at the delay: a wrong server count shows up as a biased queue length."""
    Qe = np.ravel(SolverCTMC(_cqn_model(SchedStrategy.FCFS, 3)).getAvgQLen())
    Qc = np.ravel(SolverCTMC(_cqn_model(SchedStrategy.FCFS, 3), method='cftp',
                             samples=SAMPLES, seed=11).getAvgQLen())
    assert np.abs(Qc - Qe).max() < 0.1


def test_rejects_open_model():
    model = Network('OQN')
    source = Source(model, 'Source')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    cl = OpenClass(model, 'C1')
    source.setArrival(cl, Exp(1))
    q1.setService(cl, Exp(2))
    model.link(Network.serialRouting(source, q1, sink))
    with pytest.raises((ValueError, RuntimeError), match='cftp method supports'):
        SolverCTMC(model, method='cftp', samples=100).getAvgQLen()


def test_rejects_multiclass():
    model = Network('CQN2')
    delay = Delay(model, 'Think')
    q1 = Queue(model, 'Q1', SchedStrategy.PS)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    c1 = ClosedClass(model, 'C1', 4, delay, 0)
    c2 = ClosedClass(model, 'C2', 3, delay, 0)
    delay.setService(c1, Exp(1))
    delay.setService(c2, Exp(1))
    q1.setService(c1, Exp(2))
    q1.setService(c2, Exp(3))
    q2.setService(c1, Exp(2))
    q2.setService(c2, Exp(3))
    P = model.initRoutingMatrix()
    P[c1] = Network.serialRouting(delay, q1, q2)
    P[c2] = Network.serialRouting(delay, q1, q2)
    model.link(P)
    with pytest.raises((ValueError, RuntimeError), match='single-class models only'):
        SolverCTMC(model, method='cftp', samples=100).getAvgQLen()


def test_rejects_nonexponential():
    model = Network('CQNph')
    delay = Delay(model, 'Think')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
    cl = ClosedClass(model, 'C1', 4, delay, 0)
    delay.setService(cl, Exp(2))
    q1.setService(cl, Erlang.fitMeanAndOrder(1, 2))
    q2.setService(cl, Exp(1))
    model.link(Network.serialRouting(delay, q1, q2))
    with pytest.raises((ValueError, RuntimeError), match='exponential service times'):
        SolverCTMC(model, method='cftp', samples=100).getAvgQLen()
