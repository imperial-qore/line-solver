"""
Native-Python regression test: SolverCTMC accepts and solves finite capacity
region (FCR) models under both the DROP and WAITQ rules.

Before adding 'Region' to ``SolverCTMC.getFeatureSet()``, the featset gate
rejected every FCR model at solve time (RuntimeError "... feature: Region ...")
even though the CTMC handler solves both rules exactly (DROP filters the state
space; WAITQ augments it with a per-region FIFO via ``_build_fcr_waitq_ssg``).
The MATLAB and JAR CTMC solvers already accept region models, so this closes a
cross-codebase parity gap. If the gate ever rejects an FCR model again, the
solve below raises and this test fails.

The exact means are the ground truth cross-verified across the MATLAB, JAR, and
Python-native CTMC handlers (a deterministic exact solve, so the tolerance is
tight):

    DROP, open Source(1.0)->Q1/PS(3.0)->Q2/PS(1.5)->Sink, region {Q2} cap 1:
        Q1 QLen = 0.5 (analytic: M/M/1-PS at rho=1/3 with 40% boundary loss),
        T_Q2 = 0.6.
    WAITQ, closed Think(INF,1.0)->Q1/FCFS(1.5)->Q2/FCFS(1.2), N=4, {Q2} cap 2:
        QLen = [0.966672, 1.221896, 1.379727], Tput = 0.966672 at every station.
"""

import os
import sys

_WORKTREE_PY = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if sys.path and sys.path[0] != _WORKTREE_PY:
    sys.path.insert(0, _WORKTREE_PY)

import numpy as np
import pytest

from line_solver import (
    Network, Delay, Queue, Source, Sink, Exp, ClosedClass, OpenClass,
    SchedStrategy, DropStrategy, SolverCTMC,
)

import line_solver as _ls
assert os.path.abspath(_ls.__file__).startswith(_WORKTREE_PY), (
    f'wrong line_solver imported: {_ls.__file__} (expected under {_WORKTREE_PY})')

RTOL = 1e-3


def _drop_model():
    """Open DROP FCR model: Source -> Q1/PS -> Q2/PS{region cap 1} -> Sink."""
    model = Network('ctmc_fcr_drop')
    source = Source(model, 'Source')
    q1 = Queue(model, 'Q1', SchedStrategy.PS)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    sink = Sink(model, 'Sink')
    c1 = OpenClass(model, 'C1')
    source.setArrival(c1, Exp(1.0))
    q1.setService(c1, Exp(3.0))
    q2.setService(c1, Exp(1.5))
    model.link(Network.serialRouting(source, q1, q2, sink))
    region = model.addRegion([q2])
    region.setGlobalMaxJobs(1)
    region.setDropRule(c1, DropStrategy.DROP)
    return model


def _waitq_model():
    """Closed WAITQ FCR model: Think -> Q1/FCFS -> Q2/FCFS{region cap 2} -> Think."""
    model = Network('ctmc_fcr_waitq')
    think = Delay(model, 'Think')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
    c1 = ClosedClass(model, 'C1', 4, think)
    think.setService(c1, Exp(1.0))
    q1.setService(c1, Exp(1.5))
    q2.setService(c1, Exp(1.2))
    model.link(Network.serialRouting(think, q1, q2))
    region = model.addRegion([q2])
    region.setGlobalMaxJobs(2)
    region.setDropRule(c1, DropStrategy.WaitingQueue)
    return model


def test_ctmc_accepts_and_solves_drop_region():
    """DROP FCR: previously rejected at the featset gate; now solves exactly."""
    sv = SolverCTMC(_drop_model(), cutoff=16, keep=False)
    q = np.atleast_2d(sv.getAvgQLen())
    t = np.atleast_2d(sv.getAvgTput())
    # Station rows: 0 = Source, 1 = Q1, 2 = Q2.
    assert abs(q[1, 0] - 0.5) < RTOL, f'DROP Q1 QLen: CTMC {q[1, 0]:.6f}, expected 0.5'
    assert abs(t[2, 0] - 0.6) < RTOL, f'DROP T_Q2: CTMC {t[2, 0]:.6f}, expected 0.6'


def test_ctmc_accepts_and_solves_waitq_region():
    """WAITQ FCR: previously rejected at the featset gate; now solves exactly,
    parking refused jobs in the region FIFO so the station queue lengths sum to
    below the closed population."""
    sv = SolverCTMC(_waitq_model(), cutoff=8, keep=False)
    q = np.atleast_2d(sv.getAvgQLen())
    t = np.atleast_2d(sv.getAvgTput())
    exact_q = np.array([0.966672, 1.221896, 1.379727])
    for i in range(3):
        assert abs(q[i, 0] - exact_q[i]) < RTOL, (
            f'WAITQ station {i} QLen: CTMC {q[i, 0]:.6f}, expected {exact_q[i]:.6f}')
        assert abs(t[i, 0] - 0.966672) < RTOL, (
            f'WAITQ station {i} Tput: CTMC {t[i, 0]:.6f}, expected 0.966672')
    # The parked mass lives in the region FIFO, off any station.
    assert float(np.sum(q[:3, 0])) < 3.9, (
        f'WAITQ station QLen sum {float(np.sum(q[:3, 0])):.4f} should be below N=4')
