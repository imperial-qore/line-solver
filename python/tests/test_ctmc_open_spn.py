"""
Native-Python regression test: SolverCTMC solves BOUNDED open stochastic Petri
nets and rejects only the genuinely-unbounded case.

An open SPN feeds tokens from a Source into a Place, drained by a Transition to
a Sink. Previously SolverCTMC raised a blanket RuntimeError for every open SPN
because the CTMC state machine never emitted the Source->Place arrival: the ARV
event misassigned the arriving token to a server-phase slot the ordinary place
lacks, so the marking never grew and the generator degenerated to a pure-death
chain (reported mean == cutoff/2, the uniform-fabrication signature). With the
arrival handled at the marking level (State.after_event_station NodeType.PLACE
branch) and a firing consuming exactly the arc weight (after_global_event
_apply_pre_post), the model

    Source Exp(0.5) -> Place P1 -> Transition T1 Exp(1.0) -> Sink

with the marking bounded by a finite Place capacity K (setClassCapacity) or by
the solver cutoff, is an M/M/1/K queue (rho=0.5). The exact stationary mean is

    E[P1] = sum_{n=0..K} n*rho^n / sum_{n=0..K} rho^n,

and the throughput equals the effective arrival rate lambda*(1 - P_block).
Both are matched to solver tolerance. A net whose Place has infinite capacity
and no finite cutoff is genuinely unbounded and must still error clearly.
"""

import os
import sys

_WORKTREE_PY = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if sys.path and sys.path[0] != _WORKTREE_PY:
    sys.path.insert(0, _WORKTREE_PY)

import numpy as np
import pytest

from line_solver import (
    Network, Source, Sink, Place, Transition, OpenClass, Exp, CTMC,
)

import line_solver as _ls
assert os.path.abspath(_ls.__file__).startswith(_WORKTREE_PY), (
    f'wrong line_solver imported: {_ls.__file__} (expected under {_WORKTREE_PY})')

RTOL = 1e-6
LAMBDA = 0.5
MU = 1.0
RHO = LAMBDA / MU


def _mm1k_mean(K):
    den = sum(RHO ** n for n in range(K + 1))
    num = sum(n * RHO ** n for n in range(K + 1))
    return num / den


def _mm1k_tput(K):
    den = sum(RHO ** n for n in range(K + 1))
    p_block = RHO ** K / den
    return LAMBDA * (1.0 - p_block)


def _build(cap=None):
    model = Network('ospn')
    source = Source(model, 'source')
    sink = Sink(model, 'sink')
    p1 = Place(model, 'P1')
    t1 = Transition(model, 'T1')
    jc = OpenClass(model, 'jobs')
    source.setArrival(jc, Exp.fitMean(1.0 / LAMBDA))
    if cap is not None:
        p1.setClassCapacity(jc, cap)
    m = t1.addMode('fire')
    t1.setDistribution(m, Exp.fitMean(1.0 / MU))
    # Enabling arc consumes one token from P1 on firing; the token is routed to
    # the Sink (no firing outcome on P1, which would double-consume).
    t1.setEnablingConditions(m, jc, p1, 1)
    R = model.initRoutingMatrix()
    R.set(jc, jc, source, p1, 1.0)
    R.set(jc, jc, p1, t1, 1.0)
    R.set(jc, jc, t1, sink, 1.0)
    model.link(R)
    return model


def _solve(model, cutoff):
    solver = CTMC(model)
    solver.options.cutoff = cutoff
    r = solver.getAvg()
    QN, TN = r[0], r[3]
    pidx = model.getStationNames().index('P1')
    return float(np.sum(QN[pidx, :])), float(np.sum(TN[pidx, :]))


@pytest.mark.parametrize('K', [8, 12])
def test_capacity_bounded_matches_mm1k(K):
    """Finite-capacity Place bounds the marking -> exact M/M/1/K."""
    qP1, xP1 = _solve(_build(cap=K), cutoff=K + 4)
    assert qP1 == pytest.approx(_mm1k_mean(K), rel=RTOL)
    assert xP1 == pytest.approx(_mm1k_tput(K), rel=RTOL)


@pytest.mark.parametrize('K', [8, 12])
def test_cutoff_bounded_matches_mm1k(K):
    """Solver cutoff bounds the marking -> exact M/M/1/K."""
    qP1, xP1 = _solve(_build(cap=None), cutoff=K)
    assert qP1 == pytest.approx(_mm1k_mean(K), rel=RTOL)
    assert xP1 == pytest.approx(_mm1k_tput(K), rel=RTOL)


def test_unbounded_open_spn_errors():
    """Infinite-capacity Place with a non-finite cutoff cannot build a finite
    generator and must be rejected with a clear message."""
    model = _build(cap=None)
    solver = CTMC(model)
    solver.options.cutoff = float('inf')
    with pytest.raises(RuntimeError, match='unbounded'):
        solver.getAvgTable()
