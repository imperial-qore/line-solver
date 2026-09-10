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
    Network, Delay, Queue, Router, Source, Sink, Exp, ClosedClass, OpenClass,
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


def _immediate_region_model():
    """Closed FCR model carrying vanishing states: the Router is stateful but
    not a station, so every marking holding the job there is eliminated by
    stochastic complementation.

    This is the shape LQN2QN emits (a thread-pool region plus reply
    pseudo-nodes), which is what exposed the two defects guarded below.
    """
    model = Network('ctmc_fcr_immediate')
    think = Delay(model, 'Think')
    router = Router(model, 'R')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
    c1 = ClosedClass(model, 'C1', 2, think)
    think.setService(c1, Exp(1.0))
    q1.setService(c1, Exp(1.0))
    q2.setService(c1, Exp(2.0))
    model.link(Network.serialRouting(think, router, q1, q2))
    region = model.addRegion([q1, q2])
    region.setGlobalMaxJobs(2)
    return model


@pytest.mark.skipif(
    os.environ.get('LINE_SOLVER_LANG') == 'java',
    reason="intercepts handler.ctmc_solve_reducible_blkdecomp, a native-only "
           "entry point; under lang='java' runAnalyzer delegates the whole "
           "solve to jline.jar and no native solver code runs. The same "
           "invariant is asserted on the public generator by "
           "test_ctmc_fcr_generator_is_proper_from_getinfgen, which runs "
           "under both langs")
def test_ctmc_fcr_generator_is_proper_with_vanishing_states():
    """The WAITQ-FCR builder seeds its generator with an identity placeholder and
    the caller skips _build_generator_sync on that path, so the placeholder has
    to be removed there. Left in, stochastic complementation over the vanishing
    states produced a negative off-diagonal and non-zero row sums, and the solve
    returned a degenerate distribution."""
    import importlib
    handler = importlib.import_module('line_solver.api.solvers.ctmc.handler')
    seen = {}
    # the handler now solves every chain through the block decomposition, so the
    # generator is captured at that entry point rather than at ctmc_solve
    original = handler.ctmc_solve_reducible_blkdecomp

    def _capture(Q, *args, **kwargs):
        dense = np.asarray(Q.todense() if hasattr(Q, 'todense') else Q, dtype=float)
        seen.setdefault('Q', dense)
        return original(Q, *args, **kwargs)

    handler.ctmc_solve_reducible_blkdecomp = _capture
    try:
        SolverCTMC(_immediate_region_model(), cutoff=4, keep=False).getAvg()
    finally:
        handler.ctmc_solve_reducible_blkdecomp = original

    Q = seen['Q']
    assert np.allclose(Q.sum(axis=1), 0.0, atol=1e-6), (
        f'generator rows must sum to zero, got {Q.sum(axis=1)}')
    offdiag = Q - np.diag(np.diag(Q))
    assert offdiag.min() >= -1e-9, (
        f'generator off-diagonals must be non-negative, min {offdiag.min():g}')


def test_ctmc_fcr_generator_is_proper_from_getinfgen():
    """Same invariant as the test above, read off the PUBLIC generator instead
    of an intercepted native call, so it holds under every lang. Both codebases
    return the same 6-state chain here: native eliminates the three vanishing
    Router markings by stochastic complementation, the JAR never enumerates
    them, and pi must solve the balance equations of the generator it is
    returned with."""
    sv = SolverCTMC(_immediate_region_model(), cutoff=4, keep=False)
    sv.getAvg()
    Q = np.atleast_2d(np.asarray(sv.getInfGen(), dtype=float))
    assert Q.shape[0] == Q.shape[1], f'generator is not square: {Q.shape}'
    assert np.allclose(Q.sum(axis=1), 0.0, atol=1e-9), (
        f'generator rows must sum to zero, got {Q.sum(axis=1)}')
    offdiag = Q - np.diag(np.diag(Q))
    assert offdiag.min() >= -1e-9, (
        f'generator off-diagonals must be non-negative, min {offdiag.min():g}')
    pi = np.asarray(sv.getSteadyState(), dtype=float).reshape(-1)
    assert pi.size == Q.shape[0], (
        f'pi ({pi.size}) does not match its generator {Q.shape}')
    assert abs(pi.sum() - 1.0) < 1e-9, f'pi does not sum to one: {pi.sum():g}'
    assert np.max(np.abs(pi @ Q)) < 1e-9, (
        f'pi is not stationary for the returned generator: '
        f'||pi Q||_inf = {np.max(np.abs(pi @ Q)):g}')


def test_ctmc_fcr_rates_survive_vanishing_state_elimination():
    """On the FCR path the arrival/departure rates were rebuilt from an empty
    fork-join action list whenever vanishing states were present, zeroing Tput
    and RespT; the getAvg near-zero mask then zeroed QLen and Util too, so the
    average table came back empty."""
    sv = SolverCTMC(_immediate_region_model(), cutoff=4, keep=False)
    t = np.atleast_2d(sv.getAvgTput())
    q = np.atleast_2d(sv.getAvgQLen())
    assert np.max(t) > 0.0, 'every station throughput was zero'
    # a closed tandem conserves flow: all three stations carry one chain
    assert np.allclose(t[:3, 0], t[0, 0], rtol=1e-6), (
        f'flow balance violated across the tandem: {t[:3, 0]}')
    exact_q = np.array([0.666667, 0.933333, 0.400000])
    for i in range(3):
        assert abs(q[i, 0] - exact_q[i]) < RTOL, (
            f'station {i} QLen: CTMC {q[i, 0]:.6f}, expected {exact_q[i]:.6f}')
    assert abs(float(np.sum(q[:3, 0])) - 2.0) < 1e-6, (
        f'closed population not conserved: QLen sum {float(np.sum(q[:3, 0])):.6f}')
