"""
Native-Python tests for the dedicated NRM SSA solver with POLLING scheduling.

A closed two-class cyclic network (Delay + single-server POLLING Queue) is
solved with the dedicated next-reaction-method engine and checked against the
exact SolverCTMC ground truth within a loose tolerance that accommodates
simulation noise. Every polling discipline (EXHAUSTIVE, GATED, KLIMITED,
DECREMENTING) is covered with and without exponential switchover.

Each case also asserts that the dedicated NRM handler actually ran (method ==
'nrm'): a silent fallback to the serial engine would still produce a
serial-correct number and hide a broken NRM, so the assertion is mandatory.
"""

import os
import numpy as np

import line_solver
from line_solver import (
    Network, Delay, Queue, Exp, ClosedClass, SchedStrategy, PollingType,
    SolverCTMC, SolverSSA,
)

SAMPLES = 300000
REL_TOL = 0.06   # 6% relative tolerance for non-trivial metrics
ABS_TOL = 0.03   # absolute floor for near-zero metrics


def test_imports_worktree_line_solver():
    """Guard against picking up a globally installed line_solver: this test
    must exercise the worktree copy that carries the polling NRM port."""
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    mod = os.path.abspath(line_solver.__file__)
    assert mod.startswith(here), \
        'line_solver imported from %s, not the worktree under %s' % (mod, here)


def _build(ptype, switchover, par=None, n1=2, n2=2):
    model = Network('poll')
    delay = Delay(model, 'Think')
    queue = Queue(model, 'Poll', SchedStrategy.POLLING)
    c1 = ClosedClass(model, 'C1', n1, delay)
    c2 = ClosedClass(model, 'C2', n2, delay)
    delay.setService(c1, Exp(0.5))
    delay.setService(c2, Exp(0.7))
    queue.setService(c1, Exp(1.0))
    queue.setService(c2, Exp(1.3))
    if par is not None:
        queue.setPollingType(ptype, par)
    else:
        queue.setPollingType(ptype)
    if switchover:
        queue.setSwitchover(c1, Exp(2.0))
        queue.setSwitchover(c2, Exp(3.0))
    P = model.initRoutingMatrix()
    P[c1, c1] = Network.serialRouting(delay, queue)
    P[c2, c2] = Network.serialRouting(delay, queue)
    model.link(P)
    return model


def _assert_nrm_path(model):
    """Guarantee the dedicated NRM handler runs for this model (not a fallback)."""
    from line_solver.api.solvers.ssa.handler import solver_ssa, SolverSSAOptions
    model.refresh_struct()
    sn = model._sn
    opt = SolverSSAOptions()
    opt.method = 'nrm'
    opt.samples = 200
    opt.seed = 1
    ret = solver_ssa(sn, opt, model)
    assert ret.method == 'nrm', \
        'expected dedicated NRM solver, got method=%r (fell back)' % ret.method


def _assert_close(name, got, ref):
    got = np.atleast_2d(np.asarray(got, dtype=float))
    ref = np.atleast_2d(np.asarray(ref, dtype=float))
    assert got.shape == ref.shape, '%s: shape %s != %s' % (name, got.shape, ref.shape)
    diff = np.abs(got - ref)
    tol = np.maximum(ABS_TOL, REL_TOL * np.abs(ref))
    bad = diff > tol
    assert not np.any(bad), \
        '%s mismatch vs CTMC:\n got=\n%s\n ref=\n%s\n diff=\n%s' % (name, got, ref, diff)


def _check_vs_ctmc(model, seed):
    ctmc = SolverCTMC(model)
    qc = np.asarray(ctmc.getAvgQLen(), dtype=float)
    tc = np.asarray(ctmc.getAvgTput(), dtype=float)

    solver = SolverSSA(model, 'nrm', seed=seed, samples=SAMPLES, verbose=False)
    qs = np.asarray(solver.getAvgQLen(), dtype=float)
    ts = np.asarray(solver.getAvgTput(), dtype=float)

    _assert_close('QLen', qs, qc)
    _assert_close('Tput', ts, tc)


def test_polling_exhaustive_no_switchover():
    model = _build(PollingType.EXHAUSTIVE, switchover=False)
    _assert_nrm_path(model)
    _check_vs_ctmc(model, seed=12345)


def test_polling_exhaustive_switchover():
    model = _build(PollingType.EXHAUSTIVE, switchover=True)
    _assert_nrm_path(model)
    _check_vs_ctmc(model, seed=12346)


def test_polling_gated_no_switchover():
    model = _build(PollingType.GATED, switchover=False)
    _assert_nrm_path(model)
    _check_vs_ctmc(model, seed=12347)


def test_polling_gated_switchover():
    model = _build(PollingType.GATED, switchover=True)
    _assert_nrm_path(model)
    _check_vs_ctmc(model, seed=12348)


def test_polling_klimited_switchover():
    model = _build(PollingType.KLIMITED, switchover=True, par=2)
    _assert_nrm_path(model)
    _check_vs_ctmc(model, seed=12349)


def test_polling_decrementing_switchover():
    model = _build(PollingType.DECREMENTING, switchover=True)
    _assert_nrm_path(model)
    _check_vs_ctmc(model, seed=12350)


if __name__ == '__main__':
    test_imports_worktree_line_solver()
    for fn in (test_polling_exhaustive_no_switchover,
               test_polling_exhaustive_switchover,
               test_polling_gated_no_switchover,
               test_polling_gated_switchover,
               test_polling_klimited_switchover,
               test_polling_decrementing_switchover):
        fn()
        print('%s passed.' % fn.__name__)
    print('All NRM polling tests passed.')
