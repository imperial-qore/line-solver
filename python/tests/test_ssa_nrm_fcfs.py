"""
Native-Python tests for the dedicated NRM SSA solver with FCFS scheduling.

Mirror of the JAR ``SolverSSATest`` FCFS NRM cases
(``testSolverSsaNrm_Fcfs_Open`` / ``_Closed`` / ``_Mixed``). The NRM solver
results are checked against the analytical SolverMVA ground truth within a
loose tolerance that accommodates simulation noise, and each case also asserts
that the dedicated NRM handler (not the basic-SSA fallback) was actually used.
"""

import numpy as np

from line_solver import (
    Network, Source, Sink, Delay, Queue, Exp,
    OpenClass, ClosedClass, SchedStrategy, SolverSSA, SolverMVA,
)

SAMPLES = 200000
REL_TOL = 0.06   # 6% relative tolerance for non-trivial metrics
ABS_TOL = 0.03   # absolute floor for near-zero metrics


def _assert_close(name, got, ref):
    got = np.atleast_2d(np.asarray(got, dtype=float))
    ref = np.atleast_2d(np.asarray(ref, dtype=float))
    assert got.shape == ref.shape, f'{name}: shape {got.shape} != {ref.shape}'
    diff = np.abs(got - ref)
    tol = np.maximum(ABS_TOL, REL_TOL * np.abs(ref))
    bad = diff > tol
    assert not np.any(bad), (
        f'{name} mismatch vs MVA:\n got=\n{got}\n ref=\n{ref}\n diff=\n{diff}')


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
        f'expected dedicated NRM solver, got method={ret.method!r} (fell back)'


def _check_vs_mva(model, seed, rows=None):
    """Compare NRM-SSA mean metrics to MVA. ``rows`` (0-based station indices)
    restricts the comparison; the Source/Sink generator rows of open models
    carry a degenerate (and, in the cyclic open representation, transiently
    negative) population in the raw getters and are excluded by passing the
    queueing station rows only."""
    solver = SolverSSA(model, 'nrm', seed=seed, samples=SAMPLES, verbose=False)
    qn, un, tn = solver.getAvgQLen(), solver.getAvgUtil(), solver.getAvgTput()

    mva = SolverMVA(model)
    mva.getAvgTable()
    mqn, mun, mtn = mva.getAvgQLen(), mva.getAvgUtil(), mva.getAvgTput()

    if rows is not None:
        rows = np.asarray(rows, dtype=int)
        qn, un, tn = qn[rows], un[rows], tn[rows]
        mqn, mun, mtn = mqn[rows], mun[rows], mtn[rows]

    _assert_close('QLen', qn, mqn)
    _assert_close('Util', un, mun)
    _assert_close('Tput', tn, mtn)


def test_ssa_nrm_fcfs_open():
    """NRM with FCFS on a two-class open M/M/1 network."""
    model = Network('open')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue1', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oc1 = OpenClass(model, 'Class1')
    oc2 = OpenClass(model, 'Class2')
    source.setArrival(oc1, Exp(0.4))
    source.setArrival(oc2, Exp(0.3))
    queue.setService(oc1, Exp(2.0))
    queue.setService(oc2, Exp(1.5))
    P = model.initRoutingMatrix()
    P[oc1, oc1] = Network.serialRouting(source, queue, sink)
    P[oc2, oc2] = Network.serialRouting(source, queue, sink)
    model.link(P)

    _assert_nrm_path(model)
    # Compare the queueing station only; Source/Sink rows are degenerate.
    _check_vs_mva(model, seed=23000, rows=[queue.getStationIndex() - 1])


def test_ssa_nrm_fcfs_closed():
    """NRM with FCFS on a single-class closed cyclic network."""
    model = Network('closed')
    delay = Delay(model, 'Think')
    queue = Queue(model, 'Q', SchedStrategy.FCFS)
    cc = ClosedClass(model, 'C', 3, delay)
    delay.setService(cc, Exp(1.0))
    queue.setService(cc, Exp(2.0))
    model.link(Network.serialRouting(delay, queue))

    _assert_nrm_path(model)
    _check_vs_mva(model, seed=24000)


def test_ssa_nrm_fcfs_multiserver():
    """NRM with FCFS on a closed multiserver (c=2) network."""
    model = Network('closed')
    delay = Delay(model, 'Think')
    queue = Queue(model, 'Q', SchedStrategy.FCFS)
    queue.setNumberOfServers(2)
    cc = ClosedClass(model, 'C', 4, delay)
    delay.setService(cc, Exp(1.0))
    queue.setService(cc, Exp(1.5))
    model.link(Network.serialRouting(delay, queue))

    _assert_nrm_path(model)
    _check_vs_mva(model, seed=25000)


if __name__ == '__main__':
    test_ssa_nrm_fcfs_open()
    test_ssa_nrm_fcfs_closed()
    test_ssa_nrm_fcfs_multiserver()
    print('All NRM FCFS tests passed.')
