"""SolverAG: the agent-based (RCAT) solver and its execution backends.

Two things are pinned here, and they are pinned against CLOSED FORMS rather than
against solver output.

First, the solver boundary. RCAT decomposes the MODEL into cooperating agents;
every other MAM method decomposes its TRAFFIC. They share the shape of the answer
and no machinery at all, which is why the RCAT methods moved out of SolverMAM --
and why SolverMAM must now refuse them by name rather than silently resolve them,
and must refuse a G-network outright rather than solve it with every signal
turned into an ordinary customer.

Second, the execution backends. Agent k's generator is
Q_k(x) = L_k + sum_{c passive at k} x_c Pb_c, so an agent reads the rest of the
model only through the scalar reversed rates x and writes only its own slot. The
sweep is Jacobi, so the agent order is immaterial and 'parallel' must be
BIT-IDENTICAL to 'serial' -- not merely close. An approximate assertion here
would pass on a backend that had quietly started racing.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np
import pytest

from line_solver import (Network, Source, Queue, Sink, OpenClass, SchedStrategy,
                         Exp, SolverAG, SolverMAM)
from line_solver.solvers.solver_ag import SolverAGOptions


def _tandem(lam=1.0, mu1=2.0, mu2=3.0):
    model = Network('Tandem')
    source = Source(model, 'Source')
    q1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    cls = OpenClass(model, 'C')
    source.setArrival(cls, Exp(lam))
    q1.setService(cls, Exp(mu1))
    q2.setService(cls, Exp(mu2))
    P = model.initRoutingMatrix()
    P.set(cls, cls, source, q1, 1.0)
    P.set(cls, cls, q1, q2, 1.0)
    P.set(cls, cls, q2, sink, 1.0)
    model.link(P)
    return model


def _solve(method='inap', config=None, **kw):
    opts = SolverAGOptions(method=method)
    if config:
        opts.config = dict(config)
    s = SolverAG(_tandem(**kw), method=method, options=opts)
    s.runAnalyzer()
    return s.result


@pytest.mark.parametrize('lam,mu1,mu2', [(1.0, 2.0, 3.0), (0.5, 1.0, 4.0)])
def test_tandem_is_two_isolated_mm1_queues(lam, mu1, mu2):
    """Burke: the departure stream of an M/M/1 in equilibrium is Poisson at the
    arrival rate, so each queue is M/M/1 with QLen rho/(1-rho). RCAT reproduces
    it because the tandem IS a product form, which is the theorem's own case."""
    r = _solve(lam=lam, mu1=mu1, mu2=mu2)
    QN = np.asarray(r.QN, dtype=float)
    for row, mu in ((1, mu1), (2, mu2)):
        rho = lam / mu
        assert QN[row, 0] == pytest.approx(rho / (1.0 - rho), abs=1e-6)


def test_threads_backend_is_bit_identical_to_serial():
    """The property the whole parallel design rests on. Jacobi has no ordering,
    and 'parallel' runs the same agent code in the same process, so anything less
    than exact equality means an agent has started reading another agent's slot."""
    ref = _solve()
    got = _solve(config={'exec': 'parallel'})
    for name in ('QN', 'UN', 'RN', 'TN'):
        a = np.asarray(getattr(ref, name), dtype=float)
        b = np.asarray(getattr(got, name), dtype=float)
        assert np.array_equal(a, b), '%s differs by %.3e' % (name, np.max(np.abs(a - b)))


def test_threads_backend_honours_an_explicit_pool_size():
    """A pinned pool size must not change the answer either: the partition is
    over agents, and every agent is solved wherever it lands."""
    ref = _solve()
    for nworkers in (1, 2, 8):
        got = _solve(config={'exec': 'parallel', 'nworkers': nworkers})
        assert np.array_equal(np.asarray(ref.QN, dtype=float),
                              np.asarray(got.QN, dtype=float))


def test_unknown_execution_backend_is_refused_by_name():
    with pytest.raises(ValueError, match="Unknown AG execution backend"):
        _solve(config={'exec': 'gpu'})


def test_cluster_without_endpoints_is_refused_by_name():
    """A cluster run with nowhere to send the agents is a configuration error,
    not a silent local run: the caller asked for distribution and would otherwise
    be told nothing when they did not get it."""
    with pytest.raises(ValueError, match='needs worker endpoints'):
        _solve(config={'exec': 'cluster'})


def test_cluster_refuses_inapinf_rather_than_substituting_a_method():
    """The ag-worker carries the FINITE agent solve. 'inapinf' replaces it with
    the matrix-geometric tail of an open agent, which the worker does not have,
    so answering with the finite solve would silently change the method."""
    with pytest.raises(ValueError, match="does not carry the 'inapinf'"):
        _solve(method='inapinf',
               config={'exec': 'cluster', 'endpoints': ['localhost:1']})


def test_cluster_degrades_to_local_when_no_worker_answers():
    """A lost worker costs wall clock and nothing else: any agent can be solved
    anywhere given x, so an unreachable endpoint must still produce the run's
    answer rather than an error or a wrong number."""
    ref = _solve()
    with pytest.warns(UserWarning, match='AG worker'):
        got = _solve(config={'exec': 'cluster',
                             'endpoints': ['127.0.0.1:9'],
                             'worker_timeout': 2.0})
    assert np.array_equal(np.asarray(ref.QN, dtype=float),
                          np.asarray(got.QN, dtype=float))


def test_solver_mam_redirects_the_moved_rcat_methods():
    """A caller carrying an old options.method must be told where the method
    went, not that it is unknown."""
    for method in ('inap', 'inapplus', 'inapinf', 'exact'):
        with pytest.raises(RuntimeError, match='moved to SolverAG'):
            SolverMAM(_tandem(), method=method).runAnalyzer()


def test_the_g_network_features_belong_to_ag_alone():
    """The RCAT builder is the only code in LINE that reads sn.issignal, so no
    MAM method may declare the signal names: declaring them is what let a
    G-network reach a decomposition that ignores the marking."""
    ag = SolverAG.getFeatureSet()
    mam = SolverMAM.getFeatureSet()
    for name in ('OpenSignal', 'ClosedSignal', 'SignalType_NEGATIVE',
                 'SignalType_CATASTROPHE', 'SignalBatchRemoval'):
        assert name in ag
        assert name not in mam


def test_the_two_solvers_advertise_disjoint_method_sets():
    ag = set(SolverAG.listValidMethods())
    mam = set(SolverMAM.listValidMethods())
    assert ag & mam == {'default'}
    assert {'inap', 'inapplus', 'inapinf'} <= ag
    assert not ({'inap', 'inapplus', 'inapinf'} & mam)


def test_exact_reports_the_method_that_actually_ran():
    """'exact' is the vestigial AutoCAT alias and falls back to inap. It is
    classified globally as an EXACT method, so leaving the name on the result
    would banner an iterative approximation as exact."""
    assert _solve(method='exact').method == 'inap'
    assert _solve(method='default').method == 'inap'


def test_para_is_an_accepted_alias_and_threads_is_gone():
    """'para' resolves to 'parallel'; 'threads' is refused BY NAME.

    The alias pair is SolverSSA's ({'para','parallel'}), so one spelling
    convention covers both solvers. 'threads' was this backend's name until
    2026-08-19: it must not silently fall through to the generic
    unknown-backend message, or a caller with an old script is told a backend
    that still exists does not.
    """
    import pytest
    from line_solver.solvers.solver_ag import exec_backend

    assert exec_backend.create({'exec': 'para'}, 'inap').__class__ is \
        exec_backend.create({'exec': 'parallel'}, 'inap').__class__

    with pytest.raises(ValueError, match="renamed to 'parallel'"):
        exec_backend.create({'exec': 'threads'}, 'inap')
