"""AMVA must not solve a multiserver model as if every station had one server.

`pfqn_linearizer` and `pfqn_bs` carry no server-count argument, so a multiserver
model has to be transformed or rerouted before it reaches them. MATLAB
`solver_amva.m` does both: the product-form arm applies Seidmann's transform
(:111-117, L -> L/m with L(m-1)/m folded into the think time) and the linearizer
family is handed to `solver_amvald` when max(nservers)>1 (:243-253). Python's
solver-level dispatch did neither -- it computed `max_servers` only inside its
open-class arm -- so `lin` and `bs` returned the SINGLE-SERVER model's numbers,
bit for bit, on a multiserver model: 24% and 26% off exact where the JAR was 2.4%
and 6.8%. Not an approximation error; a different model.

The bit-identity is what makes this testable without a golden: a method that
honours `nservers` cannot return the same numbers for m=[1,2,3] and m=[1,1,1].
"""
import numpy as np
import pytest

from line_solver import (Network, Queue, Delay, ClosedClass, Exp, SchedStrategy,
                         SolverMVA)

D = np.array([[0.6, 0.4], [0.9, 0.3], [0.5, 0.8]])
N = [3, 2]
MULTI = [1, 2, 3]
SINGLE = [1, 1, 1]


def _model(nserv):
    m = Network('ms_amva')
    d = Delay(m, 'Think')
    qs = []
    for i in range(len(nserv)):
        q = Queue(m, 'Q%d' % (i + 1), SchedStrategy.PS)
        q.setNumberOfServers(int(nserv[i]))
        qs.append(q)
    cls = [ClosedClass(m, 'C%d' % (r + 1), N[r], d, 0) for r in range(2)]
    for r, c in enumerate(cls):
        d.setService(c, Exp(1.0))
        for i in range(len(nserv)):
            qs[i].setService(c, Exp(1.0 / D[i, r]))
    nodes = [d] + qs
    P = m.initRoutingMatrix()
    for c in cls:
        for k in range(len(nodes)):
            P.set(c, c, nodes[k], nodes[(k + 1) % len(nodes)], 1.0)
    m.link(P)
    return m


def _tput(nserv, method, **kwargs):
    s = SolverMVA(_model(nserv), method=method, **kwargs)
    s.getAvgTable()
    return np.asarray(s.getAvgTput(), dtype=float)[0, :]


@pytest.mark.parametrize('method', ['lin', 'bs', 'qd'])
def test_multiserver_answer_differs_from_the_single_server_one(method):
    """The regression that was silent: identical numbers meant nservers was dropped."""
    multi = _tput(MULTI, method)
    single = _tput(SINGLE, method)
    assert not np.allclose(multi, single, rtol=1e-9, atol=0), \
        ("method '%s' returned the single-server result %s for nservers=%s: "
         "the server counts were ignored" % (method, single, MULTI))


@pytest.mark.parametrize('method,tol', [('lin', 0.05), ('qd', 0.05), ('bs', 0.10)])
def test_multiserver_amva_is_close_to_exact(method, tol):
    """Exact load-dependent MVA is affordable at this size and is the reference."""
    exact = _tput(MULTI, 'exact')
    approx = _tput(MULTI, method)
    err = float(np.max(np.abs(approx - exact) / np.abs(exact)))
    assert err < tol, "method '%s' is %.1f%% off exact %s (got %s)" % (
        method, 100 * err, exact, approx)


def test_seidmann_transform_matches_the_hand_computation():
    """L -> L/m, with L(m-1)/m charged to the think time from the ORIGINAL demands."""
    L = np.array([[1.0, 2.0], [4.0, 8.0]])
    Z = np.array([0.5, 0.25])
    Lms, Zms = SolverMVA._amva_seidmann(L, Z, np.array([2.0, 4.0]))
    assert np.allclose(Lms, np.array([[0.5, 1.0], [1.0, 2.0]]))
    # 1*(1/2) + 4*(3/4) = 3.5 and 2*(1/2) + 8*(3/4) = 7.0
    assert np.allclose(Zms, np.array([0.5 + 3.5, 0.25 + 7.0]))


def test_single_server_model_is_untouched_by_the_transform():
    """The gate is max(nservers)>1, so a single-server model must take the old path."""
    L = np.array([[1.0, 2.0]])
    Z = np.array([0.5, 0.25])
    Lms, Zms = SolverMVA._amva_seidmann(L, Z, np.array([1.0]))
    assert np.allclose(Lms, L) and np.allclose(Zms, Z)


def test_conway_rule_uses_the_conway_multiserver_linearizer():
    """It reaches pfqn_conwayms, which takes nservers. The rule was refused while
    that routine divided by empty state sums (fixed 2026-07-31, see
    _kb/06-solver-catalog.md); dispatching it again is only safe because MATLAB,
    the JAR and python now return the same numbers for it."""
    x = _tput(MULTI, 'lin', config={'multiserver': 'conway'})
    single = _tput(SINGLE, 'lin')
    assert np.all(np.isfinite(x))
    assert not np.allclose(x, single, rtol=1e-9, atol=0)
    # MATLAB solver_amva.m with the same rule. Pinned as a golden rather than as
    # an accuracy bound: these stations are PS while Conway's correction is
    # derived for multiserver FCFS, so the rule is 6.7% off exact here (vs 1.6%
    # for krzesinski) in EVERY codebase. Agreement is the property under test.
    assert np.allclose(x, np.array([0.783226626, 0.626726626]), atol=1e-9)


def test_krzesinski_rule_uses_the_multiserver_linearizer():
    """It reaches pfqn_linearizermx, which takes nservers, so it must beat the
    single-server answer rather than reproduce it."""
    exact = _tput(MULTI, 'exact')
    x = _tput(MULTI, 'lin', config={'multiserver': 'krzesinski'})
    single = _tput(SINGLE, 'lin')
    assert not np.allclose(x, single, rtol=1e-9, atol=0)
    assert float(np.max(np.abs(x - exact) / np.abs(exact))) < 0.05
