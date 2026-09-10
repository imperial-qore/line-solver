"""SolverMVA must publish its iteration count and convergence flag in every lang.

The count answers "did this fixed point run out of budget", and it is the only
input to the non-convergence warning when the handler reports no flag. It used to
be unreadable in python: the native result dict had no `iter` key at all, and the
lang='java' delegation dropped both signals on the floor because the CLI payload
carried neither. A model that warns under one lang and is silent under the other
is worse than one that never warns, because the silence reads as convergence.

The flag is authoritative where present, and `None` is not `False`: on the
load-dependent route the counter aggregates the nested inner sweeps and saturates
the budget by construction on a solve whose outer residual is exactly zero, so
only the residual decides there. See G31(b) in line-gaps.md.
"""
import os

import numpy as np
import pytest

from line_solver import (Network, Queue, Delay, ClosedClass, Exp, SchedStrategy,
                         SolverMVA)

NSERV = [1, 2, 3]
D = np.array([[0.6, 0.4], [0.9, 0.3], [0.5, 0.8]])
N = [3, 2]


def _model():
    m = Network('mva_convergence')
    d = Delay(m, 'Think')
    qs = []
    for i in range(3):
        q = Queue(m, 'Q%d' % (i + 1), SchedStrategy.PS)
        q.setNumberOfServers(NSERV[i])
        qs.append(q)
    cls = [ClosedClass(m, 'C%d' % (r + 1), N[r], d, 0) for r in range(2)]
    for r, c in enumerate(cls):
        d.setService(c, Exp(1.0))
        for i in range(3):
            qs[i].setService(c, Exp(1.0 / D[i, r]))
    nodes = [d] + qs
    P = m.initRoutingMatrix()
    for c in cls:
        for k in range(len(nodes)):
            P.set(c, c, nodes[k], nodes[(k + 1) % len(nodes)], 1.0)
    m.link(P)
    return m


def _jar_available():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.exists(os.path.join(os.path.dirname(here), 'common', 'jline.jar'))


@pytest.mark.parametrize('method', ['lin', 'qd', 'bs'])
def test_native_result_publishes_an_iteration_count(method):
    # lang is named, not inherited: this asserts on the NATIVE result dict, whose
    # delegated counterpart is the _JavaAvgResult object the next test covers.
    s = SolverMVA(_model(), method=method, lang='python')
    s.getAvgTable()
    assert 'iter' in s._result, "the native result dict has no 'iter' key"
    assert 'converged' in s._result, "the native result dict has no 'converged' key"
    assert s._result['iter'] and s._result['iter'] > 0, \
        "method '%s' reported iter=%r" % (method, s._result['iter'])


@pytest.mark.parametrize('method', ['lin', 'qd', 'bs'])
def test_delegated_solve_carries_the_count_across_the_transport(method):
    if not _jar_available():
        pytest.skip('common/jline.jar not built')
    s = SolverMVA(_model(), method=method, lang='java')
    s.getAvgTable()
    assert getattr(s._result, 'iter', None) is not None, \
        "lang='java' returned no iteration count: the CLI payload dropped it"
    assert s._result.iter > 0, "method '%s' reported iter=%r" % (method, s._result.iter)


def test_delegated_solve_carries_the_convergence_flag():
    """The flag is what the warning keys on where the count is an aggregate."""
    if not _jar_available():
        pytest.skip('common/jline.jar not built')
    s = SolverMVA(_model(), method='lin', lang='java')
    s.getAvgTable()
    assert getattr(s._result, 'converged', None) is True, \
        'this model converges; the delegated solve reported %r' % getattr(s._result, 'converged', None)


def test_delegated_non_convergence_raises_the_same_warning(monkeypatch):
    """A delegated solve reporting converged=False must warn, as the native path does."""
    from line_solver.api.io import logging as line_logging
    seen = []
    monkeypatch.setattr(line_logging, 'line_warning_always',
                        lambda caller, msg, *a: seen.append(msg % a if a else msg))

    s = SolverMVA(_model(), method='lin')
    s._lastiter, s._lastconverged, s._lastiterbudget = 7, False, 1000
    s._warn_if_not_converged('lin')
    assert any('did not meet the convergence tolerance' in w for w in seen), seen

    seen.clear()
    s._lastconverged = True
    s._lastiter = 10 ** 6          # a saturated count must not override the flag
    s._warn_if_not_converged('lin')
    assert not seen, 'a converged solve warned anyway: %r' % seen
