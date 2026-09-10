"""An unseeded reducible CTMC must say that it is averaging the recurrent classes.

A 2-station cyclic closed network with an order-preserving discipline (LCFS-PR)
and two classes cannot let jobs overtake, so the relative service order is frozen
at time zero and the ordered chain splits into two closed communicating classes.
The stationary distribution is then not unique and is fixed only by the declared
initial state; the solver restricts the chain to the class that state reaches.

When the initial state cannot be located in the enumerated space, that pruning
cannot happen and the block decomposition falls back to weighting the recurrent
classes EQUALLY. Nothing in the model implies that weighting -- it is the answer
to no question -- so it must be announced rather than returned silently.

The mixture value 63/110 is not an arbitrary recording: it is exactly the number
python returned before the seed lookup was fixed (see G60 in line-gaps.md),
against MATLAB/JAR 3/5 for the declared order. Pinning it here keeps the
diagnostic tied to the value it exists to explain.
"""
import numpy as np
import pytest

from line_solver import (Network, Queue, ClosedClass, Exp, SchedStrategy,
                         SolverCTMC)
from line_solver.api.solvers.ctmc import handler as ctmc_handler

MU = [[2.0, 3.0], [1.5, 4.0]]
CLASS_A = 3.0 / 5.0        # the class the declared initial state reaches
MIXTURE = 63.0 / 110.0     # (3/5 + 6/11)/2, the equal-weight fallback
TOL = 1e-9


def _model():
    m = Network('lcfspr_cyclic')
    q1 = Queue(m, 'Q1', SchedStrategy.LCFSPR)
    q2 = Queue(m, 'Q2', SchedStrategy.LCFSPR)
    c1 = ClosedClass(m, 'C1', 1, q1, 0)
    c2 = ClosedClass(m, 'C2', 1, q1, 0)
    for j, q in enumerate((q1, q2)):
        q.setService(c1, Exp(MU[j][0]))
        q.setService(c2, Exp(MU[j][1]))
    P = m.initRoutingMatrix()
    for c in (c1, c2):
        P.set(c, c, q1, q2, 1.0)
        P.set(c, c, q2, q1, 1.0)
    m.link(P)
    return m


def _tput(model, lang=None):
    kw = {} if lang is None else {'lang': lang}
    return float(np.asarray(SolverCTMC(model, **kw).getAvgTput(),
                            dtype=float).ravel()[0])


def _tput_and_warnings(model, monkeypatch):
    """Solve, returning (class-1 throughput, warning texts). The warnings are
    collected at the logger, not from stdout: the logger binds its sink once at
    construction, so a stream-level capture does not see them.

    The solve is pinned to lang='python' because both the patched seed lookup and
    the recording sink are native-python objects: under LINE_SOLVER_LANG=java the
    solve leaves the process and neither patch is reached, so the assertions
    would be about a chain this test never blinded."""
    from line_solver.api.io import logging as line_logging
    seen = []

    def _record(caller, msg, *args):
        seen.append(msg % args if args else msg)

    monkeypatch.setattr(line_logging, 'line_warning_always', _record)
    return _tput(model, lang='python'), seen


def test_seeded_solve_returns_the_class_the_initial_state_selects():
    assert _tput(_model()) == pytest.approx(CLASS_A, abs=TOL)


def test_unseeded_reducible_solve_warns_that_it_averages_the_classes(monkeypatch):
    """With the seed lookup blinded the answer becomes the equal-weight mixture,
    and the solver must say so."""
    monkeypatch.setattr(ctmc_handler, '_initial_state_index',
                        lambda sn, space: None)
    X, warnings = _tput_and_warnings(_model(), monkeypatch)
    assert X == pytest.approx(MIXTURE, abs=1e-9), \
        'expected the equal-weight mixture, got %r' % X
    assert any('closed communicating classes' in w for w in warnings), \
        'the unseeded reducible solve returned a mixture with no diagnostic: %r' % warnings


def test_irreducible_model_is_not_flagged(monkeypatch):
    """FCFS on the same topology lets a job lap another, so the ordered chain is
    irreducible and no diagnostic is due: the warning must not fire on it."""
    m = Network('fcfs_cyclic')
    q1 = Queue(m, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Q2', SchedStrategy.FCFS)
    c1 = ClosedClass(m, 'C1', 1, q1, 0)
    c2 = ClosedClass(m, 'C2', 1, q1, 0)
    for j, q in enumerate((q1, q2)):
        q.setService(c1, Exp(MU[j][0]))
        q.setService(c2, Exp(MU[j][1]))
    P = m.initRoutingMatrix()
    for c in (c1, c2):
        P.set(c, c, q1, q2, 1.0)
        P.set(c, c, q2, q1, 1.0)
    m.link(P)
    _, warnings = _tput_and_warnings(m, monkeypatch)
    assert not any('closed communicating' in w for w in warnings), warnings


def test_count_bscc_counts_only_classes_with_no_way_out():
    """One transient state feeding two absorbing cycles: two BSCCs, not three."""
    # states 0 (transient) -> {1,2} cycle A, {3,4} cycle B
    Q = np.zeros((5, 5))
    Q[0, 1] = 1.0
    Q[0, 3] = 2.0
    Q[1, 2] = Q[2, 1] = 1.0
    Q[3, 4] = Q[4, 3] = 1.0
    from scipy.sparse import csc_matrix
    from scipy.sparse.csgraph import connected_components
    nscc, labels = connected_components(csc_matrix(Q > 0), directed=True,
                                        connection='strong', return_labels=True)
    assert ctmc_handler._count_bscc(Q, labels, nscc) == 2
