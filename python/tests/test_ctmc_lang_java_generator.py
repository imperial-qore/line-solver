"""lang='java' regressions for the CTMC state-probability API.

Two things the JSON dispatch has to get right and used to get wrong:

* the (generator, state space, aggregated state space, pi) group must be ONE
  triple. On a reducible chain the analyzer restricts the generator to the
  component supporting pi, and the JAR reports both the full chain and the
  restricted ("work") one. Adopting the restricted state space next to the full
  generator produced a result on which ``pi @ Q`` could not even be formed.
* ``getProb(station)`` had no JAR route at all, so it fell through to the native
  path and died on ``_result.station_col_ranges``, which only the native
  analyzer produces.

Both assertions are self-referential (a balance equation, and native as the
reference computed in the same test) rather than hardcoded numbers.

NOTE for anyone extending this file: the two codebases enumerate the state space
in DIFFERENT ROW ORDERS, and on some models they do not even build the same
chain (native eliminates immediate states by stochastic complementation where
the JAR keeps them with a large rate). Never compare pi element-by-element
across langs; join on the state vector, as ``_pi_by_state`` does here.
"""

import os

import numpy as np
import pytest

import line_solver
# The worktree copy of line_solver must be the one under test; a global install
# would silently validate the wrong code.
assert os.path.realpath(__file__).rsplit('/python/', 1)[0] in os.path.realpath(
    line_solver.__file__), (
    "test must import the worktree line_solver, got %s" % line_solver.__file__)

from line_solver import (Network, Queue, Delay, ClosedClass, Exp, SchedStrategy,
                         SolverCTMC)

TOL = 1e-9


def _require_jar():
    """Skip unless jline.jar can be located: lang='java' needs it."""
    try:
        from line_solver.solvers.jar_dispatch import find_jar
        find_jar()
    except Exception as exc:            # pragma: no cover - environment guard
        pytest.skip('lang=java needs jline.jar: %s' % exc)


def _reducible_pas():
    """Closed pass-and-swap tandem: REDUCIBLE, so the analyzer restricts the
    generator to the component reachable from the declared placement. Three
    classes, one job each, swap graph 1-2, 2-3."""
    G = np.zeros((3, 3))
    for a, b in [(1, 2), (2, 3)]:
        G[a - 1, b - 1] = 1
        G[b - 1, a - 1] = 1
    m = Network('pas_closed_tandem')
    q1 = Queue(m, 'PASQueue1', SchedStrategy.PAS)
    q2 = Queue(m, 'PASQueue2', SchedStrategy.PAS)
    jc = [ClosedClass(m, 'Class%d' % (r + 1), 1, q1) for r in range(3)]
    q1.setService(lambda c: 1.0)
    q2.setService(lambda c: 1.3)
    q1.setSwapGraph(G); q1.setNumberOfServers(1); q1.setCap(3)
    q2.setSwapGraph(G); q2.setNumberOfServers(1); q2.setCap(3)
    m.addLink(q1, q2); m.addLink(q2, q1)
    for r in range(3):
        q1.setProbRouting(jc[r], q2, 1.0)
        q2.setProbRouting(jc[r], q1, 1.0)
    m.link(m.initRoutingMatrix())
    # A closed PAS network needs an explicit placement: it selects the recurrent
    # component the stationary distribution lives on.
    q1.setState(np.array([1, 2, 3]))
    return m


def _closed_cqn():
    m = Network('cqn_delay_queue')
    d = Delay(m, 'Delay')
    q = Queue(m, 'Queue1', SchedStrategy.FCFS)
    c = ClosedClass(m, 'Class1', 2, d, 0)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(3.0))
    P = m.initRoutingMatrix()
    P.set(c, c, d, q, 1.0)
    P.set(c, c, q, d, 1.0)
    m.link(P)
    return m


def _solve(model, lang, **kwargs):
    s = SolverCTMC(model, lang=lang, **kwargs)
    s.getAvgTable()
    return s


def _pi_by_state(solver):
    """{state row -> probability}. Keyed by the STATE VECTOR because the row
    order of the state space differs between the codebases."""
    pi = np.asarray(solver.getSteadyState(), dtype=float).reshape(-1)
    space = np.atleast_2d(np.asarray(solver._result.space, dtype=float))
    assert space.shape[0] == pi.size, \
        'pi (%d) and its state space (%s) disagree' % (pi.size, space.shape)
    return dict((tuple(np.round(row, 9)), p) for row, p in zip(space, pi))


def test_java_generator_triple_is_self_consistent_on_reducible_chain():
    """generator, state space and pi come back as one triple under lang='java'."""
    _require_jar()
    java = _solve(_reducible_pas(), 'java', cutoff=3)

    pi = np.asarray(java.getSteadyState(), dtype=float).reshape(-1)
    Q = np.atleast_2d(np.asarray(java.getInfGen(), dtype=float))
    space = np.atleast_2d(np.asarray(java._result.space, dtype=float))
    space_aggr = np.atleast_2d(np.asarray(java._result.space_aggr, dtype=float))

    assert pi.size > 0, 'no stationary distribution returned'
    assert Q.shape[0] == Q.shape[1] == pi.size, \
        'generator %s does not match pi (%d)' % (Q.shape, pi.size)
    assert space.shape[0] == pi.size, \
        'state space %s does not match pi (%d)' % (space.shape, pi.size)
    assert space_aggr.shape[0] == pi.size, \
        'aggregated state space %s does not match pi (%d)' % (space_aggr.shape, pi.size)

    # The defining property: pi solves the balance equations of the generator it
    # is returned with. This is what a full-space Q beside a restricted pi broke.
    assert abs(pi.sum() - 1.0) < TOL, 'pi does not sum to one: %g' % pi.sum()
    assert np.max(np.abs(pi @ Q)) < TOL, \
        'pi is not stationary for the returned generator: ||pi Q||_inf = %g' \
        % np.max(np.abs(pi @ Q))

    # And it is the same distribution the native solver computes, matched on the
    # state vector rather than on the row index.
    native = _solve(_reducible_pas(), 'python', cutoff=3)
    dj, dn = _pi_by_state(java), _pi_by_state(native)
    assert set(dj) == set(dn), \
        'state spaces differ: %d java rows, %d native rows, %d shared' \
        % (len(dj), len(dn), len(set(dj) & set(dn)))
    for state in dn:
        assert abs(dj[state] - dn[state]) < TOL, \
            'pi differs at state %s: java %g native %g' % (state, dj[state], dn[state])


def test_java_getprob_station_matches_native():
    """getProb(station) is answerable under lang='java'."""
    _require_jar()
    model_j, model_n = _closed_cqn(), _closed_cqn()
    java = _solve(model_j, 'java', cutoff=3, seed=1)
    native = _solve(model_n, 'python', cutoff=3, seed=1)

    for ist in range(len(model_n.getStations())):
        got = float(java.getProb(model_j.getStations()[ist]))
        ref = float(native.getProb(model_n.getStations()[ist]))
        assert abs(got - ref) < TOL, \
            'getProb(station %d): java %g vs native %g' % (ist, got, ref)
        assert 0.0 <= got <= 1.0, 'getProb(station %d) is not a probability: %g' % (ist, got)
