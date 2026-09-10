"""An entry lqns never invoked has a service time of zero, not an absent one.

lqns omits ``phase1-service-time`` from ``result-entry`` exactly when the
entry's throughput is zero -- nothing was served, so there is no per-invocation
mean to report. Read verbatim that is a NaN, and it landed in the RespT column
where the LQN table says an entry HAS a response time and every other solver
reports one. The NaN mask is part of the answer (see
``_kb/06-solver-catalog.md``), and a tolerance-based parity comparison cannot
see a break in it: ``compare_values`` passes any cell where either side is NaN.

The value is derived from the activity rows and ONLY where they are unanimous:
if every activity reachable from the entry reports a zero service time then
every aggregation law agrees on zero -- the serial sum, the branch-weighted mean
of an OrFork, the order statistic of an AndFork. It is deliberately not
generalised the way ``test_lqns_entry_procutil`` generalises utilization, which
DOES add over an activity graph; the last test below is what pins that
distinction down.

Twins: ``jar/src/test/java/jline/solvers/wrappers/lqns/SolverLQNSEntryServiceTimeTest.java``,
the ``lqns entry svct`` cases in ``cpp/tests/test_solver_lqns.cpp`` and
``line-test.git/test/testsMisc/test_lqns_entry_svct.m``.
"""

import math

import pytest

from line_solver import (Activity, ActivityPrecedence, Entry, Exp, LayeredNetwork,
                         Processor, SchedStrategy, SolverLQNS, SolverLN, SolverMVA,
                         Task)

pytestmark = pytest.mark.skipif(not SolverLQNS.isAvailable(),
                                reason='no lqns binary on the PATH')


def _unreachable_entry():
    """A working two-tier model plus a component no reference task reaches.

    T3 is not a reference task and nobody calls E3, so lqns solves it at zero
    throughput and writes a `result-entry` with no `phase1-service-time`. This is
    the shape `lqn_ofbiz` carries in its USAGE_DELAY component.
    """
    model = LayeredNetwork('lqn_unreachable_entry')
    p1 = Processor(model, 'P1', 1, SchedStrategy.INF)
    t1 = Task(model, 'T1', 1, SchedStrategy.REF).on(p1)
    e1 = Entry(model, 'E1').on(t1)
    p2 = Processor(model, 'P2', 1, SchedStrategy.INF)
    t2 = Task(model, 'T2', 1, SchedStrategy.INF).on(p2)
    e2 = Entry(model, 'E2').on(t2)
    t1.set_think_time(Exp(1.0))
    Activity(model, 'A1', Exp(1.0)).on(t1).bound_to(e1).synch_call(e2, 1)
    Activity(model, 'A2', Exp(1.0)).on(t2).bound_to(e2).replies_to(e2)

    p3 = Processor(model, 'P3', 1, SchedStrategy.INF)
    t3 = Task(model, 'T3', 1, SchedStrategy.FCFS).on(p3)
    e3 = Entry(model, 'E3').on(t3)
    # TWO activities in series, so every writer emits the activity-graph
    # (task-activities) form. A single bound activity is written as
    # entry-phase-activities by the JAR, where the reader already has a fallback
    # of its own and the omission never surfaces -- the twins share this fixture
    # so that all four exercise the same path.
    a3 = Activity(model, 'A3', Exp(1.0)).on(t3).bound_to(e3)
    a3b = Activity(model, 'A3b', Exp(1.0)).on(t3).replies_to(e3)
    t3.add_precedence(ActivityPrecedence.serial([a3, a3b]))
    return model


def _branching_entry():
    """An entry whose activity graph BRANCHES, so its service time is not a sum.

    E2 runs A20 and then an OrFork to A21 or A22. lqns reports the entry's
    phase1-service-time as the branch-weighted total; summing the activity rows
    would over-count the arm not taken.
    """
    model = LayeredNetwork('lqn_branching_entry')
    p1 = Processor(model, 'P1', 1, SchedStrategy.PS)
    t1 = Task(model, 'T1', 10, SchedStrategy.REF).on(p1)
    e1 = Entry(model, 'E1').on(t1)

    p2 = Processor(model, 'P2', 1, SchedStrategy.PS)
    t2 = Task(model, 'T2', 1, SchedStrategy.INF).on(p2)
    e2 = Entry(model, 'E2').on(t2)

    t1.set_think_time(Exp(0.1))
    Activity(model, 'A1', Exp(1)).on(t1).bound_to(e1).synch_call(e2, 1)

    a20 = Activity(model, 'A20', Exp(1)).on(t2).bound_to(e2)
    a21 = Activity(model, 'A21', Exp(1)).on(t2).replies_to(e2)
    a22 = Activity(model, 'A22', Exp(1)).on(t2).replies_to(e2)
    t2.add_precedence(ActivityPrecedence.or_fork(a20, [a21, a22], [0.5, 0.5]))
    return model


def _rows(table):
    return {str(r.Node): r for r in table.itertuples()}


def test_an_entry_lqns_never_invoked_reports_zero_not_nan():
    table = _rows(SolverLQNS(_unreachable_entry(), verbose=False).get_avg_table())

    # the unreachable component solves at zero throughput
    assert table['E3'].Tput == pytest.approx(0.0, abs=1e-12)
    assert table['A3'].Tput == pytest.approx(0.0, abs=1e-12)
    assert table['A3b'].Tput == pytest.approx(0.0, abs=1e-12)

    # the regression itself: RespT was NaN, because lqns omits the attribute
    assert not math.isnan(table['E3'].RespT), 'E3.RespT must not be NaN'
    assert table['E3'].RespT == pytest.approx(0.0, abs=1e-12)

    # the reachable entries are untouched and still carry lqns' own numbers
    assert table['E1'].RespT == pytest.approx(2.0, rel=1e-3)
    assert table['E2'].RespT == pytest.approx(1.0, rel=1e-3)


def test_it_agrees_with_ln_on_the_same_model():
    """SolverLN is the reference for which cells exist; LQNS must not differ."""
    lqns = _rows(SolverLQNS(_unreachable_entry(), verbose=False).get_avg_table())
    opts = SolverLN.default_options()
    opts.verbose = 0
    ln = _rows(SolverLN(_unreachable_entry(),
                        lambda net: SolverMVA(net, verbose=False), opts).get_avg_table())

    # ResidT is excluded by decision: lqns reports no residence time at all.
    for name in ('P3', 'T3', 'E3', 'A3', 'A3b'):
        for metric in ('QLen', 'Util', 'RespT', 'ArvR', 'Tput'):
            a = math.isnan(getattr(lqns[name], metric))
            b = math.isnan(getattr(ln[name], metric))
            assert a == b, '%s.%s: LQNS %s, LN %s' % (
                name, metric, 'NaN' if a else 'value', 'NaN' if b else 'value')


def test_the_derivation_does_not_touch_a_branching_entry():
    """Where lqns DOES report the attribute, its value survives verbatim.

    Utilizations add over an activity graph and response times do not, so the
    fallback must never become a sum: on this OrFork the sum over the entry's
    activities exceeds what lqns reports, and the reported value is the one that
    must come back.
    """
    table = _rows(SolverLQNS(_branching_entry(), verbose=False).get_avg_table())

    e2 = table['E2'].RespT
    assert not math.isnan(e2)
    summed = sum(table[a].RespT for a in ('A20', 'A21', 'A22'))
    # one arm of the fork is not taken, so the sum over-counts
    assert summed > e2 + 1e-6, ('the fixture no longer branches: sum %g against '
                                'reported %g' % (summed, e2))
