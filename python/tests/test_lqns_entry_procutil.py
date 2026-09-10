"""
An entry's Util is the sum over ITS OWN activities, not the raw .lqxo attribute.

lqns credits host work to whichever level declares the host demand. In the
activity-graph (``task-activities``) form -- the only form LINE's writer ever
emits -- an entry declares none, so lqns writes a literal ``proc-utilization="0"``
on every ``result-entry`` and puts the work on the ``result-activity`` rows.
Read verbatim, every entry's Util came back 0 while the model plainly runs.

``lqn_twotasks`` is the case that pins the rule down: T2 hosts TWO entries, so
the entry sums (0.75 over A20/A21/A22, 0.25 over A3) are NOT the task's
proc-utilization (1.0). A fix that copied the task value, or that summed every
activity of the task into each of its entries, would pass a one-entry task and
fail here.

The assertions are physical, not transcribed: each entry's Util is its own
throughput times the host demand of the activities bound to it, which is the
quantity the .lqxo activity rows already carry.
"""

import pytest

from line_solver import (Activity, ActivityPrecedence, Entry, Exp, Erlang,
                         LayeredNetwork, Processor, SchedStrategy, SolverLQNS,
                         Task)

pytestmark = pytest.mark.skipif(not SolverLQNS.isAvailable(),
                                reason='no lqns binary on the PATH')


def _two_tasks():
    """T1 calls E2 and E3, both hosted by T2; E2 runs a three-activity chain."""
    model = LayeredNetwork('lqn_twotasks')
    p1 = Processor(model, 'P1', 1, SchedStrategy.PS)
    t1 = Task(model, 'T1', 100, SchedStrategy.REF).on(p1)
    e1 = Entry(model, 'E1').on(t1)

    p2 = Processor(model, 'P2', 1, SchedStrategy.PS)
    t2 = Task(model, 'T2', 1, SchedStrategy.INF).on(p2)
    e2 = Entry(model, 'E2').on(t2)
    e3 = Entry(model, 'E3').on(t2)

    t1.set_think_time(Erlang.fit_mean_and_order(10, 1))
    Activity(model, 'A1', Exp(1)).on(t1).bound_to(e1).synch_call(e2).synch_call(e3, 1)

    a20 = Activity(model, 'A20', Exp(1)).on(t2).bound_to(e2)
    a21 = Activity(model, 'A21', Exp(1)).on(t2)
    a22 = Activity(model, 'A22', Exp(1)).on(t2).replies_to(e2)
    t2.add_precedence(ActivityPrecedence.serial([a20, a21, a22]))

    Activity(model, 'A3', Exp(1)).on(t2).bound_to(e3).replies_to(e3)
    return model


def _rows(table):
    return {str(r.Node): r for r in table.itertuples()}


def test_entry_procutil_sums_its_own_activities():
    table = _rows(SolverLQNS(_two_tasks(), verbose=False).get_avg_table())

    # Every activity holds unit demand and the same throughput, so each
    # contributes Tput * 1.0 and an entry's Util is that times its chain length.
    for name in ('A1', 'A20', 'A21', 'A22', 'A3'):
        assert table[name].Util == pytest.approx(table[name].Tput, rel=1e-3)

    e2_expected = sum(table[a].Util for a in ('A20', 'A21', 'A22'))
    e3_expected = table['A3'].Util
    assert table['E2'].Util == pytest.approx(e2_expected, rel=1e-3)
    assert table['E3'].Util == pytest.approx(e3_expected, rel=1e-3)
    assert table['E1'].Util == pytest.approx(table['A1'].Util, rel=1e-3)

    # The regression itself: a verbatim read gives 0 on every entry.
    for name in ('E1', 'E2', 'E3'):
        assert table[name].Util > 0


def test_entry_procutil_is_not_the_task_value():
    """T2's two entries must split its proc-utilization, not each inherit it."""
    table = _rows(SolverLQNS(_two_tasks(), verbose=False).get_avg_table())

    t2, e2, e3 = table['T2'].Util, table['E2'].Util, table['E3'].Util
    assert e2 == pytest.approx(3 * e3, rel=1e-3)   # 3 activities against 1
    assert e2 + e3 == pytest.approx(t2, rel=1e-3)  # and together they are the task
    assert e2 != pytest.approx(t2, rel=1e-3)
    assert e3 != pytest.approx(t2, rel=1e-3)
