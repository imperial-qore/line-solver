"""Round-trip of the LINE .lqnx dialect: routed call groups.

`synch_call_rrobin` / `synch_call_jsq` state that one dispatcher issues the
activity's calls over a set of target entries. The stock LQN schema has no
element for that, so an export used to emit a valid file describing a DIFFERENT
model: three independent coins in place of a dispatcher, at the same aggregate
call rate. Nothing downstream could tell the two apart, because the call means
are exactly what a group preserves.

The member calls stay ordinary <synch-call> elements and the grouping alone
rides the extra <call-group> element. That is what keeps the file readable by
lqns and lqsim, which have no dispatcher and analyse the ungrouped twin, and it
is also the trap: a reader that re-issued the calls from the group element would
double the call rate.
"""
import os
import tempfile

import pytest

from line_solver import (Activity, Entry, Exp, LayeredNetwork, Processor,
                         SchedStrategy, Task)
from line_solver.constants import RoutingStrategy


def _group_model(kind):
    """Client dispatching one call per invocation over three identical servers.

    KIND is 'rrobin', 'jsq', or 'plain' for the ungrouped twin with the same
    per-target means.
    """
    m = LayeredNetwork('cgRT')
    pc = Processor(m, 'PC', 1, SchedStrategy.INF)
    ps = Processor(m, 'PS', 1, SchedStrategy.PS)
    tc = Task(m, 'TC', 10, SchedStrategy.REF).on(pc).set_think_time(Exp(1 / 5))
    ts = [Task(m, 'TS%d' % i, 5, SchedStrategy.FCFS).on(ps) for i in (1, 2, 3)]
    ec = Entry(m, 'EC').on(tc)
    es = [Entry(m, 'ES%d' % i).on(ts[i - 1]) for i in (1, 2, 3)]
    ac = Activity(m, 'AC', Exp(2)).on(tc).bound_to(ec)
    if kind == 'rrobin':
        ac.synch_call_rrobin(es, 1.0)
    elif kind == 'jsq':
        ac.synch_call_jsq(es, 1.0)
    else:
        for e in es:
            ac.synch_call(e, 1.0 / 3)
    for i in (1, 2, 3):
        Activity(m, 'AS%d' % i, Exp(1)).on(ts[i - 1]).bound_to(es[i - 1]).replies_to(es[i - 1])
    return m


def _roundtrip(model):
    path = os.path.join(tempfile.mkdtemp(prefix='lqnx_cg_'), 'm.lqnx')
    model.writeXML(path)
    return path, LayeredNetwork.parse_xml(path)


def _groups(model):
    """(caller name, strategy, target names) per group, off the model objects."""
    out = []
    for act in model.activities:
        for strategy, entries in getattr(act, 'call_groups', []):
            out.append((act.name, strategy, [e.name for e in entries]))
    return out


@pytest.mark.parametrize('kind,strategy', [('rrobin', RoutingStrategy.RROBIN),
                                           ('jsq', RoutingStrategy.JSQ)])
def test_group_survives_roundtrip(kind, strategy):
    src = _group_model(kind)
    _, back = _roundtrip(src)
    assert _groups(back) == _groups(src)
    assert _groups(back) == [('AC', strategy, ['ES1', 'ES2', 'ES3'])]


def test_member_calls_are_not_issued_twice():
    """The group element records the grouping, not the calls."""
    src = _group_model('rrobin')
    path, back = _roundtrip(src)
    before, after = src.getStruct(), back.getStruct()
    assert after.ncalls == before.ncalls == 3
    means = [after.callproc[c].getMean() for c in range(after.ncalls)]
    assert sum(means) == pytest.approx(1.0)

    # and the document says so: the members are ordinary synch-calls
    text = open(path).read()
    assert text.count('<synch-call') == 3
    assert text.count('<call-group') == 1


def test_group_indices_survive_into_the_struct():
    src = _group_model('jsq')
    _, back = _roundtrip(src)
    before, after = src.getStruct(), back.getStruct()
    assert len(after.callgroups) == len(before.callgroups) == 1
    (cb, sb, tb), (ca, sa, ta) = before.callgroups[0], after.callgroups[0]
    assert before.hashnames[cb] == after.hashnames[ca]
    assert sb == sa
    assert [before.hashnames[i] for i in tb] == [after.hashnames[i] for i in ta]


def test_ungrouped_twin_reads_back_with_no_group():
    """The same call means with no dispatcher must not acquire one."""
    src = _group_model('plain')
    path, back = _roundtrip(src)
    assert _groups(back) == []
    assert open(path).read().count('<call-group') == 0
    assert back.getStruct().ncalls == 3


def test_unknown_strategy_is_refused_by_name():
    src = _group_model('rrobin')
    path, _ = _roundtrip(src)
    text = open(path).read().replace('strategy="RROBIN"', 'strategy="PROB"')
    with open(path, 'w') as fd:
        fd.write(text)
    with pytest.raises(ValueError, match='RROBIN and JSQ'):
        LayeredNetwork.parse_xml(path)
