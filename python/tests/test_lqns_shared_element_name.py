"""
The .lqxo parser must key a result row by NAME AND KIND, not by name alone.

A LINE-generated layered model routinely gives a processor, its task and that
task's entry the SAME name, and lqn.names holds all three. A name-only lookup
returns whichever element it meets first, so the .lqxo rows of the other two are
written into it and one result file yields three different wrong answers. The
kind is not ambiguous in the document -- it is the tag being read -- so the
parser has it available at every assignment.

The probe builds one model twice, once with unique names and once with three
elements sharing a name. The physics is identical, so every element must get the
same numbers under both. Under the name-only lookup the processor row carries
the task's throughput and the task and entry rows come back NaN.
"""

import numpy as np
import pytest

from line_solver import (Activity, Entry, Exp, LayeredNetwork, Processor,
                         SchedStrategy, SolverLQNS, Task)

pytestmark = pytest.mark.skipif(not SolverLQNS.isAvailable(),
                                reason='no lqns binary on the PATH')


def _build(proc_name, task_name, entry_name):
    m = LayeredNetwork('lqn_shared_name')
    p1 = Processor(m, 'P1', 1, SchedStrategy.PS)
    p2 = Processor(m, proc_name, 1, SchedStrategy.PS)
    t1 = Task(m, 'T1', 5, SchedStrategy.REF).on(p1).set_think_time(Exp(1 / 2))
    t2 = Task(m, task_name, 5, SchedStrategy.FCFS).on(p2).set_think_time(Exp(1 / 3))
    e1 = Entry(m, 'E1').on(t1)
    e2 = Entry(m, entry_name).on(t2)
    Activity(m, 'AS1', Exp(2)).on(t1).bound_to(e1).synch_call(e2, 1)
    Activity(m, 'AS2', Exp(4)).on(t2).bound_to(e2).replies_to(e2)
    return m


def _solve(model):
    """Per element, keyed by (kind, position within kind), the parsed rows."""
    solver = SolverLQNS(model, method='lqns')
    res = solver.getAvg()
    lqn = solver.getStruct()
    rows = {}
    for i in range(int(lqn.nidx)):
        kind = int(lqn.type[i])
        rows.setdefault(kind, []).append(
            (float(res.UN[i]), float(res.TN[i]), float(res.PN[i])))
    return rows


def test_shared_element_name_does_not_merge_result_rows():
    unique = _solve(_build('P2', 'T2', 'E2'))
    shared = _solve(_build('c0', 'c0', 'c0'))

    assert set(unique.keys()) == set(shared.keys())
    for kind in unique:
        assert len(unique[kind]) == len(shared[kind])
        for j, (a, b) in enumerate(zip(unique[kind], shared[kind])):
            for va, vb in zip(a, b):
                assert (np.isnan(va) and np.isnan(vb)) or va == pytest.approx(vb, abs=1e-9), (
                    'element %d of kind %d disagrees between the unique-name and the '
                    'shared-name model: %r vs %r' % (j, kind, a, b))
