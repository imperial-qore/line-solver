"""
LQN features that LQN2QN must carry into the converted queueing network.

Both the AND-join quorum and the activity think time used to be read into the
LayeredNetworkStruct and then silently dropped by the conversion: a k-of-n join
became a wait-for-all Join, and a think time vanished instead of appearing as a
processor-releasing delay in series with the host demand.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import warnings

import pytest

from line_solver import (Activity, ActivityPrecedence, Entry, Exp, LayeredNetwork,
                         Processor, SchedStrategy, Task)
from line_solver.lang.base import JoinStrategy
from line_solver.io import LQN2QN


def fork_join_model(quorum, with_think):
    """Client -> Server, the server entry forking into three branches joined
    with the given quorum, optionally with a think time on the first branch."""
    m = LayeredNetwork('fj')
    p1 = Processor(m, 'P1', 1, SchedStrategy.INF)
    p2 = Processor(m, 'P2', 1, SchedStrategy.PS)
    client = Task(m, 'Client', 5, SchedStrategy.REF).on(p1)
    client.set_think_time(Exp(1.0))
    ce = Entry(m, 'CE').on(client)
    server = Task(m, 'Server', 3, SchedStrategy.FCFS).on(p2)
    se = Entry(m, 'SE').on(server)
    ca = Activity(m, 'ca', Exp(10.0)).on(client).bound_to(ce)
    ca.synch_call(se, 1.0)
    a0 = Activity(m, 'a0', Exp(1.0)).on(server).bound_to(se)
    branches = [Activity(m, 'b%d' % i, Exp(2.0)).on(server) for i in (1, 2, 3)]
    aj = Activity(m, 'aj', Exp(3.0)).on(server)
    if with_think:
        branches[0].set_think_time(Exp(4.0))
    server.add_precedence(ActivityPrecedence.AndFork(a0, branches))
    server.add_precedence(ActivityPrecedence.AndJoin(branches, aj, [quorum]))
    aj.replies_to(se)
    return m


def convert(model):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return LQN2QN(model)


def the_join(qn):
    for n in qn.getNodes():
        if n.__class__.__name__ == 'Join':
            return n
    return None


def act_think_station(qn):
    for n in qn.getNodes():
        if n.__class__.__name__ == 'Delay' and n.getName() == 'ActivityThink':
            return n
    return None


def test_quorum_join_is_carried_to_the_join_node():
    """A genuine quorum k < n reaches the Join node as a PARTIAL join requiring k."""
    qn = convert(fork_join_model(2, False))
    join = the_join(qn)
    assert join is not None, "the AND-fork must produce a Join node"
    with_quorum = [c for c in qn.getClasses() if join.get_required(c) != -1]
    assert len(with_quorum) == 1, "the quorum belongs to the class that enters the fork"
    assert join.get_required(with_quorum[0]) == 2
    assert join.get_strategy(with_quorum[0]) == JoinStrategy.PARTIAL


def test_full_quorum_leaves_the_join_standard():
    """A quorum equal to the branch count is an ordinary wait-for-all join."""
    qn = convert(fork_join_model(3, False))
    join = the_join(qn)
    for c in qn.getClasses():
        assert join.get_required(c) == -1
        assert join.get_strategy(c) == JoinStrategy.STD


def test_activity_think_time_becomes_a_delay_step():
    """The think time sits on a shared INF station, so the processor is released."""
    qn = convert(fork_join_model(2, True))
    think = act_think_station(qn)
    assert think is not None, "an activity think time must produce the ActivityThink delay"
    served = [c for c in qn.getClasses()
              if think.getServiceProcess(c) is not None
              and think.getServiceProcess(c).getMean() > 1e-8]
    assert len(served) == 1, "exactly the one think-bearing activity has a think step"
    assert think.getServiceProcess(served[0]).getMean() == pytest.approx(0.25)


def test_no_think_time_creates_no_delay_station():
    """Without any activity think time the shared delay is not created at all."""
    assert act_think_station(convert(fork_join_model(2, False))) is None


def test_api_io_lqn2qn_is_the_canonical_converter():
    """``api.io.lqn2qn`` must delegate to the full converter, not a reduced one."""
    from line_solver.api.io import lqn2qn
    qn = convert(fork_join_model(2, True))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        qn2 = lqn2qn(fork_join_model(2, True))
    assert type(qn2) is type(qn)
    assert len(qn2.getNodes()) == len(qn.getNodes())
