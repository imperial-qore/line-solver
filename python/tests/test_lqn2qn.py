"""
LQN features that LQN2QN must carry into the converted queueing network.

The AND-join quorum, the activity think time and the replication factor used to
be read into the LayeredNetworkStruct and then silently dropped by the
conversion: a k-of-n join became a wait-for-all Join, a think time vanished
instead of appearing as a processor-releasing delay in series with the host
demand, and r replicas collapsed into a single unscaled station.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import warnings

import numpy as np
import pytest

from line_solver import (Activity, ActivityPrecedence, CacheTask, DiscreteSampler, Entry, Erlang, Exp,
                         SetupTask, Immediate, ItemEntry, LayeredNetwork, Processor,
                         ReplacementStrategy, SchedStrategy, Task)
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


# ------------------------------------------------------------------ replication

def replicated_model(rep, fan_out=None, servers=1, task_mult=2):
    """Client -> Server, both replicated rep-fold, one task replica per
    processor replica. An explicit fan-out makes each client replica call every
    server replica instead of its own."""
    m = LayeredNetwork('repl')
    p1 = Processor(m, 'P1', 1, SchedStrategy.INF)
    p2 = Processor(m, 'P2', servers, SchedStrategy.FCFS)
    p1.setReplication(rep)
    p2.setReplication(rep)
    client = Task(m, 'Client', 3, SchedStrategy.REF).on(p1)
    client.setReplication(rep)
    client.set_think_time(Exp(1.0))
    ce = Entry(m, 'CE').on(client)
    server = Task(m, 'Server', task_mult, SchedStrategy.FCFS).on(p2)
    server.setReplication(rep)
    if fan_out is not None:
        client.setFanOut('Server', fan_out)
    se = Entry(m, 'SE').on(server)
    ca = Activity(m, 'ca', Exp(10.0)).on(client).bound_to(ce)
    ca.synch_call(se, 1.0)
    a0 = Activity(m, 'a0', Exp(2.0)).on(server).bound_to(se)
    a0.replies_to(se)
    return m


def convert_repl(model, mode):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return LQN2QN(model, replication=mode)


def station_named(qn, name):
    for n in qn.getNodes():
        if n.getName() == name:
            return n
    return None


def test_materialised_replicas_are_separate_stations_and_chains():
    """Three replicas give three processor stations and three closed chains."""
    qn = convert_repl(replicated_model(3), 'materialize')
    names = [n.getName() for n in qn.getNodes()]
    assert ['P2', 'P2_r2', 'P2_r3'] == [n for n in names if n.startswith('P2')]
    thinks = [c for c in qn.getClasses() if c.getName().startswith('Client_Think')]
    assert len(thinks) == 3, "one closed chain per reference task replica"
    assert all(c.population == 3 for c in thinks), "each replica keeps its own population"


def test_pairing_replicates_the_subsystem_without_coupling_it():
    """At fan-out 1 replica i calls replica i, so the chains stay disjoint."""
    qn = convert_repl(replicated_model(2, fan_out=1), 'materialize')
    sn = qn.getStruct()
    p2, p22 = station_named(qn, 'P2'), station_named(qn, 'P2_r2')
    served = [[c for c in qn.getClasses()
               if st.getServiceProcess(c) is not None
               and st.getServiceProcess(c).getMean() > 1e-8] for st in (p2, p22)]
    assert len(served[0]) == 1 and len(served[1]) == 1
    assert served[0][0] is not served[1][0], "a replica serves only its own class"
    assert float(np.nansum(np.asarray(sn.njobs))) == 6.0


def test_broadcast_fanout_couples_every_caller_replica_to_every_replica():
    """At fan-out r each caller replica reaches all r callee replicas."""
    qn = convert_repl(replicated_model(2, fan_out=2), 'materialize')
    p2, p22 = station_named(qn, 'P2'), station_named(qn, 'P2_r2')
    served = [len([c for c in qn.getClasses()
                   if st.getServiceProcess(c) is not None
                   and st.getServiceProcess(c).getMean() > 1e-8]) for st in (p2, p22)]
    assert served == [2, 2], "both caller replicas reach both server replicas"


def test_pooled_replicas_scale_servers_population_and_admission():
    """Pooling keeps one station of r times the capacity and one r-fold chain."""
    qn = convert_repl(replicated_model(3, servers=2, task_mult=2), 'pool')
    assert station_named(qn, 'P2_r2') is None, "the replicas collapse into one station"
    assert station_named(qn, 'P2').getNumberOfServers() == 6
    thinks = [c for c in qn.getClasses() if c.getName().startswith('Client_Think')]
    assert len(thinks) == 1 and thinks[0].population == 9
    assert float(np.asarray(qn.regions[0]._constraint_b).ravel()[0]) == 6.0


def test_materialised_thread_pools_are_one_admission_row_per_replica():
    """Each task replica owns its own threads, so it owns its own row."""
    qn = convert_repl(replicated_model(3, servers=2, task_mult=2), 'materialize')
    A = np.asarray(qn.regions[0]._constraint_A)
    b = np.asarray(qn.regions[0]._constraint_b).ravel()
    assert A.shape[0] == 3 and list(b) == [2.0, 2.0, 2.0]
    assert (A.sum(axis=0) <= 1).all(), "no class is admitted by two replica rows"


def test_pooling_warns_and_materialisation_does_not():
    """The pooled representation is an approximation and must say so."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        LQN2QN(replicated_model(2), replication='pool')
    assert any('pooled' in str(c.message) for c in caught)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        LQN2QN(replicated_model(2), replication='materialize')
    assert not any('pooled' in str(c.message) for c in caught)


def test_auto_materialises_a_small_model():
    """Within the instantiation budget 'auto' is materialisation."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        qn = LQN2QN(replicated_model(3))
    assert not any('pooled' in str(c.message) for c in caught)
    assert station_named(qn, 'P2_r3') is not None


def test_auto_pools_a_broadcast_expansion_over_budget():
    """A broadcast fan-out down a deep call chain falls back to pooling."""
    m = LayeredNetwork('deep')
    p0 = Processor(m, 'P0', 1, SchedStrategy.INF)
    client = Task(m, 'Client', 2, SchedStrategy.REF).on(p0)
    client.set_think_time(Exp(1.0))
    ce = Entry(m, 'CE').on(client)
    prev_task = client
    prev_act = Activity(m, 'ca', Exp(10.0)).on(client).bound_to(ce)
    for d in range(1, 5):
        pd = Processor(m, 'P%d' % d, 1, SchedStrategy.FCFS)
        pd.setReplication(4)
        t = Task(m, 'T%d' % d, 5, SchedStrategy.FCFS).on(pd)
        t.setReplication(4)
        e = Entry(m, 'E%d' % d).on(t)
        prev_task.setFanOut('T%d' % d, 4)
        prev_act.synch_call(e, 1.0)
        a = Activity(m, 'a%d' % d, Exp(2.0)).on(t).bound_to(e)
        a.replies_to(e)
        prev_task, prev_act = t, a
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        qn = LQN2QN(m)
    assert any('pooled' in str(c.message) for c in caught)
    assert station_named(qn, 'P1_r2') is None


def test_unknown_replication_mode_is_rejected():
    with pytest.raises(ValueError):
        LQN2QN(replicated_model(2), replication='surrogate')


def test_no_replication_leaves_the_conversion_untouched():
    """A model without replication converts identically in all three modes."""
    def shape(mode):
        qn = convert_repl(fork_join_model(2, True), mode)
        return ([n.getName() for n in qn.getNodes()],
                [c.getName() for c in qn.getClasses()])
    assert shape('auto') == shape('materialize') == shape('pool')


# -------------------------------------------------------------------- retrieval

def cache_model(retrieval):
    """Client -> CacheTask, the read forking into a hit and a miss activity, the
    miss being the fetch. With retrieval the fetch becomes a retrieval system."""
    m = LayeredNetwork('lcq')
    p1 = Processor(m, 'P1', 1, SchedStrategy.INF)
    t1 = Task(m, 'T1', 4, SchedStrategy.REF).on(p1)
    t1.set_think_time(Exp(1.0))
    e1 = Entry(m, 'E1').on(t1)
    pc = Processor(m, 'PC', 1, SchedStrategy.PS)
    c2 = CacheTask(m, 'C2', 4, 2, ReplacementStrategy.RR, 1).on(pc)
    if retrieval:
        c2.set_retrieval(True)
    i2 = ItemEntry(m, 'I2', 4, DiscreteSampler([0.25] * 4)).on(c2)
    Activity(m, 'A1', Immediate()).on(t1).bound_to(e1).synch_call(i2, 1)
    ac2 = Activity(m, 'AC2', Immediate()).on(c2).bound_to(i2)
    hit = Activity(m, 'AC2h', Exp(1.0)).on(c2).replies_to(i2)
    miss = Activity(m, 'AC2m', Exp(0.5)).on(c2).replies_to(i2)
    c2.add_precedence(ActivityPrecedence.cache_access(ac2, [hit, miss]))
    return m


def test_retrieval_creates_a_fetch_station():
    """A delayed-hit CacheTask gets one PS fetch station holding the miss demand."""
    qn = convert(cache_model(True))
    fetch = station_named(qn, 'C2_Cache_Fetch')
    assert fetch is not None, "the retrieval system must be a station of the QN"
    read = [c for c in qn.getClasses() if c.getName() == 'AC2'][0]
    assert fetch.getServiceProcess(read).getMean() == pytest.approx(2.0), \
        "the fetch time is the miss activity's demand"


def test_retrieval_moves_the_miss_demand_off_the_processor():
    """The fetch happens in the retrieval system, so it is not charged twice."""
    qn = convert(cache_model(True))
    pc = station_named(qn, 'PC')
    miss = [c for c in qn.getClasses() if c.getName() == 'AC2m'][0]
    hit = [c for c in qn.getClasses() if c.getName() == 'AC2h'][0]
    assert pc.getServiceProcess(hit).getMean() == pytest.approx(1.0)
    svc = pc.getServiceProcess(miss)
    assert svc is None or svc.getMean() < 1e-8, "the miss demand moved to the fetch station"


def test_no_retrieval_creates_no_fetch_station():
    """Without the retrieval flag the miss branch keeps its demand at the processor."""
    qn = convert(cache_model(False))
    assert station_named(qn, 'C2_Cache_Fetch') is None
    pc = station_named(qn, 'PC')
    miss = [c for c in qn.getClasses() if c.getName() == 'AC2m'][0]
    assert pc.getServiceProcess(miss).getMean() == pytest.approx(2.0)


# ---------------------------------------------------------------- setup task

def function_model(setup, host_sched=SchedStrategy.FCFS):
    """Client -> SetupTask: with a setup and a delay-off time the function's
    server shuts down when idle and pays a cold start on the next arrival."""
    m = LayeredNetwork('faas')
    p1 = Processor(m, 'P1', 1, SchedStrategy.INF)
    pf = Processor(m, 'PF', 1, host_sched)
    c = Task(m, 'Client', 2, SchedStrategy.REF).on(p1)
    c.set_think_time(Exp(0.5))
    ce = Entry(m, 'CE').on(c)
    f = SetupTask(m, 'F', 1, SchedStrategy.FCFS).on(pf)
    if setup:
        f.set_setup_time(Exp(2.0))
        f.set_delay_off_time(Exp(1.0))
    fe = Entry(m, 'FE').on(f)
    Activity(m, 'ca', Exp(10.0)).on(c).bound_to(ce).synch_call(fe, 1.0)
    Activity(m, 'a0', Exp(2.0)).on(f).bound_to(fe).replies_to(fe)
    return m


def armed_classes(qn, station_name):
    st = station_named(qn, station_name)
    return [(c.getName(), st.getSetupTime(c).getMean(), st.getDelayOffTime(c).getMean())
            for c in qn.getClasses() if st.getSetupTime(c) is not None]


def test_function_task_setup_reaches_the_host_station():
    """The setup/delay-off pair is set per step class of the setup task."""
    qn = convert(function_model(True))
    assert armed_classes(qn, 'PF') == [('a0', pytest.approx(0.5), pytest.approx(1.0))]


def test_no_setup_leaves_the_station_always_on():
    """An ordinary task never arms the setup/delay-off pair."""
    assert armed_classes(convert(function_model(False)), 'PF') == []


def test_setup_on_an_infinite_server_is_warned_and_dropped():
    """An infinite-server processor never shuts down, so it never sets up."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        qn = LQN2QN(function_model(True, host_sched=SchedStrategy.INF))
    assert any('never shuts down' in str(c.message) for c in caught)
    assert armed_classes(qn, 'PF') == []


def moment_model(think, call_mean):
    """Reference task with the given think time, calling a server call_mean times."""
    m = LayeredNetwork('mom')
    p1 = Processor(m, 'P1', 1, SchedStrategy.INF)
    p2 = Processor(m, 'P2', 1, SchedStrategy.PS)
    c = Task(m, 'C', 2, SchedStrategy.REF).on(p1)
    c.set_think_time(think)
    ce = Entry(m, 'CE').on(c)
    s = Task(m, 'S', 1, SchedStrategy.FCFS).on(p2)
    se = Entry(m, 'SE').on(s)
    Activity(m, 'ca', Exp(1.0)).on(c).bound_to(ce).synch_call(se, call_mean)
    Activity(m, 'a0', Exp(2.0)).on(s).bound_to(se).replies_to(se)
    return m


def test_think_time_keeps_its_scv_through_the_conversion():
    """The think delay carries the task distribution, not its mean fitted to an Exp."""
    qn = convert(moment_model(Erlang.fitMeanAndSCV(2.0, 0.25), 1.0))
    delay = station_named(qn, 'C_Think')
    svc = [delay.getService(c) for c in qn.getClasses()
           if delay.getService(c) is not None and delay.getService(c).getMean() > 1e-9]
    assert len(svc) == 1
    assert svc[0].getMean() == pytest.approx(2.0)
    assert svc[0].getSCV() == pytest.approx(0.25)


def test_call_count_below_one_is_bernoulli():
    """A call that happens with probability p has mean p and SCV (1-p)/p, which
    Geometric(1/p) cannot represent: its parameter would exceed 1."""
    lsn = moment_model(Exp(1.0), 0.4).getStruct()
    dist = lsn.callproc[0]
    assert dist.__class__.__name__ == 'Bernoulli'
    assert dist.getMean() == pytest.approx(0.4)
    assert dist.getSCV() == pytest.approx(0.6 / 0.4)


def test_call_count_above_one_is_geometric():
    """At or above one call the count is geometric with the declared mean."""
    dist = moment_model(Exp(1.0), 3.0).getStruct().callproc[0]
    assert dist.__class__.__name__ == 'Geometric'
    assert dist.getMean() == pytest.approx(3.0)
    assert dist.getSCV() == pytest.approx(1.0 - 1.0 / 3.0)


def test_call_count_of_zero_is_a_degenerate_placeholder():
    """A call declared with mean 0 is not a call: the count is the zero-mean
    placeholder, the same object MATLAB and the JAR store."""
    dist = moment_model(Exp(1.0), 0.0).getStruct().callproc[0]
    assert dist.__class__.__name__ == 'Immediate'
    assert dist.getMean() == pytest.approx(0.0)


def test_call_of_zero_leaves_the_callee_unreachable_and_is_rejected():
    """A zero-mean call reaches no step on the callee's processor, so that Queue
    would carry no class at all. MATLAB and the JAR both reject the converted
    model in sanitize; the native Python link now does too."""
    with pytest.raises(ValueError, match='no service configured'):
        convert(moment_model(Exp(1.0), 0.0))
