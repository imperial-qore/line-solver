#!/usr/bin/env python3
"""Join-the-shortest-queue call dispatch over a set of target tasks.

The twin of lqn_rrobin: a client task issues its synchronous calls to three
interchangeable server tasks, but the target is the one holding the fewest
jobs at dispatch time rather than the next one in cyclic order. The two models
carry the same aggregate call rate; what JSQ adds over round-robin is that the
dispatch reacts to the state of the servers, so it also absorbs asymmetry in
the service times, not only the variance of the branching.

Only the squashed ('flat') layering can express this: under 'srvn' each server
task lives in its own submodel and is replaced, in the client's submodel, by a
surrogate delay, so no node ever has arcs to more than one of them. The layer
solver must also implement state-dependent routing (SSA here); MVA, NC and FLD
are rejected rather than silently returning the probabilistic split.
"""

from line_solver import *
from line_solver.solvers import SolverLNOptions


def lqn_jsq():
    model = LayeredNetwork('LQN-JSQ')

    PC = Processor(model, 'PC', 1, SchedStrategy.INF)
    PS = Processor(model, 'PS', 1, SchedStrategy.PS)

    TC = Task(model, 'TC', 10, SchedStrategy.REF).on(PC).set_think_time(Exp(1 / 5))
    TS1 = Task(model, 'TS1', 5, SchedStrategy.FCFS).on(PS)
    TS2 = Task(model, 'TS2', 5, SchedStrategy.FCFS).on(PS)
    TS3 = Task(model, 'TS3', 5, SchedStrategy.FCFS).on(PS)

    EC = Entry(model, 'EC').on(TC)
    ES1 = Entry(model, 'ES1').on(TS1)
    ES2 = Entry(model, 'ES2').on(TS2)
    ES3 = Entry(model, 'ES3').on(TS3)

    # ONE call per invocation, its destination the least loaded of the three
    # servers. As in the round-robin twin, mean_calls=1 over 3 targets is where
    # the strategy bites: the probabilistic model makes 0..3 calls per
    # invocation with the same mean, JSQ makes exactly one, to the shortest
    # queue.
    Activity(model, 'AC', Exp(2)).on(TC).bound_to(EC).synch_call_jsq([ES1, ES2, ES3], 1)
    Activity(model, 'AS1', Exp(1)).on(TS1).bound_to(ES1).replies_to(ES1)
    Activity(model, 'AS2', Exp(1)).on(TS2).bound_to(ES2).replies_to(ES2)
    Activity(model, 'AS3', Exp(1)).on(TS3).bound_to(ES3).replies_to(ES3)

    return model


def lqn_jsq_asymmetric():
    """Servers of unequal speed, where JSQ and round-robin part company.

    Round-robin hands each server exactly one third of the calls whatever their
    speed; JSQ follows the queues, so the fast server absorbs more of them and
    the entry throughputs come out unequal.
    """
    model = LayeredNetwork('LQN-JSQ-Asym')

    PC = Processor(model, 'PC', 1, SchedStrategy.INF)
    PS = Processor(model, 'PS', 10, SchedStrategy.PS)

    TC = Task(model, 'TC', 4, SchedStrategy.REF).on(PC).set_think_time(Exp(1 / 2))
    TS1 = Task(model, 'TS1', 4, SchedStrategy.FCFS).on(PS)
    TS2 = Task(model, 'TS2', 4, SchedStrategy.FCFS).on(PS)

    EC = Entry(model, 'EC').on(TC)
    ES1 = Entry(model, 'ES1').on(TS1)
    ES2 = Entry(model, 'ES2').on(TS2)

    Activity(model, 'AC', Exp(10)).on(TC).bound_to(EC).synch_call_jsq([ES1, ES2], 1)
    Activity(model, 'AS1', Exp(1)).on(TS1).bound_to(ES1).replies_to(ES1)
    Activity(model, 'AS2', Exp(4)).on(TS2).bound_to(ES2).replies_to(ES2)

    return model


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    options = SolverLNOptions()
    options.config['layering'] = 'flat'
    options.verbose = False
    print(SolverLN(lqn_jsq(), SolverSSA, options).get_avg_table())
