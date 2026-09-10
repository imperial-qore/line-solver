#!/usr/bin/env python3
"""Round-robin call dispatch over a set of target tasks.

A client task issues its synchronous calls to three interchangeable server
tasks in cyclic order rather than by probabilistic branching. The two models
carry the same aggregate call rate; what differs is that round-robin removes
the variance of the branching, which smooths the server queues.

Only the squashed ('flat') layering can express this: under 'srvn' each server
task lives in its own submodel and is replaced, in the client's submodel, by a
surrogate delay, so no node ever has arcs to more than one of them.
"""

from line_solver import *
from line_solver.solvers import SolverLNOptions


def lqn_rrobin():
    model = LayeredNetwork('LQN-RRobin')

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

    # ONE call per invocation, its destination cycling over the three servers.
    # mean_calls=1 over 3 targets is where round-robin actually bites: the
    # probabilistic twin makes 0..3 calls per invocation with the same mean,
    # round-robin makes exactly one. At mean_calls=3 the two agree on
    # everything except the order.
    Activity(model, 'AC', Exp(2)).on(TC).bound_to(EC).synch_call_rrobin([ES1, ES2, ES3], 1)
    Activity(model, 'AS1', Exp(1)).on(TS1).bound_to(ES1).replies_to(ES1)
    Activity(model, 'AS2', Exp(1)).on(TS2).bound_to(ES2).replies_to(ES2)
    Activity(model, 'AS3', Exp(1)).on(TS3).bound_to(ES3).replies_to(ES3)

    return model


def lqn_rrobin_probabilistic():
    """The ungrouped twin: same per-target call means, branching at random."""
    model = LayeredNetwork('LQN-RRobin-Prob')

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

    AC = Activity(model, 'AC', Exp(2)).on(TC).bound_to(EC)
    AC.synch_call(ES1, 1 / 3).synch_call(ES2, 1 / 3).synch_call(ES3, 1 / 3)
    Activity(model, 'AS1', Exp(1)).on(TS1).bound_to(ES1).replies_to(ES1)
    Activity(model, 'AS2', Exp(1)).on(TS2).bound_to(ES2).replies_to(ES2)
    Activity(model, 'AS3', Exp(1)).on(TS3).bound_to(ES3).replies_to(ES3)

    return model


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    options = SolverLNOptions()
    options.config['layering'] = 'flat'
    options.verbose = False
    # SolverMVA is rejected here: a routed call group needs a layer solver with
    # state-dependent routing, and MVA would silently return the probabilistic
    # split under a round-robin label.
    print(SolverLN(lqn_rrobin(), SolverSSA, options).get_avg_table())
