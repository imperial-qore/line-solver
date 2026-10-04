"""
An AND fork/join on a task that ALSO receives an entry-level open arrival.

The Server task has two entries, each with its own fork/join: SE is called by
the closed Client (rendezvous), OE takes an exogenous Poisson stream of rate
0.1. OE carries no reply activity because an open-arrival entry is
send-no-reply, which lqns enforces.

External references on this model: lqns 6.2.28 (valid) gives Client throughput
0.413391, Server task throughput 0.513391, OE throughput 0.1 with open-wait
1.22917 and entry service 0.841667; lqsim (T=5e5, seed 1234) gives 0.4154,
0.52242 and open-wait 1.10546.

SolverLN does NOT solve this combination: the fork-join transform mints its own
Source and collides with the Source the open stream is routed through. The flat
counterpart that IS solved is examples/basic/forkJoin/fj_mixed_openclosed.py.
"""

from line_solver import *


def lqn_fork_open_arrival():
    model = LayeredNetwork('lqnForkOpenArrival')

    P1 = Processor(model, 'P1', float('inf'), SchedStrategy.INF)
    T1 = Task(model, 'Client', 1, SchedStrategy.REF).on(P1)
    T1.set_think_time(Exp.fit_mean(1.0))
    CE = Entry(model, 'CE').on(T1)

    P2 = Processor(model, 'P2', float('inf'), SchedStrategy.INF)
    T2 = Task(model, 'Server', 1, SchedStrategy.FCFS).on(P2)
    T2.set_think_time(Immediate())
    SE = Entry(model, 'SE').on(T2)
    OE = Entry(model, 'OE').on(T2)
    OE.set_arrival(Exp(0.1))

    Activity(model, 'CA', Exp.fit_mean(0.5)).on(T1).bound_to(CE).synch_call(SE)

    RA1 = Activity(model, 'RA1', Exp.fit_mean(0.2)).on(T2).bound_to(SE)
    RA2 = Activity(model, 'RA2', Exp.fit_mean(0.3)).on(T2)
    RA3 = Activity(model, 'RA3', Exp.fit_mean(0.4)).on(T2)
    RA4 = Activity(model, 'RA4', Exp.fit_mean(0.1)).on(T2).replies_to(SE)
    T2.add_precedence(ActivityPrecedence.AndFork(RA1, [RA2, RA3]))
    T2.add_precedence(ActivityPrecedence.AndJoin([RA2, RA3], RA4))

    OA1 = Activity(model, 'OA1', Exp.fit_mean(0.2)).on(T2).bound_to(OE)
    OA2 = Activity(model, 'OA2', Exp.fit_mean(0.3)).on(T2)
    OA3 = Activity(model, 'OA3', Exp.fit_mean(0.4)).on(T2)
    OA4 = Activity(model, 'OA4', Exp.fit_mean(0.1)).on(T2)
    T2.add_precedence(ActivityPrecedence.AndFork(OA1, [OA2, OA3]))
    T2.add_precedence(ActivityPrecedence.AndJoin([OA2, OA3], OA4))

    return model


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.STD)
    model = lqn_fork_open_arrival()

    try:
        print(SolverLN(model, verbose=False).get_avg_table())
    except Exception as err:
        print('LN refuses this model: %s' % err)
