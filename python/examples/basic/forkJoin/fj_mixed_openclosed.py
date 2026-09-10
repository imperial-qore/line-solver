"""
Fork-join traversed by an open AND a closed class at once.

Flat counterpart of examples/basic/layeredModel/lqn_fork_open_arrival.py: the
same fork carries an exogenous Poisson stream (rate 0.1) and a closed chain
(1 job, think time 1). The MVA fork-join transform solves this; the layered
builder refuses the analogous LQN, because the transform's auxiliary Source
collides with the layer's open-stream Source.

Reference (SolverJMT, seed 23000, 2e5 samples): Branch1 open RespT 0.362,
Branch1 closed RespT 0.310, Join open RespT 0.215, closed throughput 0.540.
"""

from line_solver import *


def fj_mixed_openclosed():
    model = Network('ForkJoinOpenClosed')

    source = Source(model, 'Source')
    client = Delay(model, 'Client')
    prefork = Queue(model, 'PreFork', SchedStrategy.FCFS)
    branch1 = Queue(model, 'Branch1', SchedStrategy.FCFS)
    branch2 = Queue(model, 'Branch2', SchedStrategy.FCFS)
    postjoin = Queue(model, 'PostJoin', SchedStrategy.FCFS)
    fork = Fork(model, 'Fork')
    join = Join(model, 'Join', fork)
    sink = Sink(model, 'Sink')

    oclass = OpenClass(model, 'Open')
    cclass = ClosedClass(model, 'Closed', 1, client)

    source.set_arrival(oclass, Exp(0.1))
    client.set_service(cclass, Exp.fit_mean(1.0))
    client.set_service(oclass, Disabled())
    for station, mean in ((prefork, 0.2), (branch1, 0.3), (branch2, 0.4), (postjoin, 0.1)):
        station.set_service(oclass, Exp.fit_mean(mean))
        station.set_service(cclass, Exp.fit_mean(mean))

    P = model.init_routing_matrix()
    for jobclass, entry, exit_node in ((oclass, source, sink), (cclass, client, client)):
        P.set(jobclass, jobclass, entry, prefork, 1.0)
        P.set(jobclass, jobclass, prefork, fork, 1.0)
        P.set(jobclass, jobclass, fork, branch1, 1.0)
        P.set(jobclass, jobclass, fork, branch2, 1.0)
        P.set(jobclass, jobclass, branch1, join, 1.0)
        P.set(jobclass, jobclass, branch2, join, 1.0)
        P.set(jobclass, jobclass, join, postjoin, 1.0)
        P.set(jobclass, jobclass, postjoin, exit_node, 1.0)

    model.link(P)
    return model


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.STD)
    model = fj_mixed_openclosed()

    print('\nMVA Results:')
    print(SolverMVA(model).get_avg_table())
    print('\nJMT Results:')
    print(SolverJMT(model, seed=23000, samples=200000).get_avg_table())
