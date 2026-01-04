"""
Asymmetric Fork-Join with Different Service Rates

This example demonstrates:
- Open fork-join network
- 3 parallel queues with different service rates
- Demonstrates impact of asymmetry on response time
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    source = Source(model, 'Source')
    queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
    queue3 = Queue(model, 'Queue3', SchedStrategy.PS)
    fork = Fork(model, 'Fork')
    join = Join(model, 'Join', fork)
    sink = Sink(model, 'Sink')

    jobclass1 = OpenClass(model, 'class1')

    source.set_arrival(jobclass1, Exp(0.2))
    # Different service rates create asymmetry
    queue1.set_service(jobclass1, Exp(1.0))
    queue2.set_service(jobclass1, Exp(2.0))
    queue3.set_service(jobclass1, Exp(3.0))

    P = model.init_routing_matrix()
    P.set(jobclass1, jobclass1, source, fork, 1.0)
    P.set(jobclass1, jobclass1, fork, queue1, 1.0)
    P.set(jobclass1, jobclass1, fork, queue2, 1.0)
    P.set(jobclass1, jobclass1, fork, queue3, 1.0)
    P.set(jobclass1, jobclass1, queue1, join, 1.0)
    P.set(jobclass1, jobclass1, queue2, join, 1.0)
    P.set(jobclass1, jobclass1, queue3, join, 1.0)
    P.set(jobclass1, jobclass1, join, sink, 1.0)

    model.link(P)

    solver = np.array([], dtype=object)
    solver = np.append(solver, JMT(model, seed=23000))
    solver = np.append(solver, MVA(model))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_table[s] = solver[s].avg_table()
