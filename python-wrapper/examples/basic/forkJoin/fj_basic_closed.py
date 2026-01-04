"""
Basic Closed Fork-Join Network

This example demonstrates:
- Closed network with 5 jobs
- Delay → Fork → Parallel Queues → Join → Delay
- Fork-join with PS scheduling at queues
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    delay = Delay(model, 'Delay')
    queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
    fork = Fork(model, 'Fork')
    join = Join(model, 'Join', fork)

    jobclass1 = ClosedClass(model, 'class1', 5, delay)

    delay.set_service(jobclass1, Exp(1.0))
    queue1.set_service(jobclass1, Exp(1.0))
    queue2.set_service(jobclass1, Exp(1.0))

    P = model.init_routing_matrix()
    P.set(jobclass1, jobclass1, delay, fork, 1.0)
    P.set(jobclass1, jobclass1, fork, queue1, 1.0)
    P.set(jobclass1, jobclass1, fork, queue2, 1.0)
    P.set(jobclass1, jobclass1, queue1, join, 1.0)
    P.set(jobclass1, jobclass1, queue2, join, 1.0)
    P.set(jobclass1, jobclass1, join, delay, 1.0)

    model.link(P)

    solver = np.array([], dtype=object)
    solver = np.append(solver, JMT(model, seed=23000))
    # Use AMVA method for better convergence with fork-join
    solver = np.append(solver, MVA(model, method='amva'))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_table[s] = solver[s].avg_table()
