"""
Mixed Network with Single-Server PS Queues

This example demonstrates:
- Mixed network with 100 closed jobs and open arrivals
- 4 PS (Processor Sharing) queues with single servers
- Closed class visits all 4 queues, open class visits first 3
- Exponential arrivals for open class
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    M = 4
    node = np.empty(M, dtype=object)
    node[0] = Queue(model, 'Queue1', SchedStrategy.PS)
    node[1] = Queue(model, 'Queue2', SchedStrategy.PS)
    node[2] = Queue(model, 'Queue3', SchedStrategy.PS)
    node[3] = Queue(model, 'Queue4', SchedStrategy.PS)  # only closed class

    source = Source(model, 'Source')
    sink = Sink(model, 'Sink')

    jobclass = np.empty(2, dtype=object)
    jobclass[0] = ClosedClass(model, 'ClosedClass', 100, node[0], 0)
    jobclass[1] = OpenClass(model, 'OpenClass', 0)

    for i in range(M):
        node[i].set_service(jobclass[0], Exp(i + 1))
        node[i].set_service(jobclass[1], Exp((i + 1) ** 0.5))

    source.set_arrival(jobclass[1], Exp(0.3))

    # Routing
    P = model.init_routing_matrix()

    # Closed class: serial routing through all 4 queues
    P.set(jobclass[0], jobclass[0], Network.serial_routing(node[0], node[1], node[2], node[3]))

    # Open class: source -> queue1 -> queue2 -> queue3 -> sink
    P.set(jobclass[1], jobclass[1], Network.serial_routing(source, node[0], node[1], node[2], sink))

    model.link(P)

    # Run multiple solvers
    solver = np.array([], dtype=object)
    solver = np.append(solver, JMT(model, keep=False, verbose=True, cutoff=3, seed=23000, samples=100000))
    solver = np.append(solver, SSA(model, keep=False, verbose=True, seed=23000, samples=100000))
    solver = np.append(solver, FLD(model, keep=False, verbose=True, cutoff=3, seed=23000))
    solver = np.append(solver, MVA(model, method='lin'))
    solver = np.append(solver, NC(model, keep=False, verbose=True, cutoff=3, seed=23000))
    solver = np.append(solver, DES(model, keep=False, verbose=True, cutoff=3, seed=23000, samples=100000))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_table[s] = solver[s].avg_table()
        print(avg_table[s])
