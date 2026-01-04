"""
Basic Mixed Queueing Network

This example demonstrates:
- Mixed network with 1 closed class and 1 open class
- 2 closed jobs in the closed class
- PS scheduling strategy
- Erlang and HyperExp distributions
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    node = np.empty(4, dtype=object)
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.PS)
    node[2] = Source(model, 'Source')
    node[3] = Sink(model, 'Sink')

    jobclass = np.empty(2, dtype=object)
    jobclass[0] = ClosedClass(model, 'ClosedClass', 2, node[0], 0)
    jobclass[1] = OpenClass(model, 'OpenClass', 0)

    node[0].set_service(jobclass[0], Erlang(3, 2))
    node[0].set_service(jobclass[1], HyperExp(0.5, 3.0, 10.0))

    node[1].set_service(jobclass[0], HyperExp(0.1, 1.0, 10.0))
    node[1].set_service(jobclass[1], Exp(1))

    node[2].set_arrival(jobclass[1], Exp(0.1))

    M = model.get_number_of_nodes()
    K = model.get_number_of_classes()

    P = model.init_routing_matrix()

    # Closed class routing (circular between delay and queue)
    P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
    P.set(jobclass[0], jobclass[0], node[1], node[0], 1.0)

    # Open class routing (source -> delay -> queue -> sink)
    P.set(jobclass[1], jobclass[1], node[2], node[0], 1.0)
    P.set(jobclass[1], jobclass[1], node[0], node[1], 1.0)
    P.set(jobclass[1], jobclass[1], node[1], node[3], 1.0)

    model.link(P)

    # Run multiple solvers
    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(model, keep=True, verbose=True, cutoff=3, seed=23000))
    solver = np.append(solver, JMT(model, keep=True, verbose=True, seed=23000))
    solver = np.append(solver, SSA(model, keep=True, verbose=True, seed=23000))
    solver = np.append(solver, MVA(model, keep=True, verbose=True))
    solver = np.append(solver, DES(model, keep=True, verbose=True, seed=23000))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_table[s] = solver[s].avg_table()
