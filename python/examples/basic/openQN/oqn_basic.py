"""
Basic Open Queueing Network

This example demonstrates:
- Open network: Source → Delay → Queue → Sink
- Single class with HyperExp and Exp distributions
- Multiple solver comparison
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    node = np.empty(4, dtype=object)
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.FCFS)
    node[2] = Source(model, 'Source')
    node[3] = Sink(model, 'Sink')

    jobclass = OpenClass(model, 'Class1', 0)

    node[0].set_service(jobclass, HyperExp(0.5, 3.0, 10.0))
    node[1].set_service(jobclass, Exp(1))
    node[2].set_arrival(jobclass, Exp(0.1))

    P = model.init_routing_matrix()
    pmatrix = [[0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 0, 0]]
    for i in range(len(node)):
        for j in range(len(node)):
            P.set(jobclass, jobclass, node[i], node[j], pmatrix[i][j])
    model.link(P)

    # Run multiple solvers
    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(model, cutoff=10))
    solver = np.append(solver, FLD(model))
    solver = np.append(solver, MVA(model))
    solver = np.append(solver, MAM(model))
    solver = np.append(solver, NC(model))
    solver = np.append(solver, JMT(model, seed=23000))
    solver = np.append(solver, SSA(model, seed=23000, samples=5000))
    solver = np.append(solver, DES(model, seed=23000, samples=5000))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_table[s] = solver[s].avg_table()
