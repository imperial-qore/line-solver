"""
Basic Open Queueing Network

This example demonstrates:
- Open network: Source -> Delay -> Queue -> Sink
- Single class with HyperExp and Exp distributions
- Multiple solver comparison
"""

from line_solver import *
import numpy as np


def oqn_basic():
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

    return model


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    # The model is built by the function above, not inline: an example that
    # builds it only under __main__ exposes nothing on import, so the JAVA and
    # C++ rows of parity-static cannot export it and SKIP every solver. That
    # reads as coverage while asserting nothing (see _example_model_vendor.py).
    model = oqn_basic()

    # Run multiple solvers
    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(model, cutoff=10))
    solver = np.append(solver, FLD(model))
    solver = np.append(solver, MVA(model))
    solver = np.append(solver, MAM(model))
    solver = np.append(solver, NC(model))
    solver = np.append(solver, JMT(model, seed=23000))
    # SSA needs a larger budget than the other solvers here: at 5e3 events the
    # standard error of Util/Tput on this rho=0.1 queue is ~4%, the same size as
    # the cross-codebase parity tolerance. 5e4 brings it to ~1.4%.
    solver = np.append(solver, SSA(model, seed=23000, samples=50000))
    # 1e4 was inherited from the MATLAB row, which used to have its LDES sample
    # budget clobbered by a generic options struct; both now use the SolverLDES
    # default of 2e5.
    solver = np.append(solver, LDES(model, seed=23000, samples=200000))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_table[s] = solver[s].avg_table()
        print(avg_table[s])
