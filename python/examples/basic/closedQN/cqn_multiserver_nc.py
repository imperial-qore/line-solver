"""
Closed Delay + multiserver FCFS queue.

Exercises the exact multiserver path of NC: with method 'exact' the
c-server station is converted to a load-dependent station with rate min(n,c),
solved via comomld/pfqn_comomrm_ld. CTMC is the exact ground truth; MVA
provides an additional cross-check.
"""

from line_solver import *
import numpy as np


def cqn_multiserver_nc():
    model = Network('model')

    node = np.empty(2, dtype=object)
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.FCFS)
    node[1].set_number_of_servers(3)

    jobclass = ClosedClass(model, 'Class1', 5, node[0], 0)

    node[0].set_service(jobclass, Exp.fit_mean(1.0))  # mean = 1
    node[1].set_service(jobclass, Exp.fit_mean(0.8))  # mean = 0.8

    P = model.init_routing_matrix()
    pmatrix = [[0, 1.0], [1.0, 0]]
    for i in range(len(node)):
        for j in range(len(node)):
            P.set(jobclass, jobclass, node[i], node[j], pmatrix[i][j])
    model.link(P)

    return model


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = cqn_multiserver_nc()

    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(model))
    solver = np.append(solver, MVA(model))
    solver = np.append(solver, NC(model, method='exact'))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name().replace("Solver", "")}')
        avg_table[s] = solver[s].avg_table()
        print(avg_table[s])
