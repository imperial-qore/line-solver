"""
Closed network with Blocking-After-Service (BAS) finite-buffer blocking.

Queue1 uses the BAS drop rule; Queue2 has a finite buffer (capacity 1).
The 'sqd' (Smith Queue Decomposition) approximation handles this model;
MVA also auto-routes BAS models to 'sqd' under its default method.
"""

from line_solver import *
import numpy as np


def cqn_bas_blocking():
    model = Network('cqn_bas_blocking')

    node = np.empty(2, dtype=object)
    node[0] = Queue(model, 'Queue1', SchedStrategy.FCFS)
    node[1] = Queue(model, 'Queue2', SchedStrategy.FCFS)

    jobclass = ClosedClass(model, 'Class1', 2, node[0], 0)

    node[0].set_service(jobclass, Exp(1.0))
    node[1].set_service(jobclass, Exp(0.8))

    node[1].setCap(1)
    node[0].setDropRule(jobclass, DropStrategy.BAS)

    model.link(Network.serialRouting(node[0], node[1]))
    return model


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = cqn_bas_blocking()

    print('\nSOLVER: MVA (method=sqd)')
    solver = MVA(model, 'sqd')
    avg_table = solver.getAvgTable()
    print(avg_table)
