"""
Two Queues Multi-class Closed Queueing Network

This example demonstrates:
- Closed network with 2 classes (10 jobs each)
- Delay station and 2 FCFS queues
- Serial routing for both classes
"""

from line_solver import *
import numpy as np


def cqn_twoqueues_multi():

    model = Network('model')

    # Nodes
    node = np.empty(3, dtype=object)
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.FCFS)
    node[2] = Queue(model, 'Queue2', SchedStrategy.FCFS)

    # Job classes
    jobclass = np.empty(2, dtype=object)
    jobclass[0] = ClosedClass(model, 'Class1', 10, node[0], 0)
    jobclass[1] = ClosedClass(model, 'Class2', 10, node[0], 0)

    # Service times
    node[0].set_service(jobclass[0], Exp.fit_mean(1.0))  # mean = 1
    node[1].set_service(jobclass[0], Exp.fit_mean(1.5))  # mean = 1.5
    node[2].set_service(jobclass[0], Exp.fit_mean(3.0))  # mean = 3.0

    node[0].set_service(jobclass[1], Exp.fit_mean(1.0))  # mean = 1
    node[1].set_service(jobclass[1], Exp.fit_mean(1.5))  # mean = 1.5
    node[2].set_service(jobclass[1], Exp.fit_mean(3.0))  # mean = 3.0

    # Routing: serial routing for both classes
    P = model.init_routing_matrix()
    P[jobclass[0], jobclass[0]] = Network.serial_routing(node)
    P[jobclass[1], jobclass[1]] = Network.serial_routing(node)
    model.link(P)

    return model


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.STD)
    model = cqn_twoqueues_multi()

    # Solve with MVA
    solver = MVA(model)
    print(f'\nSOLVER: {solver.get_name()}')
    avg_table = solver.avg_table()
    print(avg_table)
