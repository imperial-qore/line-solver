"""Product-form state-dependent routing, Krzesinski (1987), Perform. Eval. 7:125-143.

Section 2.5's central server with two peripheral centers and two levels of
subnetwork nesting, whose routing probabilities are Table 1 of the paper, plus
Section 5.1.1's four-peripheral network comparing SDR against the
state-independent routing of eq. (15).
"""

import numpy as np

from line_solver import (ClosedClass, Exp, Network, Queue, RoutingStrategy,
                         SchedStrategy, SolverCTMC, SolverNC)


def krzesinski_sec25():
    """Central server, branches {2} and {3}, C = (-1,-1), d_12=1, d_13=2, d_23=3."""
    mu = [1.0, 0.8, 0.5]
    N = 3
    model = Network('sdr_central_server')
    node = [Queue(model, 'CPU', SchedStrategy.FCFS),
            Queue(model, 'Disk1', SchedStrategy.FCFS),
            Queue(model, 'Disk2', SchedStrategy.FCFS)]
    jobclass = ClosedClass(model, 'Class1', N, node[0], 0)
    for i in range(3):
        node[i].set_service(jobclass, Exp(mu[i]))
    model.add_link(node[0], node[0])   # denied entry: the busy form of waiting
    model.add_link(node[0], node[1])
    model.add_link(node[0], node[2])
    model.add_link(node[1], node[0])
    model.add_link(node[2], node[0])
    node[1].set_prob_routing(jobclass, node[0], 1.0)
    node[2].set_prob_routing(jobclass, node[0], 1.0)

    d = np.zeros((2, 3))
    d[0, 1] = 1.0
    d[0, 2] = 2.0
    d[1, 2] = 3.0
    node[0].set_state_dep_routing(jobclass, node[0],
                                  [[], [node[1]], [node[2]]],
                                  [0, 1, 2], [-1, -1], d)
    return model


if __name__ == '__main__':
    model = krzesinski_sec25()
    print(SolverNC(model).get_avg_table())
    print(SolverCTMC(model).get_avg_table())
