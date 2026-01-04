#!/usr/bin/env python3
"""Gallery Example: Finite Repairmen CQN"""

import random
from line_solver import *

def gallery_repairmen(nservers=1, seed=2300):
    model = Network('Finite repairmen CQN')

    # Block 1: nodes
    M = 1
    random.seed(seed)

    node = []
    node.append(Queue(model, 'Queue1', SchedStrategy.PS))
    node[0].setNumberOfServers(nservers)
    node.append(Delay(model, 'Delay1'))

    # Block 2: classes
    jobclass = ClosedClass(model, 'Class1', round(random.random() * 10 * M + 3), node[0], 0)

    node[0].setService(jobclass, Exp.fit_mean(random.random() + 1))
    node[1].setService(jobclass, Exp.fit_mean(2.0))

    # Block 3: topology
    P = model.init_routing_matrix()
    P[jobclass, jobclass] = Network.serial_routing(node)
    model.link(P)

    return model

if __name__ == '__main__':
    model = gallery_repairmen()
    print(f"Model: {model.getName()}")
