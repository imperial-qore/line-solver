#!/usr/bin/env python3
"""Gallery Example: Multi-class CQN"""

import random
from line_solver import *

def gallery_cqn_multiclass(m=1, r=2, wantdelay=True):
    model = Network('Multi-class CQN')

    # Block 1: nodes
    node = []
    for i in range(m):
        node.append(Queue(model, f'Queue {i+1}', SchedStrategy.PS))

    if wantdelay:
        node.append(Delay(model, 'Delay 1'))

    # Block 2: classes
    jobclass = []
    for s in range(r):
        jobclass.append(ClosedClass(model, f'Class{s+1}', 5, node[0], 0))

    for s in range(r):
        for i in range(m):
            node[i].setService(jobclass[s], Exp.fit_mean(round(50 * random.random())))
        if wantdelay:
            node[-1].setService(jobclass[s], Exp.fit_mean(round(100 * random.random())))

    # Block 3: topology
    P = model.init_routing_matrix()
    for s in range(r):
        P[jobclass[s], jobclass[s]] = Network.serial_routing(node)
    model.link(P)

    return model

if __name__ == '__main__':
    model = gallery_cqn_multiclass()
    print(f"Model: {model.getName()}")
