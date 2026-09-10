#!/usr/bin/env python3
"""Gallery Example: Single-class CQN"""

import random
from line_solver import *

def gallery_cqn(M=2, use_delay=False, seed=23000):
    model = Network('Single-class CQN')

    # Block 1: nodes
    random.seed(seed)
    station = []
    for i in range(M):
        station.append(Queue(model, f'Queue{i+1}', SchedStrategy.PS))

    if use_delay:
        station.append(Delay(model, 'Delay1'))

    # Block 2: classes
    jobclass = ClosedClass(model, 'Class1', round(random.random() * M + 2), station[0], 0)

    for i in range(M):
        station[i].setService(jobclass, Exp.fit_mean(random.random() + i + 1))

    if use_delay:
        station[-1].setService(jobclass, Exp.fit_mean(2.0))

    # Block 3: topology
    P = model.init_routing_matrix()
    P[jobclass, jobclass] = Network.serial_routing(station)
    model.link(P)

    return model

if __name__ == '__main__':
    model = gallery_cqn()
    print(f"Model: {model.getName()}")
