#!/usr/bin/env python3
"""Gallery Example: Hyp/Hyp/1-Linear Queue"""

import numpy as np
from line_solver import *

def gallery_hyphyp1_linear(n=2, Umax=0.9):
    model = Network('Hyp/Hyp/1-Linear')

    # Block 1: nodes
    line = [Source(model, 'mySource')]
    for i in range(n):
        line.append(Queue(model, f'Queue{i+1}', SchedStrategy.FCFS))
    line.append(Sink(model, 'mySink'))

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    line[0].setArrival(oclass, HyperExp.fit_mean_and_scv(1, 2))

    means = list(np.linspace(0.1, Umax, n//2))
    if n % 2 == 0:
        means = means + means[::-1]
    else:
        means = means + [Umax] + means[::-1]

    for i in range(n):
        line[i+1].setService(oclass, HyperExp.fit_mean_and_scv(means[i], 1+i+1))

    # Block 3: topology
    model.link(Network.serial_routing(line))

    return model

if __name__ == '__main__':
    model = gallery_hyphyp1_linear()
    print(f"Model: {model.getName()}")
