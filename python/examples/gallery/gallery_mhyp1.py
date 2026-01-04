#!/usr/bin/env python3
"""Gallery Example: M/H/1 Queue"""

from line_solver import *

def gallery_mhyp1():
    model = Network('M/H/1')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, HyperExp.fit_mean_and_scv(0.5, 4))

    # Block 3: topology
    model.link(Network.serial_routing([source, queue, sink]))

    return model

if __name__ == '__main__':
    model = gallery_mhyp1()
    print(f"Model: {model.getName()}")
