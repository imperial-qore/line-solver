#!/usr/bin/env python3
"""Gallery Example: Cox/M/1 Queue"""

from line_solver import *

def gallery_coxm1():
    model = Network('Cox/M/1')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Coxian.fit_central(1, 0.99, 1.999))
    queue.setService(oclass, Exp(2))

    # Block 3: topology
    model.link(Network.serial_routing([source, queue, sink]))

    return model

if __name__ == '__main__':
    model = gallery_coxm1()
    print(f"Model: {model.getName()}")
