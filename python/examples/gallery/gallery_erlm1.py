#!/usr/bin/env python3
"""Gallery Example: Er/M/1 Queue"""

from line_solver import *

def gallery_erlm1():
    model = Network('Er/M/1')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Erlang.fit_mean_and_order(1, 5))
    queue.setService(oclass, Exp(2))

    # Block 3: topology
    model.link(Network.serial_routing([source, queue, sink]))

    return model

if __name__ == '__main__':
    model = gallery_erlm1()
    print(f"Model: {model.getName()}")
