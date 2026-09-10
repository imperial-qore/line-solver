#!/usr/bin/env python3
"""Gallery Example: Er/M/1-PS Queue"""

from line_solver import *

def gallery_erlm1ps():
    model = Network('Er/M/1-PS')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.PS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Erlang.fit_mean_and_order(1, 5))
    queue.setService(oclass, Exp(2))

    # Block 3: topology
    model.link(Network.serial_routing([source, queue, sink]))

    return model

if __name__ == '__main__':
    model = gallery_erlm1ps()
    print(f"Model: {model.getName()}")
