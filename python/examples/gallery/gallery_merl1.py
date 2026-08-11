#!/usr/bin/env python3
"""Gallery Example: M/E/1 Queue"""

from line_solver import *

def gallery_merl1():
    model = Network('M/E/1')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, Erlang.fit_mean_and_order(0.5, 2))

    # Block 3: topology
    model.link(Network.serial_routing([source, queue, sink]))

    return model

if __name__ == '__main__':
    model = gallery_merl1()
    print(f"Model: {model.getName()}")
