#!/usr/bin/env python3
"""Gallery Example: D/M/1 Queue"""

from line_solver import *

def gallery_dm1():
    model = Network('D/M/1')

    # Block 1: nodes
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    # Block 2: classes
    oclass1 = OpenClass(model, 'Class1')
    source.setArrival(oclass1, Det(1))
    queue.setService(oclass1, Exp(2))

    # Block 3: topology
    model.link(Network.serial_routing([source, queue, sink]))

    return model

if __name__ == '__main__':
    model = gallery_dm1()
    print(f"Model: {model.getName()}")
