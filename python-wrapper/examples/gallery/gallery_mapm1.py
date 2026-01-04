#!/usr/bin/env python3
"""Gallery Example: MAP/M/1 Queue"""

from line_solver import *

def gallery_mapm1(map=None):
    if map is None:
        map = MAP.rand()

    model = Network('MAP/M/1')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, map)
    queue.setService(oclass, Exp(2))

    # Block 3: topology
    model.link(Network.serial_routing([source, queue, sink]))

    return model

if __name__ == '__main__':
    model = gallery_mapm1()
    print(f"Model: {model.getName()}")
