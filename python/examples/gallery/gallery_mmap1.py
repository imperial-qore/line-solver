#!/usr/bin/env python3
"""Gallery Example: M/MAP/1 Queue"""

from line_solver import *

def gallery_mmap1(map=None):
    if map is None:
        map = MAP.rand()
        map = map.set_mean(0.5)

    model = Network('M/MAP/1')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp(1))
    queue.setService(oclass, map)

    # Block 3: topology
    model.link(Network.serial_routing([source, queue, sink]))

    return model

if __name__ == '__main__':
    model = gallery_mmap1()
    print(f"Model: {model.getName()}")
