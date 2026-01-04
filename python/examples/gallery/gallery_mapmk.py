#!/usr/bin/env python3
"""Gallery Example: MAP/M/k Queue"""

from line_solver import *

def gallery_mapmk(map=None, k=2):
    if map is None:
        map = MAP.rand()

    model = Network('MAP/M/k')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, map)
    queue.setService(oclass, Exp(2))
    queue.setNumberOfServers(k)

    # Block 3: topology
    model.link(Network.serial_routing([source, queue, sink]))

    return model

if __name__ == '__main__':
    model = gallery_mapmk()
    print(f"Model: {model.getName()}")
