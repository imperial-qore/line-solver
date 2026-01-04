#!/usr/bin/env python3
"""Gallery Example: M/MAP/1 Multiclass Queue"""

from line_solver import *

def gallery_mmap1_multiclass(map1=None, map2=None):
    if map1 is None:
        map1 = MAP.rand(2)
        map1 = map1.set_mean(0.5)
    if map2 is None:
        map2 = MAP.rand(3)
        map2 = map2.set_mean(0.5)

    model = Network('M/MAP/1')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass1 = OpenClass(model, 'myClass1')
    source.setArrival(oclass1, Exp(0.35 / map1.get_mean()))
    queue.setService(oclass1, map1)

    oclass2 = OpenClass(model, 'myClass2')
    source.setArrival(oclass2, Exp(0.15 / map2.get_mean()))
    queue.setService(oclass2, map2)

    # Block 3: topology
    P = model.init_routing_matrix()
    P[1] = Network.serial_routing([source, queue, sink])
    P[2] = Network.serial_routing([source, queue, sink])
    model.link(P)

    return model

if __name__ == '__main__':
    model = gallery_mmap1_multiclass()
    print(f"Model: {model.getName()}")
