#!/usr/bin/env python3
"""Gallery Example: M[2]/M[2]/1-PS Queue"""

from line_solver import *

def gallery_mm1_ps_multiclass():
    model = Network('M[2]/M[2]/1-PS')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.PS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass1 = OpenClass(model, 'myClass1')
    source.setArrival(oclass1, Exp(1))
    queue.setService(oclass1, Exp(4))

    oclass2 = OpenClass(model, 'myClass2')
    source.setArrival(oclass2, Exp(0.5))
    queue.setService(oclass2, Exp(4))

    # Block 3: topology
    P = model.init_routing_matrix()
    P[1] = Network.serial_routing([source, queue, sink])
    P[2] = Network.serial_routing([source, queue, sink])
    model.link(P)

    return model

if __name__ == '__main__':
    model = gallery_mm1_ps_multiclass()
    print(f"Model: {model.getName()}")
