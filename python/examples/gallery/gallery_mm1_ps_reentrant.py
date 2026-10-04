#!/usr/bin/env python3
"""Gallery Example: M/M/1-PS-Reentrant Queue"""

from line_solver import *

def gallery_mm1_ps_reentrant():
    model = Network('M/M/1')

    # Block 1: nodes
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.PS)
    sink = Sink(model, 'Sink')

    # Block 2: classes
    oclass1 = OpenClass(model, 'Class1')
    oclass2 = OpenClass(model, 'Class2')
    source.setArrival(oclass1, Exp(1))
    source.setArrival(oclass2, Disabled())
    queue.setService(oclass1, Exp(2))
    queue.setService(oclass2, Exp(3))

    # Block 3: topology
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass2, queue, queue, 1)
    P.set(oclass2, oclass2, queue, sink, 1)
    model.link(P)

    return model

if __name__ == '__main__':
    model = gallery_mm1_ps_reentrant()
    print(f"Model: {model.getName()}")
