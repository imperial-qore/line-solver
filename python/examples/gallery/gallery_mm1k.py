#!/usr/bin/env python3
"""Gallery Example: gallery_mm1k"""

from line_solver import *

def gallery_mm1k(K=3):
    """M/M/1/K queue with finite capacity K (blocking / loss)."""
    model = Network('M/M/1/K')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    queue.setNumberOfServers(1)
    queue.setCapacity(K)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1', 0)
    source.setArrival(oclass, Exp(0.8))
    queue.setService(oclass, Exp(1.0))
    P = model.init_routing_matrix()
    P.set(oclass, oclass, source, queue, 1.0)
    P.set(oclass, oclass, queue, sink, 1.0)
    model.link(P)
    return model



if __name__ == '__main__':
    model = gallery_mm1k()
    print('Model built:', type(model).__name__)
