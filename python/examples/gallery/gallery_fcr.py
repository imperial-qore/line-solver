#!/usr/bin/env python3
"""Gallery Example: gallery_fcr"""

from line_solver import *

def gallery_fcr(K=3):
    """Finite capacity region with dropping around a single queue (like M/M/1/K)."""
    model = Network('FCR-Dropping')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1', 0)
    source.setArrival(oclass, Exp(0.8))
    queue.setService(oclass, Exp(1.0))
    P = model.init_routing_matrix()
    P.set(oclass, oclass, source, queue, 1.0)
    P.set(oclass, oclass, queue, sink, 1.0)
    model.link(P)
    fcr = model.add_region([queue])
    fcr.setGlobalMaxJobs(K)
    fcr.setDropRule(oclass, True)
    return model



if __name__ == '__main__':
    model = gallery_fcr()
    print('Model built:', type(model).__name__)
