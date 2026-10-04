#!/usr/bin/env python3
"""Gallery Example: gallery_fj_open"""

from line_solver import *

def gallery_fj_open():
    """Open fork-join network (single class, two parallel tasks)."""
    model = Network('Fork-Join-Open')
    source = Source(model, 'Source')
    queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    fork = Fork(model, 'Fork')
    join = Join(model, 'Join', fork)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'class1')
    source.setArrival(oclass, Exp(0.05))
    queue1.setService(oclass, Exp(1.0))
    queue2.setService(oclass, Exp(2.0))
    P = model.init_routing_matrix()
    P.set(oclass, oclass, source, fork, 1.0)
    P.set(oclass, oclass, fork, queue1, 1.0)
    P.set(oclass, oclass, fork, queue2, 1.0)
    P.set(oclass, oclass, queue1, join, 1.0)
    P.set(oclass, oclass, queue2, join, 1.0)
    P.set(oclass, oclass, join, sink, 1.0)
    model.link(P)
    return model



if __name__ == '__main__':
    model = gallery_fj_open()
    print('Model built:', type(model).__name__)
