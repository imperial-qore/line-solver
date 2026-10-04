#!/usr/bin/env python3
"""Gallery Example: gallery_fj_closed"""

from line_solver import *

def gallery_fj_closed():
    """Closed fork-join network (single class, two parallel tasks)."""
    model = Network('Fork-Join-Closed')
    delay = Delay(model, 'Delay')
    queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
    fork = Fork(model, 'Fork')
    join = Join(model, 'Join', fork)
    oclass = ClosedClass(model, 'class1', 5, delay)
    delay.setService(oclass, Exp(1.0))
    queue1.setService(oclass, Exp(1.0))
    queue2.setService(oclass, Exp(1.0))
    P = model.init_routing_matrix()
    P.set(oclass, oclass, delay, fork, 1.0)
    P.set(oclass, oclass, fork, queue1, 1.0)
    P.set(oclass, oclass, fork, queue2, 1.0)
    P.set(oclass, oclass, queue1, join, 1.0)
    P.set(oclass, oclass, queue2, join, 1.0)
    P.set(oclass, oclass, join, delay, 1.0)
    model.link(P)
    return model



if __name__ == '__main__':
    model = gallery_fj_closed()
    print('Model built:', type(model).__name__)
