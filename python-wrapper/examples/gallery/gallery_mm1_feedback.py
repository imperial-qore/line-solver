#!/usr/bin/env python3
"""Gallery Example: M/M/1-Feedback Queue"""

from line_solver import *

def gallery_mm1_feedback(p=1/3):
    model = Network('M/M/1-Feedback')

    # Block 1: nodes
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    # Block 2: classes
    oclass1 = OpenClass(model, 'Class1')
    source.setArrival(oclass1, Exp.fit_mean(1))
    queue.setService(oclass1, Exp.fit_mean(0.5))

    # Block 3: topology
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass1, queue, queue, p)
    P.set(oclass1, oclass1, queue, sink, 1-p)
    model.link(P)

    return model

if __name__ == '__main__':
    model = gallery_mm1_feedback()
    print(f"Model: {model.getName()}")
