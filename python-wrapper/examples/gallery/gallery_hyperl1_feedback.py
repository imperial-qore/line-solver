#!/usr/bin/env python3
"""Gallery Example: Hyper/Erl/1-Feedback Queue"""

from line_solver import *

def gallery_hyperl1_feedback():
    model = Network('Hyper/Erl/1-Feedback')

    # Block 1: nodes
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    # Block 2: classes
    oclass1 = OpenClass(model, 'Class1')
    source.setArrival(oclass1, HyperExp.fit_mean_and_scv(1, 64))
    queue.setService(oclass1, Erlang.fit_mean_and_order(0.5, 5))

    # Block 3: topology
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass1, queue, queue, 0.9)
    P.set(oclass1, oclass1, queue, sink, 0.1)
    model.link(P)

    return model

if __name__ == '__main__':
    model = gallery_hyperl1_feedback()
    print(f"Model: {model.getName()}")
