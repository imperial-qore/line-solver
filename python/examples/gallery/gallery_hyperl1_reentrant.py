#!/usr/bin/env python3
"""Gallery Example: Hyper/Erl/1-Reentrant Queue"""

from line_solver import *

def gallery_hyperl1_reentrant():
    model = Network('Hyper/Erl/1-Reentrant')

    # Block 1: nodes
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    # Block 2: classes
    oclass1 = OpenClass(model, 'Class1')
    oclass2 = OpenClass(model, 'Class2')
    source.setArrival(oclass1, HyperExp.fit_mean_and_scv(1, 64))
    source.setArrival(oclass2, Disabled())
    queue.setService(oclass1, Erlang.fit_mean_and_order(0.5, 5))
    queue.setService(oclass2, Exp(3))

    # Block 3: topology
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass2, queue, queue, 1)
    P.set(oclass2, oclass2, queue, sink, 1)
    model.link(P)

    return model

if __name__ == '__main__':
    model = gallery_hyperl1_reentrant()
    print(f"Model: {model.getName()}")
