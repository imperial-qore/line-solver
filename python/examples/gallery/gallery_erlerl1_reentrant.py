#!/usr/bin/env python3
"""Gallery Example: Erl/Erl/1-Reentrant Queue"""

from line_solver import *

def gallery_erlerl1_reentrant():
    model = Network('Erl/Erl/1-Reentrant')
    n = 5

    # Block 1: nodes
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    # Block 2: classes
    oclass1 = OpenClass(model, 'Class1')
    oclass2 = OpenClass(model, 'Class2')
    source.setArrival(oclass1, Erlang.fit_mean_and_order(1, n))
    source.setArrival(oclass2, Disabled())
    queue.setService(oclass1, Erlang.fit_mean_and_order(0.5, n))
    queue.setService(oclass2, Exp(3))

    # Block 3: topology
    P = model.init_routing_matrix()
    P.set(oclass1, oclass1, source, queue, 1)
    P.set(oclass1, oclass2, queue, queue, 0.5)
    P.set(oclass2, oclass2, queue, sink, 0.5)
    model.link(P)

    return model

if __name__ == '__main__':
    model = gallery_erlerl1_reentrant()
    print(f"Model: {model.getName()}")
