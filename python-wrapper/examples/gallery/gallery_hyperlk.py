#!/usr/bin/env python3
"""Gallery Example: Hyper/Erl/k Queue"""

from line_solver import *

def gallery_hyperlk(k=2):
    model = Network('Hyper/Erl/k')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, HyperExp.fit_mean_and_scv_balanced(1/1.8, 4))
    queue.setService(oclass, Erlang.fit_mean_and_scv(1, 0.25))
    queue.setNumberOfServers(k)

    # Block 3: topology
    model.link(Network.serial_routing([source, queue, sink]))

    return model

if __name__ == '__main__':
    model = gallery_hyperlk()
    print(f"Model: {model.getName()}")
