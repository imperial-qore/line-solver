#!/usr/bin/env python3
"""Gallery Example: M/D/k Queue"""

from line_solver import *

def gallery_mdk(k=2):
    model = Network('M/D/k')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, Exp.fit_mean(1))
    queue.setService(oclass, Det(2/k))
    queue.setNumberOfServers(k)

    # Block 3: topology
    model.link(Network.serial_routing([source, queue, sink]))

    return model

if __name__ == '__main__':
    model = gallery_mdk()
    print(f"Model: {model.getName()}")
