#!/usr/bin/env python3
"""Gallery Example: MAP/M/1 Queue"""

from line_solver import *

def gallery_mapm1(map=None, seed=23000):
    if map is None:
        import numpy as np
        D0 = np.array([[-0.6984901916396979, 0.45234650636128054],
                        [0.34690024319398277, -0.8194057961021199]])
        D1 = np.array([[0.2125067546435463, 0.033636930634871165],
                        [0.4441520099524867, 0.028353542955650513]])
        map = MAP(D0, D1)

    model = Network('MAP/M/1')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    source.setArrival(oclass, map)
    queue.setService(oclass, Exp(2))

    # Block 3: topology
    model.link(Network.serial_routing([source, queue, sink]))

    return model

if __name__ == '__main__':
    model = gallery_mapm1()
    print(f"Model: {model.getName()}")
