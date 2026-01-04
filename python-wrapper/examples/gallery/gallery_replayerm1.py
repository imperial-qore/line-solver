#!/usr/bin/env python3
"""Gallery Example: Trace/M/1 Queue"""

import os
from line_solver import *

def gallery_replayerm1(filename=None):
    if filename is None:
        # Find example_trace.txt in the package
        import line_solver
        pkg_dir = os.path.dirname(line_solver.__file__)
        filename = os.path.join(pkg_dir, 'example_trace.txt')
        # If not found, try the matlab directory
        if not os.path.exists(filename):
            filename = '/home/gcasale/Dropbox/code/line-dev.git/matlab/example_trace.txt'

    model = Network('Trace/M/1')

    # Block 1: nodes
    source = Source(model, 'mySource')
    queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'mySink')

    # Block 2: classes
    oclass = OpenClass(model, 'myClass')
    replayer = Replayer(filename)
    source.setArrival(oclass, replayer)
    queue.setService(oclass, Exp(3 / replayer.get_mean()))

    # Block 3: topology
    model.link(Network.serial_routing([source, queue, sink]))

    return model

if __name__ == '__main__':
    model = gallery_replayerm1()
    print(f"Model: {model.getName()}")
