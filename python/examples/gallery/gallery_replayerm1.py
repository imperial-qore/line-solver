#!/usr/bin/env python3
"""Gallery Example: Trace/M/1 Queue"""

import os
from line_solver import *

def gallery_replayerm1(filename=None):
    if filename is None:
        # Try to find example_trace.txt in common locations
        script_dir = os.path.dirname(os.path.abspath(__file__))
        possible_paths = [
            os.path.join(script_dir, '..', 'gettingstarted', 'example_trace.txt'),
            os.path.join(script_dir, '..', 'basic', 'openQN', 'example_trace.txt'),
        ]
        for path in possible_paths:
            if os.path.exists(path):
                filename = path
                break
        if filename is None:
            raise FileNotFoundError("example_trace.txt not found in expected locations")

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
