"""
Trace-Driven Open Network

This example demonstrates:
- Trace-driven arrivals using Replayer
- Service times from trace file
- Simple open network: Source → Queue → Sink
"""

from line_solver import *
import os

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    jobclass = OpenClass(model, 'OpenClass', 0)

    source.set_arrival(jobclass, Exp(1))

    # Get path to trace file in the same directory as the MATLAB example
    trace_file = '/home/gcasale/Dropbox/code/line-dev.git/matlab/examples/basic/openQN/example_trace.txt'
    queue.set_service(jobclass, Replayer(trace_file))

    model.link(Network.serial_routing([source, queue, sink]))

    # Run solvers
    avg_table_1 = JMT(model, seed=23000).avg_table()
    print('JMT Result:')
    print(avg_table_1)

    avg_table_2 = DES(model, seed=23000).avg_table()
    print('\nDES Result:')
    print(avg_table_2)
