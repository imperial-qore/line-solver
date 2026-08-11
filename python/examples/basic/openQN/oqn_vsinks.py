"""
Open Network with Virtual Sinks (Multiple Sinks)

This example demonstrates:
- Multiple virtual sinks using Router nodes
- Different routing probabilities to different sinks
- Two classes with different exit probabilities
- Monitoring individual sink throughputs
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    source = Source(model, 'Source')
    queue = Queue(model, 'Queue1', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    vsink1 = Router(model, 'VSink1')
    vsink2 = Router(model, 'VSink2')

    ocl1 = OpenClass(model, 'Class1')
    ocl2 = OpenClass(model, 'Class2')

    source.set_arrival(ocl1, Exp(1.0))
    queue.set_service(ocl1, Exp(100.0))

    source.set_arrival(ocl2, Exp(1.0))
    queue.set_service(ocl2, Exp(100.0))

    P = model.init_routing_matrix()

    # Class 1 routing (60% to vsink1, 40% to vsink2)
    P.set(ocl1, ocl1, source, queue, 1.0)
    P.set(ocl1, ocl1, queue, vsink1, 0.6)
    P.set(ocl1, ocl1, queue, vsink2, 0.4)
    P.set(ocl1, ocl1, vsink1, sink, 1.0)
    P.set(ocl1, ocl1, vsink2, sink, 1.0)

    # Class 2 routing (10% to vsink1, 90% to vsink2)
    P.set(ocl2, ocl2, source, queue, 1.0)
    P.set(ocl2, ocl2, queue, vsink1, 0.1)
    P.set(ocl2, ocl2, queue, vsink2, 0.9)
    P.set(ocl2, ocl2, vsink1, sink, 1.0)
    P.set(ocl2, ocl2, vsink2, sink, 1.0)

    model.link(P)

    # Use getAvgNodeTable to see the throughputs of vsink1 and vsink2
    print('MVA Results:')
    avg_table = MVA(model).avg_table()
    print(avg_table)

    avg_node_table = MVA(model).avg_node_table()
    print('\nNode-level metrics (including virtual sinks):')
    print(avg_node_table)

    print('\nMAM Results:')
    avg_table_mam = MAM(model).avg_table()
    print(avg_table_mam)

    avg_node_table_mam = MAM(model).avg_node_table()
    print('\nMAM Node-level metrics:')
    print(avg_node_table_mam)

    print('\nNC Results:')
    avg_table_nc = NC(model).avg_table()
    print(avg_table_nc)

    avg_node_table_nc = NC(model).avg_node_table()
    print('\nNC Node-level metrics:')
    print(avg_node_table_nc)
