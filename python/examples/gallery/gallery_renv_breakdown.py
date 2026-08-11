#!/usr/bin/env python3
"""Gallery Example: gallery_renv_breakdown"""

from line_solver import *

def gallery_renv_breakdown():
    """Random environment: single server with breakdown/repair (UP/DOWN stages).

    Returns an Environment whose base model is an M/M/1 queue alternating
    between an UP stage (fast service) and a DOWN stage (degraded service).
    """
    model = Network('ServerWithFailures')
    source = Source(model, 'Arrivals')
    queue = Queue(model, 'Server', SchedStrategy.FCFS)
    sink = Sink(model, 'Departures')
    jobclass = OpenClass(model, 'Jobs')
    source.setArrival(jobclass, Exp(0.8))
    queue.setService(jobclass, Exp(2.0))
    queue.setNumberOfServers(1)
    P = model.init_routing_matrix()
    P.set(jobclass, jobclass, source, queue, 1.0)
    P.set(jobclass, jobclass, queue, sink, 1.0)
    model.link(P)
    env = Environment('ServerEnv')
    env.add_node_failure_repair(model, queue, Exp(0.1), Exp(1.0), Exp(0.5))
    env.init()
    return env

if __name__ == '__main__':
    model = gallery_renv_breakdown()
    print('Model built:', type(model).__name__)
