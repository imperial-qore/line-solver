import numpy as np
from line_solver import Network, Source, Queue, Sink, Exp
from line_solver import ClosedClass, OpenClass, SchedStrategy


def infer_quick_model_rnn(is_open, stations, classes, servers=None, jobs=None, routing=None):
    """Generate simple queueing network for RNN-based estimation.

    Similar to infer_quick_model but initializes queue states randomly.

    Args:
        is_open: True for open network, False for closed
        stations: list of SchedStrategy values
        classes: 2-D array of service demands (num_classes x num_stations)
        servers: list of server counts per station (default: all 1)
        jobs: list of job counts per class (default: all 1)
        routing: routing matrix (optional)

    Returns:
        model: LINE Network model

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    classes = np.atleast_2d(np.asarray(classes, dtype=float))
    num_stations = len(stations)
    num_classes = classes.shape[0]

    if servers is None:
        servers = [1] * num_stations
    if jobs is None:
        jobs = [1] * num_classes

    model = Network('quickModel')

    if is_open:
        source = Source(model, 'mySource')
        sink = Sink(model, 'mySink')

    queue_nodes = []
    for i in range(num_stations):
        q = Queue(model, f'QueueStation{i + 1}', stations[i])
        q.setNumberOfServers(servers[i])
        queue_nodes.append(q)

    jobclasses = []
    for c in range(num_classes):
        if not is_open:
            jc = ClosedClass(model, f'Class{c + 1}', int(jobs[c]), queue_nodes[0])
        else:
            jc = OpenClass(model, f'Class{c + 1}')
        jobclasses.append(jc)

        for i in range(num_stations):
            queue_nodes[i].setService(jc, Exp.fitMean(classes[c, i]))

    if is_open:
        all_nodes = [source] + queue_nodes + [sink]
        P = model.initRoutingMatrix()
        for c in range(num_classes):
            for idx in range(len(all_nodes) - 1):
                P.set(jobclasses[c], jobclasses[c], all_nodes[idx], all_nodes[idx + 1], 1.0)
        model.link(P)
    else:
        P = model.initRoutingMatrix()
        for c in range(num_classes):
            if routing is not None:
                for i in range(num_stations):
                    for j in range(num_stations):
                        if routing[c][i][j] > 0:
                            P.set(jobclasses[c], jobclasses[c], queue_nodes[i], queue_nodes[j], routing[c][i][j])
            else:
                for i in range(num_stations):
                    next_i = (i + 1) % num_stations
                    P.set(jobclasses[c], jobclasses[c], queue_nodes[i], queue_nodes[next_i], 1.0)
        model.link(P)

    return model
