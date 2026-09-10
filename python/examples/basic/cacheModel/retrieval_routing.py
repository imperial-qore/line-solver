"""Cache with probabilistic per-item routing through a retrieval system."""

import numpy as np
from line_solver import *


def retrieval_routing():
    access_prob = [0.6, 0.3, 0.1]          # per-item access probabilities (3 items)

    model = Network('Probabilistic Routing')

    n = len(access_prob)                   # number of items
    capacity = [2]                         # per-level cache capacity

    source = Source(model, 'Source')
    cache_node = Cache(model, 'Cache', n, capacity, ReplacementStrategy.FIFO)
    is_queue = Queue(model, 'IS Queue', SchedStrategy.INF)
    queue1 = Queue(model, 'Queue 1', SchedStrategy.FCFS)
    queue2 = Queue(model, 'Queue 2', SchedStrategy.FCFS)
    queues = [is_queue, queue1, queue2]
    sink = Sink(model, 'Sink')

    job_class = OpenClass(model, 'InitClass', 0)
    hit_class = OpenClass(model, 'HitClass', 0)
    miss_class = OpenClass(model, 'MissClass', 0)

    source.set_arrival(job_class, Exp(1))

    # Read class service per queue = default per-item fetch service (uniform over items).
    is_queue.set_service(job_class, Exp(2.0))
    queue1.set_service(job_class, Exp(3.0))
    queue2.set_service(job_class, Exp(3.0))

    cache_node.set_read(job_class, DiscreteSampler(access_prob))
    cache_node.set_hit_class(job_class, hit_class)
    cache_node.set_miss_class(job_class, miss_class)

    cache_node.set_retrieval_system(job_class, miss_class, queues)

    # Per-item routing over [IS(0), Queue1(1), Queue2(2), Cache(3)], applied via the
    # per-item routing methods. Index n_queues (last row/col, 0-based) is the cache;
    # row = from, col = to.
    routing_matrices = [
        np.array([[0.00, 0.50, 0.00, 0.50],   # from IS
                  [0.00, 0.00, 0.70, 0.30],    # from Queue1
                  [0.00, 0.00, 0.00, 1.00],    # from Queue2
                  [0.70, 0.30, 0.00, 0.00]]),  # from Cache
        np.array([[0.00, 0.30, 0.00, 0.70],
                  [0.00, 0.00, 0.50, 0.50],
                  [0.00, 0.00, 0.00, 1.00],
                  [0.20, 0.80, 0.00, 0.00]]),
        np.array([[0.00, 0.60, 0.00, 0.40],
                  [0.00, 0.00, 0.40, 0.60],
                  [0.00, 0.00, 0.00, 1.00],
                  [0.50, 0.50, 0.00, 0.00]]),
    ]

    n_queues = len(queues)
    for item in range(n):
        R = routing_matrices[item]
        for a in range(n_queues):
            cache_node.set_item_routing_prob(job_class, item, cache_node, queues[a], R[n_queues, a])
            cache_node.set_item_routing_prob(job_class, item, queues[a], cache_node, R[a, n_queues])
            for b in range(n_queues):
                cache_node.set_item_routing_prob(job_class, item, queues[a], queues[b], R[a, b])

    P = model.init_routing_matrix()
    P.set(job_class, job_class, source, cache_node, 1.0)
    P.set(hit_class, hit_class, cache_node, sink, 1.0)
    P.set(miss_class, miss_class, cache_node, sink, 1.0)
    model.link(P)

    return model, cache_node


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)
    model, cache_node = retrieval_routing()

    # Simulation
    print('SOLVER: SSA')
    print(SSA(model, samples=100000, method='serial', seed=1).get_avg_cache_table().to_string(index=False))
    print('SOLVER: LDES')
    print(LDES(model, samples=1000000, seed=1).get_avg_cache_table().to_string(index=False))

    # Analytical
    print('SOLVER: MVA')
    print(MVA(model).get_avg_cache_table().to_string(index=False))
    print('SOLVER: NC')
    print(NC(model).get_avg_cache_table().to_string(index=False))

    # Item-level cache occupancy
    print(MVA(model).get_avg_item_table())
    print(NC(model).get_avg_item_table())
