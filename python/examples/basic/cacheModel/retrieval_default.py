"""Retrieval system whose per-item miss routing AND service are taken, by default,
from the read class: routing from the read class's edges among the retrieval queues
in the top-level routing matrix P, and service from the read class's service
distribution at each queue. Per-item overrides (set_item_routing_prob (with the cache as source/dest) / Queue.set_item_service_rate) then reconfigure one
item at finer granularity."""

from line_solver import *


def retrieval_default():
    access_prob = [0.6, 0.3, 0.1]          # per-item access probabilities (3 items)

    model = Network('DelayedHits')

    n = len(access_prob)                   # number of items
    capacity = [1]                         # per-level cache capacity

    source = Source(model, 'Source')
    cache_node = Cache(model, 'Cache', n, capacity, ReplacementStrategy.FIFO)
    # PS retrieval stations: per-item (class-dependent) service rates are admissible
    # in the analytical retrieval algorithm (FCFS/SIRO would require identical rates).
    queue1 = Queue(model, 'Queue_1', SchedStrategy.PS)
    queue2 = Queue(model, 'Queue_2', SchedStrategy.PS)
    sink = Sink(model, 'Sink')

    job_class = OpenClass(model, 'InitClass', 0)
    hit_class = OpenClass(model, 'HitClass', 0)
    miss_class = OpenClass(model, 'MissClass', 0)

    source.set_arrival(job_class, Exp(1))

    # Read class service at each retrieval queue = default per-item fetch service.
    queue1.set_service(job_class, Exp(2.0))
    queue2.set_service(job_class, Exp(3.0))

    cache_node.set_read(job_class, DiscreteSampler(access_prob))
    cache_node.set_hit_class(job_class, hit_class)
    cache_node.set_miss_class(job_class, miss_class)

    # No service_rates and no routing_matrices: both inherited from the read class.
    cache_node.set_retrieval_system(job_class, miss_class, [queue1, queue2])

    # Item-level overrides for item 0: skip Queue_2 and fetch faster at Queue_1.
    cache_node.set_item_routing_prob(job_class, 0, queue1, queue2, 0.0)  # delete default edge
    cache_node.set_item_routing_prob(job_class, 0, queue1, cache_node, 1.0)  # exit after Queue_1
    queue1.set_item_service_rate(cache_node, job_class, 0, 5.0)                 # faster item-0 fetch

    P = model.init_routing_matrix()
    P.set(job_class, job_class, source, cache_node, 1.0)
    # Default retrieval topology, drawn once for the read class: cache -> Q1 -> Q2 -> cache
    P.set(job_class, job_class, cache_node, queue1, 1.0)   # entry into retrieval
    P.set(job_class, job_class, queue1, queue2, 1.0)       # Queue_1 -> Queue_2
    P.set(job_class, job_class, queue2, cache_node, 1.0)   # exit back to cache
    P.set(hit_class, hit_class, cache_node, sink, 1.0)
    P.set(miss_class, miss_class, cache_node, sink, 1.0)
    model.link(P)

    return model, cache_node


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)
    model, cache_node = retrieval_default()

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
