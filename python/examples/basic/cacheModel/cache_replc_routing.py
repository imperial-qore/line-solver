"""
Cache with State-Dependent Routing

This example demonstrates:
- Cache with routing based on hit/miss classes
- Different delays for hits and misses
- Router node for probabilistic routing
- Open network with FIFO replacement
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    n = 5  # number of items
    m = 2  # cache capacity

    source = Source(model, 'Source')
    cache_node = Cache(model, 'Cache', n, m, ReplacementStrategy.FIFO)
    router_node = Router(model, 'Router')
    delay1 = Delay(model, 'Delay1')
    delay2 = Delay(model, 'Delay2')
    sink = Sink(model, 'Sink')

    job_class = OpenClass(model, 'InitClass', 0)
    hit_class = OpenClass(model, 'HitClass', 0)
    miss_class = OpenClass(model, 'MissClass', 0)

    source.set_arrival(job_class, Exp(2))

    # Different service times for hits and misses
    delay1.set_service(hit_class, Exp(10))
    delay1.set_service(miss_class, Exp(1))

    delay2.set_service(hit_class, Exp(20))
    delay2.set_service(miss_class, Exp(2))

    # Uniform item references
    p_access = [1.0 / n] * n
    cache_node.set_read(job_class, DiscreteSampler(p_access))

    cache_node.set_hit_class(job_class, hit_class)
    cache_node.set_miss_class(job_class, miss_class)

    # Set up routing
    model.add_link(source, cache_node)
    model.add_link(cache_node, router_node)
    model.add_link(router_node, delay1)
    model.add_link(router_node, delay2)
    model.add_link(delay1, sink)
    model.add_link(delay2, sink)

    source.set_prob_routing(job_class, cache_node, 1.0)

    cache_node.set_prob_routing(hit_class, router_node, 1.0)
    cache_node.set_prob_routing(miss_class, router_node, 1.0)

    router_node.set_routing(hit_class, RoutingStrategy.RAND)
    router_node.set_routing(miss_class, RoutingStrategy.RAND)

    delay1.set_prob_routing(hit_class, sink, 1.0)
    delay1.set_prob_routing(miss_class, sink, 1.0)

    delay2.set_prob_routing(hit_class, sink, 1.0)
    delay2.set_prob_routing(miss_class, sink, 1.0)

    # Run multiple solvers
    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(model, keep=False, cutoff=1))

    model.reset()
    solver = np.append(solver, SSA(model, samples=10000, verbose=True, seed=23000))

    model.reset()
    solver = np.append(solver, MVA(model))

    model.reset()
    solver = np.append(solver, NC(model))

    avg_node_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_node_table[s] = solver[s].avg_node_table()

    print(f'\nHit Ratio: {cache_node.get_hit_ratio()}')
    print(f'Miss Ratio: {cache_node.get_miss_ratio()}')
