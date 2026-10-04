"""
Cache with h-LRU / LRU(m) Replacement Policy

h LRU lists of capacities m[0..h-1]; a miss inserts the item at the head of
list 1, a hit in list l exchanges the item with the tail of list l+1. Exact
in CTMC, simulated in SSA/LDES; MVA uses the characteristic-time (TTL)
approximation of Gast and Van Houdt (SIGMETRICS 2015), which reduces to the
Che approximation for h=1.
"""

from line_solver import *
import numpy as np


def cache_replc_hlru():
    model = Network('model')

    n = 6         # number of items
    m = [2, 1]    # list capacities: list 1 holds 2 items, list 2 holds 1

    source = Source(model, 'Source')
    cache_node = Cache(model, 'Cache', n, m, ReplacementStrategy.HLRU)
    sink = Sink(model, 'Sink')

    job_class = OpenClass(model, 'InitClass', 0)
    hit_class = OpenClass(model, 'HitClass', 0)
    miss_class = OpenClass(model, 'MissClass', 0)

    source.set_arrival(job_class, Exp(1))

    p_access = Zipf(1.2, n)  # Zipf-like item references
    cache_node.set_read(job_class, p_access)

    cache_node.set_hit_class(job_class, hit_class)
    cache_node.set_miss_class(job_class, miss_class)

    P = model.init_routing_matrix()
    P.set(job_class, job_class, source, cache_node, 1.0)
    P.set(hit_class, hit_class, cache_node, sink, 1.0)
    P.set(miss_class, miss_class, cache_node, sink, 1.0)

    model.link(P)
    return model, cache_node


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)
    model, cache_node = cache_replc_hlru()

    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(model, keep=False, cutoff=1))  # exact

    model.reset()
    solver = np.append(solver, MVA(model))  # TTL approximation

    model.reset()
    solver = np.append(solver, SSA(model, samples=10000, seed=23000))

    avg_node_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_node_table[s] = solver[s].avg_node_table()
        print(avg_node_table[s])
