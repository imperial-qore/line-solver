"""
Cache with per-item storage costs (sizes) and per-list cost caps.

Each item i carries a storage cost sigma_i and list j may hold items of total
cost at most k_j. A promotion that would breach a cap serves the request
without changing the cache state. SolverNC evaluates the constrained
normalizing constant E(m,k) of Casale-Gast (IEEE/ACM ToN 29(2), 2021), Sec. IX;
SolverLDES simulates the same rule directly.
"""

from line_solver import *
import numpy as np


def cache_itemsize_costcap():

    model = Network('model')

    n = 6                          # number of items
    m = np.array([1, 1])           # two lists, one item each
    sizes = [1, 1, 1, 2, 2, 2]     # small and large items
    caps = [2, 1]                  # list 2 admits small items only

    source = Source(model, 'Source')
    cache_node = Cache(model, 'Cache', n, m, ReplacementStrategy.RR)
    sink = Sink(model, 'Sink')

    job_class = OpenClass(model, 'InitClass', 0)
    hit_class = OpenClass(model, 'HitClass', 0)
    miss_class = OpenClass(model, 'MissClass', 0)

    source.set_arrival(job_class, Exp(2))
    cache_node.set_read(job_class, DiscreteSampler((1.0 / n) * np.ones(n)))
    cache_node.set_item_sizes(sizes)
    cache_node.set_cost_caps(caps)

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
    model, cache_node = cache_itemsize_costcap()

    solver = NC(model, 'exact')
    solver.getAvgNodeTable()
    solver.getAvgCacheTable()
    solver.getAvgItemTable()

    print("Hit Ratio:", cache_node.get_hit_ratio())
    print("Mean per-list storage cost:", cache_node.get_list_cost())
