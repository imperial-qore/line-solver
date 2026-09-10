"""
Cache with q-LRU Replacement Policy

On a miss the item is admitted (LRU head insert) with probability q, otherwise
it passes through uncached. Admission filtering can raise the hit ratio over
plain LRU under skewed popularity. Exact in CTMC, simulated in SSA/LDES. Not
product-form, so MVA/NC/FLD reject it (use CTMC/SSA/LDES).
"""

from line_solver import *
import numpy as np


def cache_replc_qlru():
    model = Network('model')

    n = 5  # number of items
    m = 2  # cache capacity

    delay = Delay(model, 'Delay')
    cache_node = Cache(model, 'Cache', n, m, ReplacementStrategy.QLRU)
    cache_node.set_admission_prob(0.5)  # admit a missed item with probability q=0.5

    job_class = ClosedClass(model, 'JobClass', 1, delay, 0)
    hit_class = ClosedClass(model, 'HitClass', 0, delay, 0)
    miss_class = ClosedClass(model, 'MissClass', 0, delay, 0)

    delay.set_service(job_class, Exp(1))
    cache_node.set_read(job_class, Zipf(1.2, n))
    cache_node.set_hit_class(job_class, hit_class)
    cache_node.set_miss_class(job_class, miss_class)

    P = model.init_routing_matrix()
    P.set(job_class, job_class, delay, cache_node, 1.0)
    P.set(hit_class, job_class, cache_node, delay, 1.0)
    P.set(miss_class, job_class, cache_node, delay, 1.0)
    model.link(P)
    return model, cache_node


if __name__ == "__main__":
    model, cache_node = cache_replc_qlru()
    print(SolverCTMC(model, 'exact', keep=False).get_avg_node_table())
