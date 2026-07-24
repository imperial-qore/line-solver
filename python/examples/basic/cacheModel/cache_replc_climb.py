"""
Cache with CLIMB (transposition) Replacement Policy

On a hit an item moves up one position; on a miss it enters at the tail.
Exact in CTMC, simulated in SSA/LDES. Not product-form, so MVA/NC/FLD reject
it (use CTMC/SSA/LDES).
"""

import sys
import os
_native_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
if _native_path not in sys.path:
    sys.path.insert(0, _native_path)

from line_solver import *
import numpy as np


def cache_replc_climb():
    model = Network('model')

    n = 5  # number of items
    m = 2  # cache capacity

    delay = Delay(model, 'Delay')
    cache_node = Cache(model, 'Cache', n, m, ReplacementStrategy.CLIMB)

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
    model, cache_node = cache_replc_climb()
    print(SolverCTMC(model, keep=False).get_avg_node_table())
