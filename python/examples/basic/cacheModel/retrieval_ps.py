"""PS variant of retrieval_simple: a single processor-sharing retrieval station."""

# Ensure native line_solver is used (not python-wrapper)
import sys
import os
_native_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
if _native_path not in sys.path:
    sys.path.insert(0, _native_path)

import numpy as np
from line_solver import *


def retrieval_ps():
    access_prob = list(np.array([49, 49, 49, 49, 7, 1, 1]) / 205)  # per-item arrival probabilities (7 items)

    model = Network('PS Model')

    n = len(access_prob)               # number of items
    capacity = [6]                     # per-level cache capacity (total m = 6)

    source = Source(model, 'Source')
    cache_node = Cache(model, 'Cache', n, capacity, ReplacementStrategy.RR)
    queue = Queue(model, 'Queue', SchedStrategy.PS)
    sink = Sink(model, 'Sink')

    job_class = OpenClass(model, 'InitClass', 0)
    hit_class = OpenClass(model, 'HitClass', 0)
    miss_class = OpenClass(model, 'MissClass', 0)

    source.set_arrival(job_class, Exp(1))

    # Read class service at the retrieval queue = default per-item fetch service.
    queue.set_service(job_class, Exp(1.0))

    cache_node.set_read(job_class, DiscreteSampler(access_prob))
    cache_node.set_hit_class(job_class, hit_class)
    cache_node.set_miss_class(job_class, miss_class)

    # No service_rates and no routing_matrices: both inherited from the read class.
    cache_node.set_retrieval_system(job_class, miss_class, queue)

    P = model.init_routing_matrix()
    P.set(job_class, job_class, source, cache_node, 1.0)
    # Retrieval topology drawn once for the read class: cache <-> queue.
    P.set(job_class, job_class, cache_node, queue, 1.0)
    P.set(job_class, job_class, queue, cache_node, 1.0)
    P.set(hit_class, hit_class, cache_node, sink, 1.0)
    P.set(miss_class, miss_class, cache_node, sink, 1.0)
    model.link(P)

    return model, cache_node


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)
    model, cache_node = retrieval_ps()

    # Simulation
    print('SOLVER: SSA')
    print(SSA(model, samples=5000, method='serial', seed=1).get_avg_cache_table().to_string(index=False))
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
