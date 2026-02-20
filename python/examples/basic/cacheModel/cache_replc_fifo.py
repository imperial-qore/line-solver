"""
Cache with FIFO Replacement Policy

This example demonstrates:
- Closed network with cache node
- FIFO (First In First Out) replacement strategy
- 5 items, cache capacity of 2
- Uniform access pattern
"""

# Ensure native line_solver is used (not python-wrapper)
import sys
import os
_native_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
if _native_path not in sys.path:
    sys.path.insert(0, _native_path)

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    n = 5  # number of items
    m = 2  # cache capacity

    delay = Delay(model, 'Delay')
    cache_node = Cache(model, 'Cache', n, m, ReplacementStrategy.FIFO)

    job_class = ClosedClass(model, 'JobClass', 1, delay, 0)
    hit_class = ClosedClass(model, 'HitClass', 0, delay, 0)
    miss_class = ClosedClass(model, 'MissClass', 0, delay, 0)

    delay.set_service(job_class, Exp(1))

    # Uniform item references
    p_access = [1.0 / n] * n
    cache_node.set_read(job_class, DiscreteSampler(p_access))

    cache_node.set_hit_class(job_class, hit_class)
    cache_node.set_miss_class(job_class, miss_class)

    P = model.init_routing_matrix()
    P.set(job_class, job_class, delay, cache_node, 1.0)
    P.set(hit_class, job_class, cache_node, delay, 1.0)
    P.set(miss_class, job_class, cache_node, delay, 1.0)

    model.link(P)

    # Run multiple solvers
    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(model, keep=False))

    model.reset()
    solver = np.append(solver, SSA(model, samples=10000, verbose=True, seed=23000))

    solver = np.append(solver, MVA(model))

    avg_node_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_node_table[s] = solver[s].avg_node_table()
        print(avg_node_table[s])

    print(f'\nHit Ratio: {cache_node.get_hit_ratio()}')
    print(f'Miss Ratio: {cache_node.get_miss_ratio()}')
