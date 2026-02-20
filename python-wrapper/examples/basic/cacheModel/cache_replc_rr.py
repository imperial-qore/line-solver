"""
Cache with Random Replacement (RR) Policy

This example demonstrates:
- Open network with cache node
- RR (Random Replacement) strategy
- 5 items, cache capacity of 2
- Zipf access pattern (skewed distribution)
- Hit ratio analysis
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

    source = Source(model, 'Source')
    cache_node = Cache(model, 'Cache', n, m, ReplacementStrategy.RR)
    sink = Sink(model, 'Sink')

    job_class = OpenClass(model, 'InitClass', 0)
    hit_class = OpenClass(model, 'HitClass', 0)
    miss_class = OpenClass(model, 'MissClass', 0)

    source.set_arrival(job_class, Exp(2))

    # Zipf-like item references
    p_access = Zipf(1.4, n)
    cache_node.set_read(job_class, p_access)

    cache_node.set_hit_class(job_class, hit_class)
    cache_node.set_miss_class(job_class, miss_class)

    P = model.init_routing_matrix()
    P.set(job_class, job_class, source, cache_node, 1.0)
    P.set(hit_class, hit_class, cache_node, sink, 1.0)
    P.set(miss_class, miss_class, cache_node, sink, 1.0)

    model.link(P)

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
        print(f'\nSOLVER: {solver[s].get_name().replace("Solver", "")}')
        avg_node_table[s] = solver[s].avg_node_table()
        print(avg_node_table[s])

    print(f'\nHit Ratio: {cache_node.get_hit_ratio()}')
    print(f'Miss Ratio: {cache_node.get_miss_ratio()}')
