"""
Compare Cache Replacement Strategies

This example demonstrates:
- Comparison of RR, FIFO, and LRU replacement strategies
- Open network with Zipf access pattern
- Hit ratio comparison across different strategies
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    n = 5
    m = 2
    alpha = 1.0

    repl_strat = [ReplacementStrategy.RR, ReplacementStrategy.FIFO, ReplacementStrategy.LRU]

    ctmc_hit_ratio = []
    mva_hit_ratio = []
    nc_hit_ratio = []

    for s in range(len(repl_strat)):
        model = Network('model')
        source = Source(model, 'Source')
        cache_node = Cache(model, 'Cache', n, m, repl_strat[s])
        sink = Sink(model, 'Sink')

        job_class = OpenClass(model, 'InitClass', 0)
        hit_class = OpenClass(model, 'HitClass', 0)
        miss_class = OpenClass(model, 'MissClass', 0)

        source.set_arrival(job_class, Exp(2))

        p_access = Zipf(alpha, n)
        cache_node.set_read(job_class, p_access)

        cache_node.set_hit_class(job_class, hit_class)
        cache_node.set_miss_class(job_class, miss_class)

        P = model.init_routing_matrix()
        P.set(job_class, job_class, source, cache_node, 1.0)
        P.set(hit_class, hit_class, cache_node, sink, 1.0)
        P.set(miss_class, miss_class, cache_node, sink, 1.0)

        model.link(P)

        # CTMC solver
        solver_ctmc = CTMC(model, keep=False, cutoff=1, verbose=False)
        solver_ctmc.avg_node_table()
        hr = cache_node.get_hit_ratio()
        ctmc_hit_ratio.append(hr if hr is not None else float('nan'))

        # MVA solver
        if MVA.supports(model):
            model.reset()
            solver_mva = MVA(model, verbose=False)
            solver_mva.avg_node_table()
            hr = cache_node.get_hit_ratio()
            mva_hit_ratio.append(hr if hr is not None else float('nan'))
        else:
            mva_hit_ratio.append(float('nan'))

        # NC solver
        if NC.supports(model):
            model.reset()
            solver_nc = NC(model, verbose=False)
            solver_nc.avg_node_table()
            hr = cache_node.get_hit_ratio()
            nc_hit_ratio.append(hr if hr is not None else float('nan'))
        else:
            nc_hit_ratio.append(float('nan'))

        print(f'{ReplacementStrategy.to_string(repl_strat[s])}: {ctmc_hit_ratio[s]:.8f}, {mva_hit_ratio[s]:.8f}, {nc_hit_ratio[s]:.8f}')
