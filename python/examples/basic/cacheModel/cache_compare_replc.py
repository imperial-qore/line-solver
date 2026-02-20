"""
Compare Cache Replacement Strategies

This example demonstrates:
- Comparison of RR, FIFO, and LRU replacement strategies
- Open network with Zipf access pattern
- Hit ratio comparison across different strategies
"""

# Ensure native line_solver is used (not python-wrapper)
import sys
import os
_native_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
if _native_path not in sys.path:
    sys.path.insert(0, _native_path)

from line_solver import *
import numpy as np

if __name__ == "__main__":
    # Use STD verbosity to show warnings (like MATLAB behavior)
    GlobalConstants.set_verbose(VerboseLevel.STD)

    n = 5
    m = [2, 1]  # Two-level cache matching MATLAB
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
        # Extract scalar from array (first requesting class)
        if isinstance(hr, (float, np.floating, np.integer)):
            hr_val = float(hr)
        elif isinstance(hr, np.ndarray):
            if hr.ndim == 0:
                hr_val = float(hr)
            elif hr.ndim == 2 and hr.shape[0] > 0 and hr.shape[1] > 0:
                hr_val = float(hr[0][0])
            elif hr.ndim == 1 and len(hr) > 0:
                hr_val = float(hr[0])
            else:
                hr_val = float('nan')
        else:
            hr_val = float('nan')
        ctmc_hit_ratio.append(hr_val)

        # MVA solver
        try:
            model.reset()
            solver_mva = MVA(model, verbose=False)
            solver_mva.avg_node_table()
            hr = cache_node.get_hit_ratio()
            if isinstance(hr, (float, np.floating, np.integer)):
                hr_val = float(hr)
            elif isinstance(hr, np.ndarray):
                if hr.ndim == 0:
                    hr_val = float(hr)
                elif hr.ndim == 2 and hr.shape[0] > 0 and hr.shape[1] > 0:
                    hr_val = float(hr[0][0])
                elif hr.ndim == 1 and len(hr) > 0:
                    hr_val = float(hr[0])
                else:
                    hr_val = float('nan')
            else:
                hr_val = float('nan')
            mva_hit_ratio.append(hr_val)
        except:
            mva_hit_ratio.append(float('nan'))

        # NC solver - check if model is supported first (prints warnings for unsupported features)
        if NC.supports(model):
            try:
                model.reset()
                solver_nc = NC(model, verbose=False)
                solver_nc.avg_node_table()
                hr = cache_node.get_hit_ratio()
                if isinstance(hr, (float, np.floating, np.integer)):
                    hr_val = float(hr)
                elif isinstance(hr, np.ndarray):
                    if hr.ndim == 0:
                        hr_val = float(hr)
                    elif hr.ndim == 2 and hr.shape[0] > 0 and hr.shape[1] > 0:
                        hr_val = float(hr[0][0])
                    elif hr.ndim == 1 and len(hr) > 0:
                        hr_val = float(hr[0])
                    else:
                        hr_val = float('nan')
                else:
                    hr_val = float('nan')
                nc_hit_ratio.append(hr_val)
            except:
                nc_hit_ratio.append(float('nan'))
        else:
            nc_hit_ratio.append(float('nan'))

        print(f'{ReplacementStrategy.to_string(repl_strat[s])}: {ctmc_hit_ratio[s]:.8f}, {mva_hit_ratio[s]:.8f}, {nc_hit_ratio[s]:.8f}')
