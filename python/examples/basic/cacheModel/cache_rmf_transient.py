"""
Cache RMF Transient Analysis

Demonstrates refined mean field (RMF) transient and steady-state
analysis of a multi-list cache with RANDOM(m) replacement.

Uses CacheRMF directly to compute:
  - Steady-state hit/miss probabilities with 1/N correction
  - Transient evolution of hit rates via coupled ODE system
"""

import numpy as np
from line_solver.solvers.solver_fld.methods.cache import CacheRMF


def cache_rmf_transient():
    # Cache parameters
    n = 10        # number of items
    m = [3, 2]    # list capacities (2-list cache)
    alpha = 0.8   # Zipf exponent

    # Zipf popularity distribution
    p = np.array([(i + 1) ** (-alpha) for i in range(n)])
    p = p / p.sum()

    print(f'Cache parameters: n={n}, m=[{", ".join(str(x) for x in m)}], Zipf({alpha:.1f})')

    # Build DDPP model
    model = CacheRMF(p, m)

    # Steady-state analysis with 1/N correction
    pi, V, (V_full, W) = model.meanFieldExpansionSteadyState(order=1)
    pi_refined = pi + V / n

    print('\nSteady-state results (refined mean field):')
    total_hit = 0.0
    for k in range(1, len(m) + 1):
        hr = model.hit_rate(pi_refined, k)
        total_hit += hr
        print(f'  Hit rate (list {k}): {hr:.6f}')
    miss_rate = model.hit_rate(pi_refined, 0)
    print(f'  Miss rate:         {miss_rate:.6f}')
    print(f'  Total hit prob:    {total_hit:.6f}')
    print(f'  Total miss prob:   {miss_rate:.6f}')

    # Transient analysis
    T, X, Vt, Wt = model.meanFieldExpansionTransient(time=50.0, n_points=200, order=1)

    print(f'\nTransient hit rates (refined, N={n}):')
    time_indices = [0, 20, 50, 100, 199]  # t=0, ~5, ~12.5, ~25, 50
    for idx in time_indices:
        xt = X[idx, :] + Vt[idx, :] / n
        hr = sum(model.hit_rate(xt, k) for k in range(1, len(m) + 1))
        print(f'  t={T[idx]:7.3f}: hit_rate={hr:.6f}')


if __name__ == "__main__":
    cache_rmf_transient()
