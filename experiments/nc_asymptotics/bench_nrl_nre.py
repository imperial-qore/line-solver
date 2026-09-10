"""Runtime of nre against nrl on load-independent closed models.

Wall clock is useless on a shared host, so cost is CPU time (time.process_time)
for ONE lG evaluation, min over repeats. Calls into the load-dependent kernel
pfqn_gld / pfqn_gldsingle are counted and timed separately, because that is
where both methods spend almost everything and it is what separates them.
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS'):
    os.environ[_v] = '1'

import json
import sys
import time

import numpy as np

sys.path.insert(0, '/data/gcasale/line-dev.git/python')
from line_solver.api.pfqn import ncld as ncld_mod                # noqa: E402
from line_solver.api.pfqn.ncld import pfqn_ncld                  # noqa: E402

COUNT = {'gld': 0, 'gldsingle': 0, 'kernel_s': 0.0}
DEPTH = [0]
_gld, _gldsingle = ncld_mod.pfqn_gld, ncld_mod.pfqn_gldsingle


def _wrap(fn, key):
    # pfqn_gld calls pfqn_gldsingle, so only the OUTERMOST entry is timed
    def inner(*a, **k):
        COUNT[key] += 1
        if DEPTH[0]:
            return fn(*a, **k)
        DEPTH[0] = 1
        t = time.process_time()
        try:
            return fn(*a, **k)
        finally:
            DEPTH[0] = 0
            COUNT['kernel_s'] += time.process_time() - t
    return inner


ncld_mod.pfqn_gld = _wrap(_gld, 'gld')
ncld_mod.pfqn_gldsingle = _wrap(_gldsingle, 'gldsingle')


def split(Ntot, R):
    N = np.full(R, Ntot // R, dtype=float)
    N[:Ntot - int(N.sum())] += 1
    return N


def measure(method, M, R, Ntot, repeats):
    rng = np.random.default_rng(1000 * M + 10 * R + Ntot)
    L = rng.uniform(0.1, 1.0, (M, R))
    N = split(Ntot, R)
    Z = np.zeros(R)
    mu = np.ones((M, Ntot))
    best, prof = None, None
    for _ in range(repeats):
        COUNT.update(gld=0, gldsingle=0, kernel_s=0.0)
        DEPTH[0] = 0
        t = time.process_time()
        pfqn_ncld(L, N, Z, mu, {'method': method})
        cpu = time.process_time() - t
        if best is None or cpu < best:
            best, prof = cpu, dict(COUNT)
    return dict(method=method, M=M, R=R, Ntot=Ntot, cpu=best,
                gld=prof['gld'], gldsingle=prof['gldsingle'],
                kernel_s=prof['kernel_s'])


def cases():
    out = []
    for Ntot in [5, 10, 20, 50, 100, 200]:                # scaling in population
        out.append(('N', 5, 2, Ntot))
    for M in [2, 5, 10, 20, 50, 100]:                     # scaling in stations
        out.append(('M', M, 2, 50))
    for R in [1, 2, 3, 4]:                                # scaling in classes
        out.append(('R', 5, R, 48))
    return out


def main():
    res = []
    for sweep, M, R, Ntot in cases():
        for method in ('nrl', 'nre'):
            reps = 3 if M * Ntot <= 500 else (2 if M * Ntot <= 2500 else 1)
            r = measure(method, M, R, Ntot, reps)
            r['sweep'] = sweep
            r['repeats'] = reps
            res.append(r)
            print(f"{sweep} M={M:3d} R={R} N={Ntot:3d} {method}  "
                  f"cpu={r['cpu']:9.3f}s  kernel={r['kernel_s']:8.3f}s "
                  f"({100 * r['kernel_s'] / r['cpu'] if r['cpu'] else 0:4.1f}%)  "
                  f"gld={r['gld']:6d} gldsingle={r['gldsingle']:6d}", flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, 'bench_nrl_nre.json'), 'w') as fh:
        json.dump(res, fh, indent=1)
    print('wrote bench_nrl_nre.json')


if __name__ == '__main__':
    main()
