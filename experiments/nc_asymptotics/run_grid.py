"""Compare the six asymptotic normalizing-constant methods on load-independent
closed models over a grid of stations (M), classes (R) and populations (N).

kt / le / ble reach the load-independent dispatcher pfqn_nc directly.
nrl / nrp / nre are load-dependent evaluators, so they are driven the way
SolverNC drives them on a load-independent model: pfqn_ncld with mu == 1.

Reference is pfqn_ca (convolution), which is exact for this model class.
"""
import itertools
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, '/data/gcasale/line-dev.git/python')

LI_METHODS = ['kt', 'le', 'ble']
LD_METHODS = ['nrl', 'nrp', 'nre']
METHODS = LI_METHODS + LD_METHODS


def split_pop(Ntot, R):
    """Population vector for R classes, as even as the total allows."""
    base = Ntot // R
    N = np.full(R, base, dtype=float)
    N[:Ntot - int(N.sum())] += 1
    return N


def eval_lG(method, L, N, Z):
    """log G by one method; returns (lG, seconds) or (nan, seconds) on failure."""
    from line_solver.api.pfqn.nc import pfqn_nc
    from line_solver.api.pfqn.ncld import pfqn_ncld
    Nt = int(round(float(np.sum(N))))
    t0 = time.perf_counter()
    try:
        if method in LI_METHODS:
            _, lG = pfqn_nc(L, N, Z, method=method)
        else:
            mu = np.ones((L.shape[0], max(Nt, 1)))
            lG = pfqn_ncld(L, N, Z, mu, {'method': method}).lG
        lG = float(np.real(lG))
    except Exception as exc:                       # noqa: BLE001 - recorded, not swallowed
        return float('nan'), time.perf_counter() - t0, type(exc).__name__
    return lG, time.perf_counter() - t0, ''


def make_demands(rng, M, R, regime):
    if regime == 'balanced':
        # near-identical demands: the degenerate case for saddle-point methods
        return np.ones((M, R)) * (1.0 + 0.01 * rng.standard_normal((M, R)))
    return rng.uniform(0.1, 1.0, (M, R))


def run_case(case):
    from line_solver.api.pfqn.nc import pfqn_ca
    block, M, R, Ntot, rep, regime, ztime, want_tput = case
    # deterministic seed: PYTHONHASHSEED randomises str hashing across workers
    seed = int(np.ravel_multi_index(
        (ord(block) - 65, M, R, Ntot, rep), (3, 21, 4, 101, 5)))
    rng = np.random.default_rng(seed)
    L = make_demands(rng, M, R, regime)
    N = split_pop(Ntot, R)
    Z = np.full(R, float(ztime))

    t0 = time.perf_counter()
    _, lG_exact = pfqn_ca(L, N, Z)
    t_exact = time.perf_counter() - t0
    lG_exact = float(np.real(lG_exact))

    rows = []
    # exact and approximate log-constants at N - e_r, for the throughputs
    exact_dn, approx_dn = {}, {m: {} for m in METHODS}
    if want_tput:
        for r in range(R):
            Nr = N.copy()
            Nr[r] -= 1
            _, lg = pfqn_ca(L, Nr, Z)
            exact_dn[r] = float(np.real(lg))

    for m in METHODS:
        lG, secs, err = eval_lG(m, L, N, Z)
        tput_err = float('nan')
        if want_tput and np.isfinite(lG):
            errs = []
            for r in range(R):
                Nr = N.copy()
                Nr[r] -= 1
                lg_r, s_r, e_r = eval_lG(m, L, Nr, Z)
                secs += s_r
                if not np.isfinite(lg_r):
                    errs = [float('nan')]
                    break
                x_hat = np.exp(lg_r - lG)
                x_ex = np.exp(exact_dn[r] - lG_exact)
                errs.append(abs(x_hat - x_ex) / x_ex if x_ex > 0 else float('nan'))
            tput_err = float(np.max(errs)) if errs else float('nan')
        rows.append(dict(block=block, regime=regime, M=M, R=R, Ntot=Ntot, rep=rep,
                         Z=ztime, method=m, lG_exact=lG_exact, lG=lG,
                         dlG=lG - lG_exact, abs_dlG=abs(lG - lG_exact),
                         tput_err=tput_err, secs=secs, t_exact=t_exact, error=err))
    return rows


def build_cases():
    cases = []
    # Block A: the main sweep -- pure queueing network, random demands
    for M, R, Ntot in itertools.product([2, 3, 5, 10, 20], [1, 2, 3],
                                        [2, 5, 10, 20, 50, 100]):
        if Ntot < R:
            continue
        for rep in range(5):
            cases.append(('A', M, R, Ntot, rep, 'random', 0.0, Ntot <= 50))
    # Block B: same methods with a think time (delay station present)
    for M, R, Ntot in itertools.product([2, 5, 20], [1, 2, 3], [5, 20, 50]):
        for rep in range(3):
            cases.append(('B', M, R, Ntot, rep, 'random', 5.0, True))
    # Block C: near-balanced demands, the degenerate saddle point
    for M, R, Ntot in itertools.product([2, 5, 20], [1, 2, 3], [5, 20, 50]):
        for rep in range(3):
            cases.append(('C', M, R, Ntot, rep, 'balanced', 0.0, True))
    # cheapest first is wrong for a pool: put the long poles at the front
    cases.sort(key=lambda c: -(c[1] * c[3] ** 2 * c[2]))
    return cases


def main():
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'grid_results.jsonl')
    cases = build_cases()
    print(f'{len(cases)} models, {len(cases) * len(METHODS)} method evaluations', flush=True)
    done = 0
    with open(out, 'w') as fh, Pool(8) as pool:
        for rows in pool.imap_unordered(run_case, cases, chunksize=1):
            for row in rows:
                fh.write(json.dumps(row) + '\n')
            fh.flush()
            done += 1
            if done % 25 == 0:
                print(f'{done}/{len(cases)}', flush=True)
    print(f'wrote {out}', flush=True)


if __name__ == '__main__':
    main()
