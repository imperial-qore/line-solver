"""
Exact MDD-stored solve of a single-class closed exponential network.

Solve a single-class closed exponential queueing network whose CTMC state
space (reachable occupancy vectors) is stored in a Multi-valued Decision
Diagram instead of an explicit state list. The reachable set is generated
with mdd_reachset and the generator matrix is assembled using the MDD's O(K)
state indexing (MDD.index), so no explicit (|S| x width) state matrix is ever
materialised during assembly -- the diagram is the store.

For single-class exponential stations the aggregated (occupancy) chain is
exact: the rate from n to n-e_i+e_j is mu_i * min(n_i, c_i) * P(i,j) for
n_i > 0, matching SolverCTMC on the same model, which makes this the live
exact oracle the mdd_mcd aggregation is validated against.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

import time
from typing import Dict, Optional

import numpy as np

from ..mc import ctmc_makeinfgen, ctmc_solve
from .mdd import MDD
from .reachset import mdd_reachset


def mdd_closedqn(mu, P, servers, N: int, options: Optional[Dict] = None) -> Dict[str, object]:
    """Exact solve of a closed exponential network over an MDD-stored state space.

    Parameters
    ----------
    mu : length-M per-station exponential service rates
    P : (M x M) Markovian routing matrix (row-stochastic, irreducible)
    servers : length-M number of servers per station (np.inf for a delay/IS)
    N : closed population
    options : dict with optional keys verbose (False), mdd (reuse an
        already-built reachable set as returned in out['mdd']), ctmcmethod
        (forwarded to ctmc_solve)

    Returns
    -------
    dict with keys mdd, Q, pi, states, QLen, U, X, stats and times (a dict of
    phase timings in seconds: reach, gen, solve, metrics; reach is 0 when
    options['mdd'] was supplied).
    """
    if options is None:
        options = {}
    verbose = bool(options.get('verbose', False))

    mu = np.asarray(mu, dtype=float).ravel()
    P = np.asarray(P, dtype=float)
    servers = np.asarray(servers, dtype=float).ravel()
    M = mu.size
    if P.shape != (M, M):
        raise ValueError('mdd_closedqn: routing matrix must be M x M')
    if servers.size != M:
        raise ValueError('mdd_closedqn: servers must have one entry per station')
    if N < 1:
        raise ValueError('mdd_closedqn: the closed population must be positive')
    domain = [N + 1] * M

    # events: a job completes at station i (rate mu_i * min(n_i, c_i)) and
    # routes to station j with probability P(i, j); self-routing i == j leaves
    # the occupancy vector unchanged and is skipped.
    ii, jj = np.nonzero(P)
    keep = ii != jj
    ii, jj = ii[keep], jj[keep]
    pij = P[ii, jj]
    E = ii.size

    # all jobs start at station 1; an irreducible routing chain makes every
    # composition of N over M stations reachable.
    init = [0] * M
    init[0] = N

    def nextfun(s):
        s = np.asarray(s, dtype=int)
        rows = []
        for a in range(E):
            if s[ii[a]] > 0:
                t = s.copy()
                t[ii[a]] -= 1
                t[jj[a]] += 1
                rows.append(t)
        return np.array(rows, dtype=int).reshape(len(rows), M)

    times = {}
    t0 = time.perf_counter()
    if options.get('mdd') is not None:
        mdd = options['mdd']  # caller already built the reachable set
        times['reach'] = 0.0
    else:
        mdd = mdd_reachset(domain, init, nextfun)
        times['reach'] = time.perf_counter() - t0

    t0 = time.perf_counter()
    n = mdd.cardinality()
    S = np.asarray(mdd.enumerate(), dtype=int)  # n x M states, MDD index order

    # assemble the generator directly from the MDD state indexing
    Q = np.zeros((n, n))
    for s in range(n):
        st = S[s]
        row = mdd.index(st)
        for a in range(E):
            i = ii[a]
            if st[i] > 0:
                busy = st[i] if np.isinf(servers[i]) else min(st[i], int(servers[i]))
                t = st.copy()
                t[i] -= 1
                t[jj[a]] += 1
                Q[row, mdd.index(t)] += mu[i] * busy * pij[a]
    Q = ctmc_makeinfgen(Q)
    times['gen'] = time.perf_counter() - t0

    t0 = time.perf_counter()
    p = np.asarray(ctmc_solve(Q, options.get('ctmcmethod'))).ravel()
    times['solve'] = time.perf_counter() - t0

    # performance metrics
    t0 = time.perf_counter()
    QLen = p @ S
    X = np.zeros(M)
    U = np.zeros(M)
    for i in range(M):
        busy = np.minimum(S[:, i], servers[i])  # busy servers in each state
        X[i] = mu[i] * (p @ busy)               # throughput = mean completion rate
        if np.isinf(servers[i]):
            U[i] = p @ S[:, i]                  # mean number busy (IS station)
        else:
            U[i] = (p @ busy) / servers[i]      # server utilisation
    times['metrics'] = time.perf_counter() - t0

    stats = mdd.stats()
    out = {
        'mdd': mdd,
        'times': times,
        'Q': Q,
        'pi': p,
        'states': S,
        'QLen': QLen,
        'U': U,
        'X': X,
        'stats': stats,
    }

    if verbose:
        from ..io.logging import line_printf
        line_printf('\nMDD-stored closed QN: %d stations, N=%d\n' % (M, N))
        line_printf('  reachable states |S| = %d\n' % stats['numStates'])
        line_printf('  MDD nodes            = %d  (%s per level)\n'
                    % (stats['numNodes'], ' '.join(str(v) for v in stats['nodesPerLevel'])))
        line_printf('  storage              = %d ints vs %d explicit (%.1fx)\n'
                    % (stats['mddInts'], stats['explicitInts'], stats['compression']))
    return out
