"""
Passage-time (response time) distribution of the open MMAP[K]/PH[K]/1 queue.

Native twin of MATLAB ``solver_mam_passage_time.m`` (JAR
``Solver_mam_passage_time.java``, C++ ``solver_mam_passage_time.h``): a two
station open model (Source + one queue) is analyzed exactly, FCFS/HOL through
the BuTools MMAPPH1FCFS sojourn-time PH law, PS through the Masuyama-Takine
MAP/M/1-PS sojourn distribution. Any other topology gets a warning and no
result, exactly as the reference.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import warnings

import numpy as np

from ...api.mam.map_analysis import (
    map_cdf,
    map_lambda,
    map_mean,
    map_pie,
    map_scale,
    map_var,
)
from ...api.mam.mmap_ops import mmap_super
from ...api.mam.mapm1ps import map_m1ps_cdf_respt
from ...constants import GlobalConstants
from .algorithms.transient_qbd import _read_map


def _sched_name(sn, ist):
    """The scheduling strategy of station ist as an upper-case name string."""
    sched = sn.sched
    if isinstance(sched, dict):
        s = sched.get(ist, None)
    else:
        s = sched[ist] if ist < len(sched) else None
    if s is None:
        return ''
    return s.name if hasattr(s, 'name') else str(s).upper()


def _proc_map(sn, ist, r):
    """(D0, D1) of the service/arrival process of (station, class), from the
    native sn.proc in any of its storage layouts."""
    proc = sn.proc
    entry = None
    if isinstance(proc, list):
        row = proc[ist] if ist < len(proc) else None
        if isinstance(row, (list, tuple)):
            entry = row[r] if r < len(row) else None
        elif isinstance(row, dict):
            entry = row.get(r, None)
        else:
            entry = row
    elif isinstance(proc, dict):
        entry = proc.get((ist, r), None)
        if entry is None:
            row = proc.get(ist, None)
            if isinstance(row, (list, tuple)) and r < len(row):
                entry = row[r]
            elif isinstance(row, dict):
                entry = row.get(r, None)
    if entry is None:
        return None, None
    return _read_map(entry)


def _priority_passage_time(sn, A, pies, S, idx_arv, idx_q, sched_q, prio, n_cdf_pts):
    """Sojourn CDF of the MMAP[K]/PH[K]/1 priority queue, all priorities
    distinct: MMAPPH1PRPR (preemptive resume, FCFSPRPRIO) or MMAPPH1NPPR
    (non-preemptive, HOL).

    Neither vendored analyzer exports 'stDistrPH', so the law is TABULATED: one
    solve for two sojourn moments sizes a shared grid (mean + 10 sigma over the
    classes), a second tabulates the CDF on it. BUTools wants D1=lowest ..
    DK=highest priority where LINE's convention is lower value = higher
    priority, so the class order flips on the way in and the outputs map back.
    """
    from line_solver.lib.thirdparty.butools.queues import MMAPPH1NPPR, MMAPPH1PRPR

    K = int(sn.nclasses)
    # ascending priority VALUE = LINE-highest first; BUTools wants lowest first
    order = np.argsort(prio, kind='stable')[::-1]
    D_marked = [np.matrix(np.asarray(A[0], dtype=float))]
    for c in order:
        D_marked.append(np.matrix(np.asarray(A[2 + int(c)], dtype=float)))
    sigma_in = [np.matrix(np.asarray(pies[int(c)], dtype=float).reshape(1, -1)) for c in order]
    S_in = [np.matrix(S[int(c)]) for c in order]
    solver = MMAPPH1PRPR if sched_q == 'FCFSPRPRIO' else MMAPPH1NPPR

    moms = solver(D_marked, sigma_in, S_in, 'stMoms', 2)
    xmax = 0.0
    for b in range(K):
        mk = np.asarray(moms[b], dtype=float).ravel()
        m1 = float(mk[0])
        sig = np.sqrt(max(float(mk[1]) - m1 ** 2, 0.0))
        xmax = max(xmax, m1 + 10.0 * sig)
    X = np.linspace(0.0, xmax, n_cdf_pts)
    # The Erlangization inside stDistr divides by t, so t = 0 is not evaluable;
    # a sojourn time is strictly positive, so F(0) = 0 exactly and the grid's
    # first point is prepended rather than asked for. A plain list, because the
    # vendored option scan compares every argument against the option names and
    # a numpy array makes that comparison ambiguous.
    dist = solver(D_marked, sigma_in, S_in, 'stDistr', [float(x) for x in X[1:]])

    RD = []
    for b in range(K):
        korig = int(order[b])
        F = np.concatenate([[0.0], np.asarray(dist[b], dtype=float).ravel()])
        RD.append({
            'station': idx_q + 1,
            'class': korig + 1,
            't': X.copy(),
            'p': F,
        })
    RD.sort(key=lambda e: e['class'])
    return RD


def solver_mam_passage_time(sn, options=None):
    """Response time CDF of the open MMAP[K]/PH[K]/1 queue, flat native contract.

    Returns:
        List of dicts with 'station', 'class' (1-based), 't', 'p', one per
        class at the queueing station -- the Source, as in the reference, gets
        no law. An unsupported topology returns [] after a warning.
    """
    M = int(sn.nstations)
    K = int(sn.nclasses)
    njobs = np.asarray(sn.njobs, dtype=float).flatten()

    # The reference's GLOBAL default (MATLAB SolverOptions.m sets
    # options.config.num_cdf_pts = 200, and the C++ MamOptions carries 200);
    # the 100-point fallback inside solver_mam_passage_time.m is dead code
    # there because the option always arrives set
    n_cdf_pts = 200
    config = getattr(options, 'config', None) if options is not None else None
    if isinstance(config, dict):
        n_cdf_pts = int(config.get('num_cdf_pts', n_cdf_pts) or n_cdf_pts)
    elif config is not None and getattr(config, 'num_cdf_pts', None):
        n_cdf_pts = int(config.num_cdf_pts)

    if not (M == 2 and np.all(np.isinf(njobs))):
        warnings.warn('This model is not supported by SolverMAM yet. '
                      'Returning with no result.')
        return []

    nservers = np.asarray(sn.nservers, dtype=float).flatten()

    A = None
    idx_arv = None
    idx_q = None
    is_ps = False
    pies = [None] * K
    S = [None] * K

    for ist in range(M):
        s = _sched_name(sn, ist)
        if s == 'EXT':
            for k in range(K):
                D0k, D1k = _proc_map(sn, ist, k)
                if D0k is None or D1k is None:
                    warnings.warn('This model is not supported by SolverMAM yet. '
                                  'Returning with no result.')
                    return []
                mk = [np.asarray(D0k, dtype=float), np.asarray(D1k, dtype=float),
                      np.asarray(D1k, dtype=float)]
                A = mk if A is None else mmap_super(A, mk)
            idx_arv = ist
        elif s in ('FCFS', 'HOL', 'FCFSPRPRIO'):
            for k in range(K):
                D0k, D1k = _proc_map(sn, ist, k)
                if D0k is None or D1k is None:
                    warnings.warn('This model is not supported by SolverMAM yet. '
                                  'Returning with no result.')
                    return []
                # A c-server station is folded into a c-times-faster server,
                # as the reference scales it
                D0k, D1k = map_scale(np.asarray(D0k, dtype=float),
                                     np.asarray(D1k, dtype=float),
                                     map_mean(np.asarray(D0k, dtype=float),
                                              np.asarray(D1k, dtype=float)) / nservers[ist])
                pies[k] = map_pie(D0k, D1k)
                S[k] = D0k
            idx_q = ist
            is_ps = False
        elif s == 'PS':
            for k in range(K):
                D0k, D1k = _proc_map(sn, ist, k)
                if D0k is None or D1k is None:
                    warnings.warn('This model is not supported by SolverMAM yet. '
                                  'Returning with no result.')
                    return []
                S[k] = np.asarray(D0k, dtype=float)
                pies[k] = map_pie(np.asarray(D0k, dtype=float), np.asarray(D1k, dtype=float))
            idx_q = ist
            is_ps = True
        else:
            raise RuntimeError('Unsupported scheduling strategy')

    if A is None or idx_q is None:
        warnings.warn('This model is not supported by SolverMAM yet. '
                      'Returning with no result.')
        return []

    sched_q = _sched_name(sn, idx_q)
    classprio = getattr(sn, 'classprio', None)
    if (sched_q in ('HOL', 'FCFSPRPRIO') and classprio is not None
            and np.asarray(classprio).size >= K > 1):
        # Priorities select the law only under a priority DISCIPLINE: a plain
        # FCFS or PS queue serves in arrival or processor order whatever the
        # classprio column says, exactly as the reference gates these analyzers
        prio = np.asarray(classprio, dtype=float).flatten()[:K]
        if np.any(prio != prio[0]):
            if len(np.unique(prio)) != K:
                raise RuntimeError('SolverMAM requires either identical priorities '
                                   'or all distinct priorities')
            return _priority_passage_time(sn, A, pies, S, idx_arv, idx_q,
                                          sched_q, prio, n_cdf_pts)

    RD = []

    if is_ps:
        # MAP/M/1-PS: Masuyama-Takine sojourn time distribution, which needs
        # exponential service
        for k in range(K):
            if S[k].shape[0] != 1:
                raise RuntimeError('PS queue requires exponential (Markovian) service times')
        mu_vec = np.array([-S[k][0, 0] for k in range(K)])
        if np.any(np.abs(mu_vec - mu_vec[0]) > GlobalConstants.FineTol):
            raise RuntimeError('Multi-class PS currently requires identical service rates')
        mu = float(mu_vec[0])

        C_map = np.asarray(A[0], dtype=float)
        if K == 1:
            D_map = np.asarray(A[1], dtype=float)
        else:
            # Aggregate all class arrival streams into a single MAP
            D_map = np.sum(np.stack([np.asarray(A[j], dtype=float)
                                     for j in range(2, len(A))]), axis=0)

        lam = map_lambda(C_map, D_map)
        rho = lam / mu
        approx_mean = 1.0 / (mu * (1.0 - rho))
        x_vals = np.linspace(0.0, approx_mean * 10.0, n_cdf_pts)
        W_bar = map_m1ps_cdf_respt(C_map, D_map, mu, x_vals)
        F = 1.0 - np.asarray(W_bar, dtype=float).ravel()
        for k in range(K):
            RD.append({
                'station': idx_q + 1,
                'class': k + 1,
                't': x_vals.copy(),
                'p': F.copy(),
            })
        return RD

    # FCFS/HOL: MMAPPH1FCFS sojourn time as a PH law per class
    from line_solver.lib.thirdparty.butools.queues import MMAPPH1FCFS

    D_marked = [np.matrix(np.asarray(A[0], dtype=float))]
    for j in range(2, len(A)):
        D_marked.append(np.matrix(np.asarray(A[j], dtype=float)))
    sigma_in = [np.matrix(np.asarray(pies[k], dtype=float).reshape(1, -1)) for k in range(K)]
    S_in = [np.matrix(S[k]) for k in range(K)]

    ret = MMAPPH1FCFS(D_marked, sigma_in, S_in, 'stDistrPH')

    for k in range(K):
        alpha = np.asarray(ret[2 * k], dtype=float).reshape(1, -1)
        D0 = np.asarray(ret[2 * k + 1], dtype=float)
        n_ph = D0.shape[0]
        # The PH law (alpha, D0) as a MAP: D1 = (-D0) 1 alpha
        D1 = (-D0) @ np.ones((n_ph, 1)) @ alpha
        mean_k = map_mean(D0, D1)
        # The PH law is defective on states alpha never starts in; the pie
        # weighting inside map_mean handles it, but the SUPPORT of the CDF grid
        # still has to reach the tail, expanding until it holds 1 - FineTol
        sigma_k = np.sqrt(max(map_var(D0, D1), 0.0))
        n = 5
        while map_cdf(D0, D1, np.array([mean_k + n * sigma_k]))[0] < 1 - GlobalConstants.FineTol:
            n += 1
        X = np.linspace(0.0, mean_k + n * sigma_k, n_cdf_pts)
        F = np.asarray(map_cdf(D0, D1, X), dtype=float).ravel()
        RD.append({
            'station': idx_q + 1,
            'class': k + 1,
            't': X.copy(),
            'p': F.copy(),
        })
    return RD
