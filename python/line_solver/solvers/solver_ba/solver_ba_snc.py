"""Stochastic network calculus UPPER bound for a feed-forward open network.

Native port of matlab/src/solvers/BA/solver_ba_snc_analyzer.m, cross-checked
against jline.solvers.ba.analyzers.Solver_ba_snc_analyzer. The api family in
:mod:`line_solver.api.snc` supplies the envelope algebra; this analyzer maps the
LINE model onto it, propagates envelopes hop by hop, and reads the bound back
per station and class. The result is valid for EVERY work-conserving scheduling
policy at every station.

UNITS ARE JOBS, NOT WORK. The arrival envelope counts jobs and the service
element is ``snc_srv_exp``, the counting process of an Exp(mu) server. That is
what lets a departure envelope from one station be the arrival envelope of the
next: a service-time work unit differs from station to station, a job does not.
On a single M/M/1 the resulting backlog bound decays as (lambda/mu)^n and the
delay bound as exp(-(mu-lambda)*d), both exact rates.

BOUND CONVENTION. R(i,r) is ``snc_mean_delay`` of the (arrival, service)
envelope pair at that station, i.e. the integral of the delay tail bound, so
each entry is a valid upper bound on its own. Q follows by Little's law from the
bounded R and the EXACT throughput T, and so does C. U is exact for the same
reason.

Reference: M. Fidler, A. Rizk (2015). A Guide to the Stochastic Network
Calculus. IEEE Communications Surveys and Tutorials 17(1), 92-105.
"""

import numpy as np

from ...api.snc import (
    snc_env_poisson, snc_env_map, snc_srv_exp,
    snc_leftover, snc_output, snc_mean_delay,
)

__all__ = ['snc_envelopes', 'snc_bound']

_TOL = 1e-10


def _sum_env(theta, parts):
    """Superposition of independent flows: the exponential forms multiply."""
    sigma = 0.0
    rho = 0.0
    for part in parts:
        s, r = part(theta)
        sigma += s
        rho += r
    return sigma, rho


def _make_sum(parts):
    return lambda theta: _sum_env(theta, parts)


def _make_poisson(lam):
    return lambda theta: snc_env_poisson(lam, theta)


def _make_map(D0, D1):
    return lambda theta: snc_env_map(D0, D1, theta)


def _make_leftover(mu, cross):
    def env(theta):
        sS, rS = snc_srv_exp(mu, theta)
        if not cross:
            return sS, rS
        sX, rX = _sum_env(theta, cross)
        return snc_leftover(sS, rS, sX, rX)
    return env


def _make_output(arv, srv):
    def env(theta):
        sA, rA = arv(theta)
        sS, rS = srv(theta)
        return snc_output(sA, rA, sS, rS, theta)
    return env


def _topo_order(adj):
    """Kahn's algorithm. Returns None when the graph has a cycle."""
    n = adj.shape[0]
    indeg = adj.sum(axis=0).astype(int)
    done = np.zeros(n, dtype=bool)
    order = []
    while True:
        cand = [i for i in range(n) if not done[i] and indeg[i] == 0]
        if not cand:
            break
        i = cand[0]
        order.append(i)
        done[i] = True
        indeg -= adj[i, :].astype(int)
        indeg[i] = 1  # keep it out of the candidate set
    if len(order) < n:
        return None
    return order


def snc_envelopes(sn):
    """Build the per-pair (arrival, service) envelopes of a feed-forward model.

    :param sn: the NetworkStruct
    :return: dict with ``arv``/``srv`` mapping ``(station, class) -> callable``,
        ``lam``/``mu`` mapping the same key to the exact rate and service rate,
        and ``M``/``K``
    """
    from ...api.sn.transforms import sn_rt_stations
    from ...api.sn.network_struct import NodeType
    from ...api.sn.proc_form import proc_to_map
    from ...constants import ProcessType, SchedStrategy

    M = int(sn.nstations)
    K = int(sn.nclasses)

    # ----- model gates -----
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    if np.any(np.isfinite(njobs)):
        raise ValueError("Method 'snc.upper' supports fully open networks only "
                         "(no closed classes).")
    sched_dict = sn.sched if sn.sched else {}
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    isSource = np.zeros(M, dtype=bool)
    for i in range(M):
        isSource[i] = (sn.nodetype[int(sn.stationToNode[i])] == NodeType.SOURCE)
    srcList = np.where(isSource)[0]
    qstat = np.where(~isSource)[0]
    if srcList.size == 0:
        raise ValueError("Method 'snc.upper' requires an open network with a Source station.")
    for i in qstat:
        if sched_dict.get(int(i), SchedStrategy.FCFS) == SchedStrategy.INF:
            raise ValueError("Method 'snc.upper' does not support delay (infinite-server) "
                             "stations: the service envelope is that of a single busy server.")
        if nservers[int(i)] > 1:
            raise ValueError("Method 'snc.upper' does not support multi-server stations.")

    # ----- station-space routing, with the Source absorbed into the injections -----
    rtst = np.asarray(sn_rt_stations(sn)[0], dtype=float)
    rates = np.asarray(sn.rates, dtype=float)

    pair_station, pair_class, pair_flat = [], [], []
    for i in qstat:
        for r in range(K):
            pair_station.append(int(i))
            pair_class.append(r)
            pair_flat.append(int(i) * K + r)
    npairs = len(pair_flat)
    pair_station = np.asarray(pair_station, dtype=int)
    pair_class = np.asarray(pair_class, dtype=int)
    pair_flat = np.asarray(pair_flat, dtype=int)

    # Injections are kept per (source, class) so that the exogenous process of
    # each stream is still identifiable once the pairs are known.
    cols = [(int(s), r0) for s in srcList for r0 in range(K)]
    inject = np.zeros((npairs, len(cols)))
    lambda0 = np.zeros(npairs)
    for col, (s, r0) in enumerate(cols):
        arr = rates[s, r0]
        if not np.isfinite(arr) or arr <= 0:
            continue
        inject[:, col] = arr * rtst[s * K + r0, pair_flat]
        lambda0 += inject[:, col]

    P = rtst[np.ix_(pair_flat, pair_flat)]

    # ----- restrict to the pairs that actually carry traffic -----
    lam_all = np.linalg.solve(np.eye(npairs) - P.T, lambda0)
    keep = np.where(lam_all > 1e-12 * max(1.0, float(lam_all.max())))[0]
    if keep.size == 0:
        raise ValueError("The model carries no open traffic.")
    lam = lam_all[keep]
    inject = inject[keep, :]
    P = P[np.ix_(keep, keep)]
    pair_station = pair_station[keep]
    pair_class = pair_class[keep]
    nk = keep.size

    mu = np.zeros(nk)
    for a in range(nk):
        i, r = int(pair_station[a]), int(pair_class[a])
        mu[a] = rates[i, r]
        if not np.isfinite(mu[a]) or mu[a] <= 0:
            raise ValueError("Station %d has no service rate for class %d but carries "
                             "its traffic." % (i + 1, r + 1))
        if sn.procid[i, r] != ProcessType.EXP:
            raise ValueError("Method 'snc.upper' requires exponential service: station %d "
                             "class %d is %s." % (i + 1, r + 1, sn.procid[i, r]))

    # ----- routing restrictions: no split downstream of the Source -----
    for a in range(nk):
        succ = np.where(P[a, :] > _TOL)[0]
        if succ.size > 1:
            raise ValueError(
                "Method 'snc.upper' requires deterministic routing downstream of the "
                "Source: station %d class %d splits its flow over %d destinations."
                % (pair_station[a] + 1, pair_class[a] + 1, succ.size))
        if succ.size == 1 and abs(P[a, succ[0]] - 1.0) > 1e-8:
            raise ValueError(
                "Method 'snc.upper' requires deterministic routing downstream of the "
                "Source: station %d class %d routes onward with probability %g."
                % (pair_station[a] + 1, pair_class[a] + 1, P[a, succ[0]]))

    # A Source that splits is exact only when its process is Poisson, since a
    # Bernoulli thinning of a Poisson stream is again Poisson.
    for col, (s, r0) in enumerate(cols):
        if inject[:, col].sum() <= 0:
            continue
        dest = np.where(inject[:, col] > _TOL)[0]
        if dest.size > 1 and sn.procid[s, r0] != ProcessType.EXP:
            raise ValueError(
                "Method 'snc.upper' can split only a Poisson Source: source %d class %d "
                "is %s and feeds %d stations."
                % (s + 1, r0 + 1, sn.procid[s, r0], dest.size))

    # ----- one service rate per station, and a feed-forward station graph -----
    stations_used = list(dict.fromkeys(int(x) for x in pair_station))
    for i in stations_used:
        r_at = mu[pair_station == i]
        if r_at.max() - r_at.min() > 1e-8 * max(1.0, float(r_at.max())):
            raise ValueError(
                "Method 'snc.upper' requires the classes sharing a station to have equal "
                "service rates: station %d carries rates in [%g, %g]."
                % (i + 1, r_at.min(), r_at.max()))

    ns = len(stations_used)
    st_idx = {i: a for a, i in enumerate(stations_used)}
    adj = np.zeros((ns, ns), dtype=bool)
    for a in range(nk):
        succ = np.where(P[a, :] > _TOL)[0]
        if succ.size == 0:
            continue
        adj[st_idx[int(pair_station[a])], st_idx[int(pair_station[succ[0]])]] = True
    order = _topo_order(adj)
    if order is None:
        raise ValueError(
            "Method 'snc.upper' requires a feed-forward network: the station graph has a "
            "cycle, so a station's cross traffic is not determined upstream of it.")

    # ----- envelope propagation, station by station in feed-forward order -----
    arvH = [None] * nk
    srvH = [None] * nk
    outH = [None] * nk
    for oi in order:
        i = stations_used[oi]
        here = [a for a in range(nk) if int(pair_station[a]) == i]
        for a in here:
            parts = []
            for col, (s, r0) in enumerate(cols):
                if inject[a, col] > _TOL:
                    if sn.procid[s, r0] == ProcessType.EXP:
                        parts.append(_make_poisson(float(inject[a, col])))
                    else:
                        D0, D1 = proc_to_map(sn.proc[s][r0])
                        parts.append(_make_map(np.atleast_2d(D0), np.atleast_2d(D1)))
            for b in range(nk):
                if P[b, a] > _TOL:
                    parts.append(outH[b])
            if not parts:
                raise ValueError("Station %d class %d carries traffic with no identifiable "
                                 "source." % (i + 1, pair_class[a] + 1))
            arvH[a] = _make_sum(parts)
        for a in here:
            cross = [arvH[b] for b in here if b != a]
            srvH[a] = _make_leftover(float(mu[a]), cross)
            outH[a] = _make_output(arvH[a], srvH[a])

    env = {'arv': {}, 'srv': {}, 'lam': {}, 'mu': {}, 'M': M, 'K': K,
           'srcList': srcList, 'rates': rates}
    for a in range(nk):
        key = (int(pair_station[a]), int(pair_class[a]))
        env['arv'][key] = arvH[a]
        env['srv'][key] = srvH[a]
        env['lam'][key] = float(lam[a])
        env['mu'][key] = float(mu[a])
    return env


def snc_bound(sn):
    """The (Q,U,R,T,C,X) block of the 'snc.upper' family.

    :param sn: the NetworkStruct
    :return: ``(Q, U, R, T, C, X)`` with R the mean-delay bound and Q its
        Little's-law image on the exact throughput
    """
    env = snc_envelopes(sn)
    M, K = env['M'], env['K']
    Q = np.zeros((M, K))
    U = np.zeros((M, K))
    R = np.zeros((M, K))
    T = np.zeros((M, K))
    C = np.zeros((1, K))
    X = np.zeros((1, K))

    for (i, r), arv in env['arv'].items():
        R[i, r] = snc_mean_delay(arv, env['srv'][(i, r)])[0]
        T[i, r] = env['lam'][(i, r)]
        U[i, r] = env['lam'][(i, r)] / env['mu'][(i, r)]

    rates = env['rates']
    for s in env['srcList']:
        for r in range(K):
            arr = rates[int(s), r]
            if np.isfinite(arr) and arr > 0:
                T[int(s), r] += arr
                X[0, r] += arr
    Q = T * R
    for r in range(K):
        if X[0, r] > 0:
            C[0, r] = float(np.sum(Q[:, r])) / X[0, r]
    return Q, U, R, T, C, X
