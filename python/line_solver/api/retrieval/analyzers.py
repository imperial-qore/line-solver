"""NC and MVA analyzers for delayed-hit (retrieval-system) caches (native Python).

Mirror of matlab/src/solvers/{NC,MVA}/solver_{nc,mva}_retrieval_analyzer.m and
java/.../Solver_{nc,mva}_retrieval_analyzer.java.
"""
import time
import warnings
from types import SimpleNamespace

import numpy as np

from ..sn import NodeType
from .cache_retrieval_inputs import cache_retrieval_inputs
from .retrieval_nc import retrieval_nc
from .retrieval_metrics import retrieval_metrics
from .retrieval_rayint import retrieval_rayint
from .retrieval_fpi import retrieval_fpi
from .retrieval_fpi_latency import retrieval_fpi_latency


def _station_metrics(sn, inp, pdh, pmiss, QN, UN, RN, TN):
    """QN=UN=phi_s, TN=fetch throughput, RN=Little.

    ``pmiss`` is the per-item miss probability pi_{i,0}: only a MISS is fetched
    through a retrieval station, so it weights the throughput.
    """
    n = len(inp['lambda'])
    S = len(inp['station_type'])
    ps_idx = [s for s in range(S) if inp['station_type'][s] != "IS"]
    for s in range(S):
        ist = inp['queue_station'][s]
        if inp['station_type'][s] == "IS":
            phi_s = float(pdh[0, :].sum())
        else:
            p = ps_idx.index(s)
            phi_s = float(pdh[1 + p, :].sum()) if (1 + p) < pdh.shape[0] else 0.0
        # MATLAB solver_nc_retrieval_analyzer.m:
        #   tput_s += sourceRate(jobinClass) * (lambda(i)/sum(lambda)) * pi0(i) * vis(s)
        # Only a MISS is fetched through a retrieval station, so the fetch rate
        # is the request rate times the ACCESS-WEIGHTED MISS probability, not the
        # request rate itself. Omitting pi0 (and the sourceRate and the lambda
        # normalisation) counted every request as a fetch, which on a model whose
        # access probabilities already sum to one returns exactly the source
        # rate: 1.0 where MATLAB gives 0.462797.
        tot_lambda = float(np.sum(inp['lambda']))
        tput = 0.0
        for i in range(n):
            R = inp['R'][i]
            a = R[0, 1:S + 1]
            P = R[1:S + 1, 1:S + 1]
            visits = np.linalg.solve(np.eye(S) - P.T, a)
            w = (inp['lambda'][i] / tot_lambda) if tot_lambda > 0 else 0.0
            miss_i = float(pmiss[i]) if pmiss is not None and i < len(pmiss) else 1.0
            tput += inp['source_rate'] * w * miss_i * visits[s]
        QN[ist, inp['jobin_class']] = phi_s
        UN[ist, inp['jobin_class']] = phi_s
        TN[ist, inp['jobin_class']] = tput
        if tput > 0:
            RN[ist, inp['jobin_class']] = phi_s / tput


def _agg(inp, pmiss, phit, pdh):
    n = len(inp['lambda'])
    tot = inp['lambda'].sum()
    hit = miss = delayed = 0.0
    for i in range(n):
        w = inp['lambda'][i] / tot
        hit += w * phit[:, i].sum()
        miss += w * pmiss[i]
        delayed += w * pdh[:, i].sum()
    return hit, miss, delayed


def _result(sn, inp, hit, miss, delayed, pdh, lG, latency_val, method, phit=None, pmiss=None):
    K = sn.nclasses
    M = sn.nstations
    XN = np.zeros((1, K))
    QN = np.zeros((M, K)); UN = np.zeros((M, K)); RN = np.zeros((M, K)); TN = np.zeros((M, K))

    ch = sn.nodeparam[inp['cache_idx']]
    hc = int(ch.hitclass[inp['jobin_class']])
    mc = int(ch.missclass[inp['jobin_class']])
    if hc >= 0:
        XN[0, hc] = inp['source_rate'] * (hit + delayed)
    if mc >= 0:
        XN[0, mc] = inp['source_rate'] * miss

    # source throughput (station-indexed)
    for i, nt in enumerate(sn.nodetype):
        if nt == NodeType.SOURCE or int(nt) == int(NodeType.SOURCE):
            sist = int(sn.nodeToStation[i])
            for k in range(K):
                rt = sn.rates[sist, k]
                if not np.isnan(rt):
                    TN[sist, k] = rt
            break

    _station_metrics(sn, inp, pdh, pmiss, QN, UN, RN, TN)

    hitprob = np.full((1, K), np.nan)
    missprob = np.full((1, K), np.nan)
    delayedprob = np.full((1, K), np.nan)
    expected_latency = np.full((1, K), np.nan)
    hitprob[0, inp['jobin_class']] = hit
    missprob[0, inp['jobin_class']] = miss
    delayedprob[0, inp['jobin_class']] = delayed
    expected_latency[0, inp['jobin_class']] = latency_val
    # per-list (per-level) hit fractions for the read class: phit is (h x n);
    # access-weighted over items gives the per-list hit probability.
    h = len(inp['m'])
    hitproblist = np.full((K, h), np.nan)
    n = len(inp['lambda'])
    itemprob = None
    if phit is not None:
        tot = float(np.sum(inp['lambda']))
        for l in range(h):
            acc = 0.0
            for i in range(n):
                acc += (inp['lambda'][i] / tot) * phit[l, i]
            hitproblist[inp['jobin_class'], l] = acc
        # per-item occupancy [n x (h+1)]: col 0 = miss, cols 1.. = per-list
        itemprob = np.zeros((n, h + 1))
        itemprob[:, 0] = pmiss
        for l in range(h):
            itemprob[:, l + 1] = phit[l, :]
    return SimpleNamespace(QN=QN, UN=UN, RN=RN, TN=TN, XN=XN, lG=lG,
                           method=method, runtime=0.0, it=1, pij=None,
                           hitprob=hitprob, missprob=missprob,
                           delayedprob=delayedprob, hitproblist=hitproblist,
                           itemprob=itemprob,
                           expected_latency=expected_latency)


def solver_nc_retrieval_analyzer(sn, options=None):
    """Exact delayed-hit analysis via retrieval_nc + retrieval_metrics.

    Those recurrences are exponential in the item count, so ``options.method =
    'rayint'`` selects instead the ray (WKB) approximation of ``retrieval_rayint``,
    which is polynomial.  It applies only when every fetch station is
    infinite-server, where the delayed-hit constant factorizes exactly as
    ``prod_k D_k`` times the plain cache constant with access factors
    ``gamma_{k,j}/D_k``; elsewhere it warns and falls back to the exact path.
    """
    t0 = time.time()
    inp = cache_retrieval_inputs(sn)
    r = inp['eta'].shape[1] - 1

    method = str(getattr(options, 'method', '') or '').lower()
    useray = method in ('rayint', 'ray')
    if useray:
        # The ray expansion needs the delayed-hit constant to factorize.  It does,
        # EXACTLY, when every fetch station is infinite-server: dividing the
        # retrieval_nc recurrence by prod_k D_k with D_k = 1 + lambda_k eta_{0,k}
        # collapses it onto cache_erec with theta_{k,j} = gamma_{k,j}/D_k, so the
        # delayed-hit cache IS a plain cache with fetch-inflated access factors.
        # A queueing (PS) fetch station breaks this: the (v_s+1) multiplicity ties
        # E(0,m) to the whole moment tower E(1_s,m), E(2_s,m), ..., and replacing it
        # by the retrieval_fpi mean field overestimates E by 13%/140%/830% at
        # n=6/8/10 (measured), growing with n.  Refuse rather than return a
        # confident wrong number.
        reason = None
        nitems = len(inp['lambda'])
        if r > 0 and np.any(inp['eta'][:, 1:] != 0):
            reason = "the retrieval system has a queueing (non infinite-server) fetch station"
        elif float(np.sum(inp['m'])) >= nitems:
            reason = "the cache is full (sum(m) >= n), where the saddle point escapes to infinity"
        if reason is not None:
            warnings.warn("solver_nc_retrieval_analyzer: method 'rayint' does not apply because "
                          "%s; falling back to the exact recurrences." % reason, RuntimeWarning)
            useray = False

    if useray:
        # --- ray (WKB) approximation, infinite-server fetch ---
        lam = np.asarray(inp['lambda'], dtype=float).ravel()
        D = 1.0 + lam * np.asarray(inp['eta'], dtype=float)[:, 0]
        theta = np.asarray(inp['gamma'], dtype=float) / D[:, None]
        _, lGcache, rayout = retrieval_rayint(theta, inp['m'])
        lG = float(np.sum(np.log(D)) + lGcache)

        # Same saddle as the constant, so the ratios are consistent with lG:
        # pi_{i,j} = theta_{i,j} xi_j / (1 + sum_l theta_{i,l} xi_l), and the
        # out-of-cache mass 1 - sum_j pi_{i,j} splits between a true miss (weight 1)
        # and an outstanding fetch (weight lambda_i eta_{0,i}) in proportion 1:D_i-1.
        txi = theta * rayout.xi
        phit = (txi / (1.0 + txi.sum(axis=1))[:, None]).T
        pmiss = (1.0 - phit.sum(axis=0)) / D
        pdh = (lam * np.asarray(inp['eta'], dtype=float)[:, 0] * pmiss)[None, :]
        used = "rayint"
    else:
        E = retrieval_nc(np.zeros(r), inp['m'], inp['lambda'], inp['eta'], inp['gamma'])
        lG = float(np.log(E))
        pmiss, phit, pdh = retrieval_metrics(inp['m'], inp['lambda'], inp['eta'], inp['gamma'])
        used = "exact"

    hit, miss, delayed = _agg(inp, pmiss, phit, pdh)
    res = _result(sn, inp, hit, miss, delayed, pdh, lG, np.nan, used, phit=phit, pmiss=pmiss)
    res.runtime = time.time() - t0
    return res


def solver_mva_retrieval_analyzer(sn, options=None):
    t0 = time.time()
    inp = cache_retrieval_inputs(sn)
    pmiss, phit, pdh = retrieval_fpi(inp['m'], inp['lambda'], inp['eta'], inp['gamma'])
    Z, _, _, _ = retrieval_fpi_latency(inp['m'], inp['lambda'], inp['gamma'],
                                       inp['alpha'], inp['T'], inp['R'], inp['station_type'])
    hit, miss, delayed = _agg(inp, pmiss, phit, pdh)
    res = _result(sn, inp, hit, miss, delayed, pdh, np.nan, Z, "fpi", phit=phit, pmiss=pmiss)
    res.runtime = time.time() - t0
    return res


def _has_source(sn):
    for nt in sn.nodetype:
        if nt == NodeType.SOURCE or int(nt) == int(NodeType.SOURCE):
            return True
    return False


def _closed_retrieval_result(sn, cache_idx, res, hitprob, missprob, delayedprob, it, method, t0):
    """Wrap a da_cacheqn_retrieval result into the retrieval-analyzer namespace."""
    K = sn.nclasses
    M = sn.nstations
    hp = np.full((1, K), np.nan)
    mp = np.full((1, K), np.nan)
    dp = np.full((1, K), np.nan)
    lat = np.full((1, K), np.nan)
    hitprob = np.asarray(hitprob, dtype=float).ravel()
    missprob = np.asarray(missprob, dtype=float).ravel()
    delayedprob = np.asarray(delayedprob, dtype=float).ravel()
    for r in range(K):
        if hitprob[r] != 0.0 or missprob[r] != 0.0:
            hp[0, r] = hitprob[r]
            mp[0, r] = missprob[r]
            dp[0, r] = delayedprob[r]
    ch = sn.nodeparam[cache_idx]
    h = len(np.atleast_1d(np.asarray(ch.itemcap)).ravel())
    hitproblist = np.full((K, h), np.nan)
    QN = getattr(res, 'QN', None); UN = getattr(res, 'UN', None)
    RN = getattr(res, 'RN', None); TN = getattr(res, 'TN', None)
    XN = getattr(res, 'XN', None); CN = getattr(res, 'CN', None)
    lG = getattr(res, 'lG', np.nan)
    return SimpleNamespace(QN=QN, UN=UN, RN=RN, TN=TN, XN=XN, CN=CN, lG=lG,
                           method=method, runtime=time.time() - t0, it=it, pij=None,
                           hitprob=hp, missprob=mp, delayedprob=dp,
                           hitproblist=hitproblist, itemprob=None,
                           expected_latency=lat, cache_idx=cache_idx)


def solver_mva_cacheqn_retrieval_analyzer(sn, options=None):
    """MVA analyzer for a CLOSED integrated cache-queueing model with a delayed-hit
    retrieval system. Delegates to da_cacheqn_retrieval with a load-dependent MVA
    network solve. Native-Python port of
    matlab/src/solvers/MVA/solver_mva_cacheqn_retrieval_analyzer.m."""
    t0 = time.time()
    from ..da import da_cacheqn_retrieval
    from ..solvers.mva.analyzers import solver_mva_analyzer

    cache_idx = -1
    for i, nt in enumerate(sn.nodetype):
        if nt == NodeType.CACHE or int(nt) == int(NodeType.CACHE):
            cache_idx = i
            break

    def netfun(snit):
        # Exact (load-independent) MVA aggregation. solver_mva_analyzer's load-
        # dependent amvald path is incompatible with the mutated class-switch struct;
        # the exact MVA still yields the correct cache hit/miss probabilities (the
        # fetch-station coalescing scaling only affects the throughput magnitude, which
        # is EXPERIMENTAL / understated - see da_cacheqn_retrieval LIMITATIONS).
        return solver_mva_analyzer(snit, _copy_options(options))

    res, hitprob, missprob, delayedprob, it, sn = da_cacheqn_retrieval(sn, netfun, options)
    return _closed_retrieval_result(sn, cache_idx, res, hitprob, missprob, delayedprob, it, "fpi", t0)


def solver_nc_cacheqn_retrieval_analyzer(sn, options=None):
    """NC analyzer for a CLOSED integrated cache-queueing model with a delayed-hit
    retrieval system. Delegates to da_cacheqn_retrieval with an NC (load-dependent
    ncld when the fetch station is present) network solve. Native-Python port of
    matlab/src/solvers/NC/solver_nc_cacheqn_retrieval_analyzer.m."""
    t0 = time.time()
    from ..da import da_cacheqn_retrieval
    from ..solvers.nc.analyzers import solver_nc_analyzer, solver_ncld_analyzer

    cache_idx = -1
    for i, nt in enumerate(sn.nodetype):
        if nt == NodeType.CACHE or int(nt) == int(NodeType.CACHE):
            cache_idx = i
            break

    def netfun(snit):
        lld = getattr(snit, 'lldscaling', None)
        cds = getattr(snit, 'cdscaling', None)
        jds = getattr(snit, 'jdscaling', None)
        opts = _copy_options(options)
        if (lld is not None and np.size(lld) > 0) or (cds is not None and np.size(cds) > 0) \
                or (jds is not None and np.size(jds) > 0):
            return solver_ncld_analyzer(snit, opts)
        return solver_nc_analyzer(snit, opts)

    res, hitprob, missprob, delayedprob, it, sn = da_cacheqn_retrieval(sn, netfun, options)
    return _closed_retrieval_result(sn, cache_idx, res, hitprob, missprob, delayedprob, it, "fpi", t0)


def _copy_options(options):
    import copy as _copy
    if options is None:
        return None
    try:
        return _copy.copy(options)
    except Exception:
        return options


def has_retrieval_cache(sn):
    """True if the model contains a Cache with a delayed-hit retrieval system."""
    if getattr(sn, 'nodeparam', None) is None:
        return False
    for idx, nt in enumerate(sn.nodetype):
        if (nt == NodeType.CACHE or int(nt) == int(NodeType.CACHE)) and idx in sn.nodeparam:
            if int(getattr(sn.nodeparam[idx], 'retrieval_system_capacity', 0)) > 0:
                return True
    return False
