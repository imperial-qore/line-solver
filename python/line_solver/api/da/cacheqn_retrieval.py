"""Decomposition-aggregation driver for a CLOSED integrated cache-queueing model
whose Cache node has a delayed-hit retrieval system (Cache.set_retrieval_system).

Native-Python port of matlab/src/api/da/da_cacheqn_retrieval.m and
jar/.../jline/api/da/Da_cacheqn_retrieval.java.

The cache is relabelled as a ClassSwitch (read -> hit / miss) and the finite-
population delayed-hit coalescing is made to emerge from the closed AMVA by giving
the retrieval (fetch) station a load-dependent COALESCING service rate
lldscaling(k)=k/d(k), d(k)=n_eff*(1-(1-1/n_eff)^k), n_eff=nitems-totalcap. The
distinct-fetch throughput saturates at 1/F and the coalescing benefit emerges from
the finite population.

LIMITATIONS (EXPERIMENTAL - see the MATLAB reference): hitprob/missprob are the
true cache probabilities (hit=P(cached), miss=1-hit) and accurate; the delayed
fraction is NOT recovered (delayedprob=0, folds into miss); the coalescing
throughput benefit is captured only in direction and understated; single fetch
station (single-backend) only.

NO CLOSED-RETRIEVAL EXAMPLE ships in the suite (every retrieval_* example is
OPEN). On a hand-built closed model this driver readily hits a reducible /
singular routing and a zero read-rate denominator at ``lam[r] = Xr*nv[ci,r]/
denom`` below (NaN read rate), so a plain closed model can return an empty or
unreliable avg_table with LinAlg/RuntimeWarnings. Treat the closed path as
experimental and validate against LDES. The C++ SolverMVA (line-mp) deliberately
REFUSES this closed path by name; only the OPEN retrieval analyzer is ported
there. Verified 2026-07-24 that the open path matches across codebases while the
closed path lacks a clean, exampled reference.
"""
import numpy as np

from ..sn import NodeType
from ..cache.rrm import cache_gamma_lp
from ..cache.miss import cache_miss_fpi
from ..mc.dtmc import dtmc_stochcomp
from ..sn.transforms import sn_refresh_visits
from .fpi import da_fpi

__all__ = ['da_cacheqn_retrieval']


def _is_cache(nt):
    return nt == NodeType.CACHE or int(nt) == int(NodeType.CACHE)


def _is_classswitch_id():
    return int(NodeType.CLASSSWITCH)


def da_cacheqn_retrieval(sn, netfun, options):
    """Run the closed cache-retrieval decomposition-aggregation fixed point.

    Args:
        sn: NetworkStruct of the closed cache+retrieval model (mutated in place).
        netfun: callable (sn) -> result with .XN; the aggregation network solve
            (mvald / ncld when the load-dependent fetch station is present).
        options: solver options (uses iter_max, iter_tol).

    Returns:
        (res, hitprob, missprob, delayedprob, it, sn) where hitprob/missprob/
        delayedprob are (nclasses,) vectors set on the read class.
    """
    I = sn.nnodes
    K = sn.nclasses

    statefulNodesClasses = []
    for ind in range(I):
        if sn.isstateful[ind] != 0:
            statefulNodesClasses.extend([ind * K + r for r in range(K)])
    statefulNodesClasses = np.asarray(statefulNodesClasses, dtype=int)

    caches = [i for i in range(I) if _is_cache(sn.nodetype[i])]
    if len(caches) != 1:
        raise RuntimeError("da_cacheqn_retrieval requires exactly one Cache node.")
    ci = caches[0]
    ch = sn.nodeparam[ci]

    # retrieval configuration (single read class, single fetch station), 0-indexed
    readClass = next(iter(ch.retrieval_system_queue_indices.keys()))
    queueNodes = list(ch.retrieval_system_queue_indices[readClass])
    if len(queueNodes) != 1:
        raise RuntimeError("da_cacheqn_retrieval currently supports a single-station "
                           "(single-backend) retrieval system.")
    fetchNode = int(queueNodes[0])
    fetchStation = int(sn.nodeToStation[fetchNode])
    nitems = int(ch.nitems)
    totcap = float(np.sum(np.asarray(ch.itemcap, dtype=float)))
    n_eff = max(1.0, nitems - totcap)

    # closed population
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    Npop = int(round(np.sum(njobs[np.isfinite(njobs)])))
    Nfin = max(1, Npop)

    # mean fetch service F of the read class at the fetch station
    rcls0 = int(ch.retrieval_classes[0, readClass])
    F = 1.0 / float(sn.rates[fetchStation, rcls0])

    # load-dependent COALESCING rate on the fetch station: lldscaling(k)=k/d(k)
    alpha = np.ones(Nfin)
    for k in range(1, Nfin + 1):
        d = n_eff * (1.0 - (1.0 - 1.0 / n_eff) ** k)
        alpha[k - 1] = k / d
    lld = sn.lldscaling
    if lld is None or np.size(lld) == 0:
        lld = np.ones((sn.nstations, Nfin))
    elif lld.shape[1] < Nfin:
        grown = np.ones((sn.nstations, Nfin))
        grown[:, :lld.shape[1]] = lld
        lld = grown
    lld[fetchStation, :Nfin] = alpha
    sn.lldscaling = lld

    # relabel the cache as a class switch
    sn.nodetype[ci] = _is_classswitch_id()

    hitClass = np.atleast_1d(np.asarray(ch.hitclass)).astype(int).ravel()
    missClass = np.atleast_1d(np.asarray(ch.missclass)).astype(int).ravel()
    retrievalClasses = np.asarray(ch.retrieval_classes)
    pread = np.asarray(ch.pread[readClass], dtype=float).ravel()
    pread = pread / np.sum(pread)

    h = len(np.atleast_1d(np.asarray(ch.itemcap)).ravel())
    mMat = np.atleast_1d(np.asarray(ch.itemcap, dtype=float)).ravel()

    hitprob = np.zeros(K)
    missprob = np.zeros(K)
    delayedprob = np.zeros(K)
    res_holder = {}

    # default per-item routing R (miss chain shift) when accost is absent
    accost = getattr(ch, 'accost', None)

    def _build_R():
        if accost is not None:
            return accost
        R_cost = []
        for _ in range(nitems):
            Rmat = np.zeros((h + 1, h + 1))
            for l in range(h):
                Rmat[l, l + 1] = 1.0
            Rmat[h, h] = 1.0
            R_cost.append(Rmat)
        return [R_cost]

    R_cost = _build_R()

    def da_sweep(x, itnum):
        lam = np.array(x, dtype=float)
        r = readClass

        # isolated cache occupancy -> per-item uncached (would-be-miss) prob pi0
        lambd3d = np.zeros((1, nitems, h + 1))
        for i in range(nitems):
            lambd3d[0, i, :] = lam[r] * pread[i]
        gamma = np.asarray(cache_gamma_lp(lambd3d, R_cost)[0], dtype=float).reshape(nitems, h)
        # cache_miss_fpi without arrival rates returns the per-item miss probability
        # 1/(1+gamma*xi) in MI (index 2); pi0 (index 3) is only populated when lambd
        # is passed. Both use the same Che-approximation formula.
        pi0 = np.asarray(cache_miss_fpi(gamma, mMat)[2], dtype=float).ravel()

        nonhit = float(np.sum(pread * pi0))
        hp = 1.0 - nonhit

        # rebuild the cache class-switch routing (read -> hit self-loop; read ->
        # per-item retrieval class at the fetch station; retrieval class at the
        # fetch station -> miss at the cache)
        sn.rtnodes[ci * K + r, :] = 0.0
        sn.rtnodes[ci * K + r, ci * K + hitClass[r]] = hp
        for i in range(nitems):
            rcls = int(retrievalClasses[i, r])
            if rcls >= 0:
                sn.rtnodes[ci * K + r, fetchNode * K + rcls] = pread[i] * pi0[i]
                sn.rtnodes[fetchNode * K + rcls, :] = 0.0
                sn.rtnodes[fetchNode * K + rcls, ci * K + missClass[r]] = 1.0
                sn.rtnodes[ci * K + rcls, :] = 0.0

        sn.rt = dtmc_stochcomp(sn.rtnodes, statefulNodesClasses)
        sn.rt_visits = sn.rt   # force sn_refresh_visits to use the rebuilt routing
        sn_refresh_visits(sn)

        res = netfun(sn)
        res_holder['res'] = res

        # aggregate nodevisits over chains
        nv = None
        for val in sn.nodevisits.values():
            nv = np.asarray(val, dtype=float) if nv is None else nv + np.asarray(val, dtype=float)

        # throughput -> new read arrival rate at the cache
        chains = np.asarray(sn.chains)
        c = int(np.where(chains[:, r] != 0)[0][0])
        inchain = np.where(chains[c, :] != 0)[0]
        refnode = int(sn.stationToNode[int(sn.refstat.ravel()[r])])
        refclass_c = int(np.asarray(sn.refclass).ravel()[c])
        denom = nv[refnode, refclass_c] if refclass_c > -1 else nv[refnode, r]
        Xr = float(np.sum(np.asarray(res.XN, dtype=float).ravel()[inchain]))
        lam[r] = Xr * nv[ci, r] / denom

        hitprob[r] = hp
        missprob[r] = nonhit
        delayedprob[r] = 0.0

        return lam, np.array(x, dtype=float)

    lambda0 = np.zeros(K)
    lambda0[readClass] = 1.0
    iter_max = getattr(options, 'iter_max', 1000) or 1000
    iter_tol = getattr(options, 'iter_tol', 1e-6) or 1e-6
    _, it, _ = da_fpi(da_sweep, lambda0, iter_max, iter_tol,
                      norm=lambda a, b: float(np.sum(np.abs(np.asarray(a) - np.asarray(b)))))

    return res_holder.get('res'), hitprob, missprob, delayedprob, it, sn
