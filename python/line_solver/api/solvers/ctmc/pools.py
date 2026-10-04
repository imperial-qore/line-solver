"""
Heterogeneous server pools (Queue.add_server_type) for the CTMC state space.

Port of MATLAB solver_ctmc_pools.m, following the LDES semantics. A class-r job
in service at a pooled station occupies one server of ONE compatible pool t and
is served by that pool's law (set_hetero_service), or by the station's own law
for r when the pool declares none. The per-class service process is replaced by
a block-diagonal phase-type law with one block per compatible pool, in
ascending pool order, so a phase of the class-r server block identifies both the
pool and the service phase. The bookkeeping read by
api.state.after_event_station_pool is stored in sn.nodeparam[ind]['ctmcpool']:

  ntypes, count[t], compat[t, r]      pools, their sizes, class compatibility
  policy                              HeteroSchedPolicy of the station
  pools[r]                            compatible pools of r, ascending (0-based)
  off[r][k], len[r][k]                offset and length of block k of class r
  alpha[r][k], exit[r][k], D0[r][k]   entry vector, exit rates, D0 of block k
  fsfrate[t, r]                       1/mean of the law of r at pool t (FSF)
  alfsorder                           pools by ascending number of classes (ALFS)
  rotate, perms, permindex, varpos    ALIS/FAIRNESS pool order, see below

ALIS and FAIRNESS keep a global pool order: a job picks the first pool of that
order among those with a free compatible server, and that pool moves to the back
only when it had more than one candidate. When some class has two or more
compatible pools the order is part of the state, as a 1-based index into perms
held in the station's local-variable block (position varpos), which ctmc_ssg
enumerates. The rewrite is idempotent, and it rebuilds sn.state for the pooled
stations in the new layout.
"""

import copy
import itertools

import numpy as np

from ....constants import ProcessType, EventType, HeteroSchedPolicy, GlobalConstants
from ....lang.base import SchedStrategy

_POOL_SCHED = (SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRIO,
               SchedStrategy.LCFS, SchedStrategy.LCFSPRIO, SchedStrategy.SIRO)


def _param(sn, ind):
    npar = getattr(sn, 'nodeparam', None)
    if npar is None:
        return None
    return npar.get(ind) if isinstance(npar, dict) else (npar[ind] if ind < len(npar) else None)


def _pget(p, key, default=None):
    if isinstance(p, dict):
        return p.get(key, default)
    return getattr(p, key, default)


def _pset(p, key, value):
    if isinstance(p, dict):
        p[key] = value
    else:
        setattr(p, key, value)


def _name(sn, ind):
    names = getattr(sn, 'nodenames', None)
    return names[ind] if names is not None and ind < len(names) else str(ind)


def _cname(sn, r):
    names = getattr(sn, 'classnames', None)
    return names[r] if names is not None and r < len(names) else str(r)


def _law(proc_ir):
    """(D0, D1) of a stored law, or None when the class is not served (disabled/absent)."""
    if proc_ir is None:
        return None
    from ...sn.proc_form import proc_to_map
    D0, D1 = proc_to_map(proc_ir)
    if D0 is None:
        return None
    D0 = np.atleast_2d(np.asarray(D0, dtype=float))
    D1 = np.atleast_2d(np.asarray(D1, dtype=float))
    if np.any(np.isnan(D0)):
        return None
    return D0, D1


def _require_ph(law, procid, sn, ind, r, what):
    bad = (ProcessType.MAP, ProcessType.MMPP2, ProcessType.MMAP, ProcessType.ME, ProcessType.RAP)
    if procid is not None and any(procid == b or getattr(procid, 'value', None) == b.value for b in bad):
        raise RuntimeError(
            "SolverCTMC serves a heterogeneous server pool with a renewal phase-type law; class '%s' at station '%s' "
            "has a %s law in %s." % (_cname(sn, r), _name(sn, ind), getattr(procid, 'name', procid), what))
    D0, D1 = law
    ex = -np.sum(D0, axis=1)
    a = np.sum(D1, axis=0)
    a = a / max(np.sum(a), GlobalConstants.FineTol)
    offd = D0 - np.diag(np.diag(D0))
    if (np.any(offd < -GlobalConstants.FineTol) or np.any(ex < -GlobalConstants.FineTol)
            or np.sum(np.abs(D1 - np.outer(ex, a))) > GlobalConstants.CoarseTol * max(1.0, np.sum(np.abs(D1)))):
        raise RuntimeError(
            "SolverCTMC serves a heterogeneous server pool with a renewal phase-type law; class '%s' at station '%s' "
            "has a correlated or matrix-exponential law in %s." % (_cname(sn, r), _name(sn, ind), what))


def _refuse_features(sn, ind, ist, R, p):
    why = ''
    lld = getattr(sn, 'lldscaling', None)
    cd = getattr(sn, 'cdscaling', None)
    jd = getattr(sn, 'jdscaling', None)

    def _entry(c):
        if c is None:
            return None
        if isinstance(c, dict):
            return c.get(ist)
        return c[ist] if ist < len(c) else None

    lld = np.asarray(lld) if lld is not None else None
    if lld is not None and lld.ndim == 2 and lld.shape[0] > ist and np.any(lld[ist, :] != 1):
        why = 'load-dependent service'
    elif _entry(cd) is not None:
        why = 'class-dependent service'
    elif _entry(jd) is not None:
        why = 'joint-dependent service'
    elif getattr(sn, 'gdscaling', None) is not None:
        why = 'global dependence'
    elif getattr(sn, 'hasbreakdown', None) is not None and np.asarray(sn.hasbreakdown).ravel().size > ind \
            and int(np.asarray(sn.hasbreakdown).ravel()[ind]) == 1:
        why = 'server breakdowns'
    elif getattr(sn, 'retrialProc', None) is not None and any(x is not None for x in sn.retrialProc[ist]):
        why = 'retrial'
    elif getattr(sn, 'balkingStrategy', None) is not None and np.any(np.asarray(sn.balkingStrategy)[ist, :] > 0):
        why = 'balking'
    elif getattr(sn, 'impatienceClass', None) is not None and np.any(np.asarray(sn.impatienceClass)[ist, :] > 0):
        why = 'reneging'
    elif getattr(sn, 'isbasblocking', None) is not None and np.asarray(sn.isbasblocking).ravel().size > ind \
            and int(np.asarray(sn.isbasblocking).ravel()[ind]) == 1:
        why = 'BAS blocking'
    elif getattr(sn, 'immfeed', None) is not None and np.asarray(sn.immfeed).ndim == 2 \
            and np.any(np.asarray(sn.immfeed)[ist, :]):
        why = 'immediate feedback'
    elif getattr(sn, 'replyblock', None) is not None and np.asarray(sn.replyblock).ndim == 2 \
            and np.asarray(sn.replyblock).shape[0] > ind and np.any(np.asarray(sn.replyblock)[ind, :] > 0):
        why = 'synchronous calls'
    elif _pget(p, 'serverparallelism') is not None and np.any(np.asarray(_pget(p, 'serverparallelism')) > 1):
        why = 'server parallelism'
    else:
        rt = getattr(sn, 'routing', None)
        if rt is not None:
            rtv = rt.get(ind) if isinstance(rt, dict) else (rt[ind] if ind < len(rt) else None)
            vals = rtv.values() if isinstance(rtv, dict) else (rtv if isinstance(rtv, (list, tuple, np.ndarray)) else [rtv])
            for v in vals:
                nm = getattr(v, 'name', str(v))
                if nm in ('RROBIN', 'WRROBIN'):
                    why = 'round-robin routing'
                    break
        issig = getattr(sn, 'issignal', None)
        rtn = getattr(sn, 'rtnodes', None)
        if not why and issig is not None and np.any(issig) and rtn is not None:
            rtn = np.asarray(rtn)
            for s in np.flatnonzero(np.asarray(issig).ravel()):
                if np.any(rtn[:, ind * R + s] > 0):
                    why = 'signals'
                    break
    if why:
        raise RuntimeError("SolverCTMC does not combine heterogeneous server pools with %s (station '%s')."
                           % (why, _name(sn, ind)))


def _pool_station(sn, ind, ist, R, p):
    T = int(_pget(p, 'nservertypes'))
    compat = np.asarray(_pget(p, 'servercompat'), dtype=float) > 0
    count = np.asarray(_pget(p, 'serverspertype'), dtype=float).ravel().astype(int)
    policy = _pget(p, 'heteroschedpolicy')
    policy = HeteroSchedPolicy.ORDER if policy is None else HeteroSchedPolicy(int(policy))
    heteroproc = _pget(p, 'heteroproc') or {}
    tnames = _pget(p, 'servertypenames') or [str(t) for t in range(T)]
    if sn.sched[ist] not in _POOL_SCHED:
        raise RuntimeError(
            "SolverCTMC supports heterogeneous server pools under FCFS, HOL, FCFSPRIO, LCFS, LCFSPRIO and SIRO; "
            "station '%s' uses %s." % (_name(sn, ind), getattr(sn.sched[ist], 'name', sn.sched[ist])))
    _refuse_features(sn, ind, ist, R, p)
    procid = np.asarray(sn.procid) if getattr(sn, 'procid', None) is not None else None
    pools = [[] for _ in range(R)]
    off = [[] for _ in range(R)]
    ln = [[] for _ in range(R)]
    alpha = [[] for _ in range(R)]
    exitr = [[] for _ in range(R)]
    D0b = [[] for _ in range(R)]
    fsfrate = np.zeros((T, R))
    for r in range(R):
        base = _law(sn.proc[ist][r]) if r < len(sn.proc[ist]) else None
        has_pool_law = any(compat[t, r] and heteroproc.get((t, r)) is not None for t in range(T))
        if base is None and not has_pool_law:
            continue  # the class is not served here
        pools[r] = [int(t) for t in np.flatnonzero(compat[:, r])]
        if not pools[r]:
            raise RuntimeError("Station '%s' declares no server pool compatible with class '%s', so a job of that "
                               "class would wait forever." % (_name(sn, ind), _cname(sn, r)))
        if base is not None:
            _require_ph(base, procid[ist, r] if procid is not None else None, sn, ind, r, 'its default service')
        blocks = []
        shift = 0
        for t in pools[r]:
            law = heteroproc.get((t, r))
            if law is None:
                if base is None:
                    raise RuntimeError(
                        "Server pool '%s' of station '%s' accepts class '%s' but declares no law for it, and the "
                        "station's own service for the class is disabled." % (tnames[t], _name(sn, ind), _cname(sn, r)))
                law = base
            else:
                law = (np.atleast_2d(np.asarray(law[0], dtype=float)), np.atleast_2d(np.asarray(law[1], dtype=float)))
                _require_ph(law, None, sn, ind, r, "pool '%s'" % tnames[t])
            D0, D1 = law
            ex = -np.sum(D0, axis=1)
            a = np.sum(D1, axis=0)
            a = a / np.sum(a)
            off[r].append(shift)
            ln[r].append(D0.shape[0])
            alpha[r].append(a)
            exitr[r].append(ex)
            D0b[r].append(D0)
            fsfrate[t, r] = 1.0 / float(a @ np.linalg.solve(-D0, np.ones(D0.shape[0])))
            blocks.append(D0)
            shift += D0.shape[0]
        n = shift
        D0x = np.zeros((n, n))
        for k, B in enumerate(blocks):
            o = off[r][k]
            D0x[o:o + B.shape[0], o:o + B.shape[0]] = B
        exx = -np.sum(D0x, axis=1)
        piex = np.zeros(n)
        piex[:ln[r][0]] = alpha[r][0]
        sn.proc[ist][r] = [D0x, np.outer(exx, piex)]
        sn.pie[ist][r] = piex
        sn.mu[ist][r] = -np.diag(D0x)
        sn.phi[ist][r] = exx / (-np.diag(D0x))
        sn.phases[ist, r] = n
        sn.phasessz[ist, r] = n
        if procid is not None:
            sn.procid[ist, r] = ProcessType.PH
        if getattr(sn, 'isph', None) is not None and np.asarray(sn.isph).ndim == 2:
            sn.isph[ist, r] = True
    ncls = np.sum(compat, axis=1)
    alfsorder = [int(t) for t in np.argsort(ncls, kind='stable')]
    rotate = (policy in (HeteroSchedPolicy.ALIS, HeteroSchedPolicy.FAIRNESS)) and any(len(pl) >= 2 for pl in pools)
    perms = sorted(itertools.permutations(range(T))) if rotate else []  # identity first
    pool = {
        'ntypes': T, 'count': count, 'compat': compat, 'policy': int(policy), 'pools': pools,
        'off': off, 'len': ln, 'alpha': alpha, 'exit': exitr, 'D0': D0b, 'fsfrate': fsfrate,
        'alfsorder': alfsorder, 'rotate': rotate, 'perms': perms,
        'permindex': {pm: i for i, pm in enumerate(perms)}, 'varpos': 0,
    }
    _pset(p, 'ctmcpool', pool)
    sn.nservers[ist] = int(np.sum(count))  # the pools are the server bank, as in LDES
    shift = 0
    for r in range(R):
        sn.phaseshift[ist, r] = shift
        shift += int(sn.phasessz[ist, r])


def _rebuild_state(sn_old, sn, ind, ist, oldrow):
    """Replay the declared initial jobs as arrivals into an empty pooled station.

    In-service jobs first, by class, then the waiting jobs from the oldest to the
    newest; each takes the most likely outcome (RAIS: the first candidate pool).
    """
    from ...state.after_event_station_pool import after_event_station_pool
    R = sn.nclasses
    Kold = np.array(sn_old.phasessz[ist], dtype=int)
    Ksold = np.array(sn_old.phaseshift[ist], dtype=int)
    V = int(np.sum(sn_old.nvars[ind])) if getattr(sn_old, 'nvars', None) is not None else 0
    oldrow = np.asarray(oldrow, dtype=float).ravel()
    nsrv = int(np.sum(Kold))
    W = oldrow.size - nsrv - V
    buf = oldrow[:W]
    srv = oldrow[W:W + nsrv]
    seq = []
    for r in range(R):
        seq += [r] * int(np.sum(srv[Ksold[r]:Ksold[r] + Kold[r]]))
    is_siro = sn.sched[ist] == SchedStrategy.SIRO
    if is_siro:
        for r in range(R):
            seq += [r] * int(buf[r])
    else:
        seq += [int(b) - 1 for b in buf[buf > 0][::-1]]  # rightmost is the oldest
    pool = _pget(_param(sn, ind), 'ctmcpool')
    Wb = R if is_siro else max(1, len(seq))
    var = [1.0] if pool['rotate'] else []
    Vn = len(var)
    st = np.concatenate([np.zeros(Wb), np.zeros(int(np.sum(sn.phasessz[ist]))), var])
    saved = sn.cap
    sn.cap = np.full(np.asarray(saved).shape, np.inf)  # the declared state is placed without the cutoff
    try:
        for c in seq:
            outs, _, outp, _, _ = after_event_station_pool(sn, ind, ist, st, EventType.ARV, c, Vn)
            if outs.shape[0] == 0:
                raise RuntimeError("The initial state of station '%s' cannot be placed on its server pools."
                                   % _name(sn, ind))
            st = outs[int(np.argmax(outp.ravel()))]
    finally:
        sn.cap = saved
    return st


def ctmc_pools(sn):
    """Rewrite every heterogeneous-server station into pooled form (in place)."""
    if getattr(sn, 'nodeparam', None) is None:
        return sn
    R = int(sn.nclasses)
    changed = []
    sn_old = None
    for ind in range(int(sn.nnodes)):
        if not sn.isstation[ind]:
            continue
        p = _param(sn, ind)
        if p is None or not _pget(p, 'nservertypes') or _pget(p, 'ctmcpool') is not None:
            continue
        ist = int(sn.nodeToStation[ind])
        # PAS/OI stations model compatible servers through the OI rank rate, not through pools
        if sn.sched[ist] in (SchedStrategy.PAS, getattr(SchedStrategy, 'OI', None)):
            continue
        if sn_old is None:
            sn_old = copy.copy(sn)
            sn_old.phasessz = np.array(sn.phasessz, copy=True)
            sn_old.phaseshift = np.array(sn.phaseshift, copy=True)
            sn_old.nvars = None if sn.nvars is None else np.array(sn.nvars, copy=True)
        _pool_station(sn, ind, ist, R, p)
        changed.append(ind)
    for ind in changed:
        ist = int(sn.nodeToStation[ind])
        isf = int(sn.nodeToStateful[ind])
        state = getattr(sn, 'state', None)
        if state is None:
            continue
        old = state.get(isf) if hasattr(state, 'get') else (state[isf] if isf < len(state) else None)
        if old is None or np.asarray(old).size == 0:
            continue
        state[isf] = _rebuild_state(sn_old, sn, ind, ist, np.atleast_2d(old)[0]).reshape(1, -1)
    return sn
