"""
Heterogeneous server pools in the CTMC state space.

Port of MATLAB State.afterEventStationPool and State.fromMarginalPool. A pooled
station has been rewritten by api.solvers.ctmc.pools.ctmc_pools: the class-r
server block holds one sub-block of phases per compatible pool (ascending pool
order), so the pool of a job in service is read from its phase, and the pool
bookkeeping lives in sn.nodeparam[ind]['ctmcpool']. The local state is
[buffer | servers | local variables]: the buffer holds the WAITING jobs only
(1-based class ids right-aligned, the rightmost is the oldest; per-class counts
under SIRO). The semantics are those of LDES:

- ARV: the candidate pools are the compatible pools with a free server. None:
  the job waits. Otherwise the HeteroSchedPolicy picks the pool: ORDER the
  lowest index, ALFS the least flexible pool, FSF the highest rate for the class
  (ties by index), RAIS each candidate with equal probability, ALIS and FAIRNESS
  the first candidate of the rotating pool order, which then moves to the back
  when there was more than one candidate.
- DEP: the freed server of pool t takes the first waiting job, in the
  discipline's service order, that pool t can serve.
- PHASE: a transition inside the block of the job's pool.

Invariant: no waiting job has a free compatible server.
"""

import numpy as np

from ...constants import SchedStrategy, EventType, HeteroSchedPolicy


def pool_of(sn, ind):
    """The ctmcpool record of node ind, or None when the node is not pooled."""
    npar = getattr(sn, 'nodeparam', None)
    if npar is None:
        return None
    p = npar.get(ind) if isinstance(npar, dict) else (npar[ind] if ind < len(npar) else None)
    if isinstance(p, dict):
        return p.get('ctmcpool')
    return getattr(p, 'ctmcpool', None) if p is not None else None


def _busy(pool, srv, Ks):
    busy = np.zeros(pool['ntypes'])
    for r, plist in enumerate(pool['pools']):
        for k, t in enumerate(plist):
            c0 = Ks[r] + pool['off'][r][k]
            busy[t] += np.sum(srv[c0:c0 + pool['len'][r][k]])
    return busy


def _classcounts(buf, srv, K, Ks, is_siro, R):
    nir = np.zeros(R)
    for r in range(R):
        nir[r] = np.sum(srv[Ks[r]:Ks[r] + K[r]])
        nir[r] += buf[r] if is_siro else np.sum(buf == r + 1)
    return nir


def _physical_capacity(sn, ist, r):
    dr = getattr(sn, 'droprule', None)
    if dr is None:
        return False
    try:
        from ..sn.network_struct import DropStrategy as SnDropStrategy
        v = int(np.asarray(dr)[ist, r])
    except (IndexError, TypeError, ValueError):
        return False
    return v != int(SnDropStrategy.WAITQ) and v != 0


def _stack(rows_list):
    """Stack rows of possibly different buffer widths, left-padding the narrower ones."""
    if not rows_list:
        return np.zeros((0, 0))
    w = max(len(r) for r in rows_list)
    return np.array([np.concatenate([np.zeros(w - len(r)), r]) for r in rows_list], dtype=float)


def _arrival(sn, ist, pool, buf, srv, var, cls, Ks, K, is_siro):
    from .after_event_station import _arrival_is_lost
    R = sn.nclasses
    out = []
    if not pool['pools'][cls]:
        return out
    nir = _classcounts(buf, srv, K, Ks, is_siro, R)
    ni = np.sum(nir)
    cap = float(np.asarray(sn.cap).ravel()[ist]) if getattr(sn, 'cap', None) is not None else np.inf
    ccap = np.inf
    if getattr(sn, 'classcap', None) is not None:
        ccap = float(np.asarray(sn.classcap)[ist, cls])
    if ni >= cap or nir[cls] >= ccap:
        if _physical_capacity(sn, ist, cls) and _arrival_is_lost(sn, ist, cls):
            out.append((np.concatenate([buf, srv, var]), 1.0, -1))  # lost: state unchanged
        return out  # blocked, or beyond the state-space cutoff
    busy = _busy(pool, srv, Ks)
    cand = [t for t in pool['pools'][cls] if busy[t] < pool['count'][t]]
    if not cand:
        b = buf.copy()
        if is_siro:
            b[cls] += 1
        else:
            zeros = np.flatnonzero(b == 0)
            if zeros.size == 0:
                b = np.concatenate([[0.0], b])
                slot = 0
            else:
                slot = zeros[-1]
            b[slot] = cls + 1
        out.append((np.concatenate([b, srv, var]), 1.0, -1))
        return out
    choice = [cand[0]]
    pch = [1.0]
    varc = var.copy()
    if len(cand) > 1:
        pol = int(pool['policy'])
        if pol == int(HeteroSchedPolicy.ALFS):
            choice = [next(t for t in pool['alfsorder'] if t in cand)]
        elif pol == int(HeteroSchedPolicy.FSF):
            rates = [pool['fsfrate'][t, cls] for t in cand]
            choice = [cand[int(np.argmax(rates))]]  # first maximum wins ties
        elif pol == int(HeteroSchedPolicy.RAIS):
            choice = list(cand)
            pch = [1.0 / len(cand)] * len(cand)
        elif pol in (int(HeteroSchedPolicy.ALIS), int(HeteroSchedPolicy.FAIRNESS)) and pool['rotate']:
            vp = pool['varpos']
            order = list(pool['perms'][int(var[vp]) - 1])
            for k, t in enumerate(order):
                if t in cand:
                    choice = [t]
                    order = order[:k] + order[k + 1:] + [t]
                    break
            varc[vp] = pool['permindex'][tuple(order)] + 1
    for t, pc in zip(choice, pch):
        k = pool['pools'][cls].index(t)
        a = pool['alpha'][cls][k]
        for j, aj in enumerate(a):
            if aj <= 0:
                continue
            s2 = srv.copy()
            s2[Ks[cls] + pool['off'][cls][k] + j] += 1
            out.append((np.concatenate([buf, s2, varc]), pc * aj, cls))
    return out


def _start(pool, buf, srv, s, t, Ks):
    k = pool['pools'][s].index(t)
    res = []
    for j, aj in enumerate(pool['alpha'][s][k]):
        if aj <= 0:
            continue
        s2 = srv.copy()
        s2[Ks[s] + pool['off'][s][k] + j] += 1
        res.append((buf, s2, aj, s))
    return res


def _serve_next(sn, ist, pool, buf, srv, t, Ks, is_siro, sched):
    """The freed server of pool t takes the next compatible waiting job, if any."""
    compat = pool['compat']
    if is_siro:
        elig = [s for s in range(len(buf)) if buf[s] > 0 and compat[t, s]]
        if not elig:
            return [(buf, srv, 1.0, -1)]
        tot = sum(buf[s] for s in elig)
        res = []
        for s in elig:
            b = buf.copy()
            b[s] -= 1
            for bb, ss, p, c in _start(pool, b, srv, s, t, Ks):
                res.append((bb, ss, p * buf[s] / tot, c))
        return res
    pos = [p for p in range(len(buf)) if buf[p] > 0 and compat[t, int(buf[p]) - 1]]
    if not pos:
        return [(buf, srv, 1.0, -1)]
    if sched in (SchedStrategy.HOL, SchedStrategy.FCFSPRIO, SchedStrategy.LCFSPRIO):
        prio = np.array([sn.classprio[int(buf[p]) - 1] for p in pos])
        pos = [p for p, q in zip(pos, prio) if q == np.min(prio)]  # a lower value is a higher priority
    if sched in (SchedStrategy.LCFS, SchedStrategy.LCFSPRIO):
        p = pos[0]  # leftmost is the newest
    else:
        p = pos[-1]  # rightmost is the oldest
    s = int(buf[p]) - 1
    b = np.concatenate([[0.0], buf[:p], buf[p + 1:]])
    return _start(pool, b, srv, s, t, Ks)


def after_event_station_pool(sn, ind, ist, inspace, event, cls, V, is_simulation=False):
    """Successors of EVENT for class CLS at pooled station IND.

    Returns (outspace, outrate, outprob, outstart, outpreempt), as api.state.after_event.
    """
    R = sn.nclasses
    pool = pool_of(sn, ind)
    K = np.array(sn.phasessz[ist], dtype=int)
    Ks = np.array(sn.phaseshift[ist], dtype=int)
    sched = sn.sched[ist]
    is_siro = sched == SchedStrategy.SIRO
    nsrv = int(np.sum(K))
    rows, rates, probs, starts = [], [], [], []
    inspace = np.atleast_2d(np.asarray(inspace, dtype=float))
    for st in inspace:
        W = st.size - nsrv - V
        buf = st[:W].copy()
        srv = st[W:W + nsrv].copy()
        var = st[W + nsrv:].copy()
        if event == EventType.ARV:
            for row, p, c in _arrival(sn, ist, pool, buf, srv, var, cls, Ks, K, is_siro):
                rows.append(row); rates.append(-1.0); probs.append(p); starts.append(c)
        elif event == EventType.DEP:
            for k, t in enumerate(pool['pools'][cls]):
                ex = pool['exit'][cls][k]
                for j in range(pool['len'][cls][k]):
                    col = Ks[cls] + pool['off'][cls][k] + j
                    nj = srv[col]
                    if nj <= 0 or ex[j] <= 0:
                        continue
                    sd = srv.copy()
                    sd[col] -= 1
                    for bb, ss, p, c in _serve_next(sn, ist, pool, buf, sd, t, Ks, is_siro, sched):
                        rows.append(np.concatenate([bb, ss, var])); rates.append(ex[j] * nj * p)
                        probs.append(1.0); starts.append(c)
        elif event == EventType.PHASE:
            for k, t in enumerate(pool['pools'][cls]):
                D0 = pool['D0'][cls][k]
                c0 = Ks[cls] + pool['off'][cls][k]
                for j in range(pool['len'][cls][k]):
                    nj = srv[c0 + j]
                    if nj <= 0:
                        continue
                    for jd in range(pool['len'][cls][k]):
                        if jd == j or D0[j, jd] <= 0:
                            continue
                        sp = srv.copy()
                        sp[c0 + j] -= 1
                        sp[c0 + jd] += 1
                        rows.append(np.concatenate([buf, sp, var])); rates.append(D0[j, jd] * nj)
                        probs.append(1.0); starts.append(-1)
    if not rows:
        return (np.zeros((0, 0)), np.zeros((0, 0)), np.ones((1, 1)),
                np.zeros((0, R)), np.zeros((0, R)))
    outspace = _stack(rows)
    outrate = np.array(rates, dtype=float).reshape(-1, 1)
    outprob = np.array(probs, dtype=float).reshape(-1, 1)
    outstart = np.zeros((len(rows), R))
    for i, c in enumerate(starts):
        if c >= 0:
            outstart[i, c] = 1.0
    outpreempt = np.zeros((len(rows), R))
    if is_simulation and len(rows) > 1:
        w = outprob.ravel() if event == EventType.ARV else outrate.ravel()
        cum = np.cumsum(w) / np.sum(w)
        fc = int(np.searchsorted(cum, np.random.rand(), side='right'))
        fc = min(fc, len(rows) - 1)
        tot = float(np.sum(outrate)) if event != EventType.ARV else -1.0
        outspace = outspace[fc:fc + 1]
        outrate = np.array([[tot]])
        outprob = np.ones((1, 1))
        outstart = outstart[fc:fc + 1]
        outpreempt = outpreempt[fc:fc + 1]
    return outspace, outrate, outprob, outstart, outpreempt


def from_marginal_pool(sn, ind, n):
    """Local [buffer | servers] states with marginal n at a pooled station.

    Every split of the jobs into in-service (per class and pool, within the pool
    sizes) and waiting ones such that no waiting job has a free compatible server,
    every phase assignment of the jobs in service, and every order of the waiting
    jobs (per-class counts under SIRO). The ordered buffer is
    max(1, min(sum(n), cap)) wide. The pool-order column is appended by ctmc_ssg,
    like every other Python local variable.
    """
    from .multiset_perms import multiset_perms
    R = sn.nclasses
    ist = int(sn.nodeToStation[ind])
    pool = pool_of(sn, ind)
    K = np.array(sn.phasessz[ist], dtype=int)
    n = np.asarray(n, dtype=int).ravel()
    is_siro = sn.sched[ist] == SchedStrategy.SIRO
    cap = float(np.asarray(sn.cap).ravel()[ist]) if getattr(sn, 'cap', None) is not None else np.inf
    W = R if is_siro else int(max(1, min(np.sum(n), cap)))
    nsrv = int(np.sum(K))
    unserved = [r for r in range(R) if not pool['pools'][r]]
    if np.sum(n) > cap or any(n[r] > 0 for r in unserved):
        return np.zeros((0, W + nsrv))
    allocs = []

    def rec(r, k, x, load):
        if r >= R:
            for s in range(R):
                if not pool['pools'][s]:
                    continue
                if n[s] - sum(x[s]) > 0 and any(load[t] < pool['count'][t] for t in pool['pools'][s]):
                    return
            allocs.append([list(v) for v in x])
            return
        if not pool['pools'][r] or k >= len(pool['pools'][r]):
            rec(r + 1, 0, x, load)
            return
        t = pool['pools'][r][k]
        used = sum(x[r][:k])
        for c in range(0, int(min(n[r] - used, pool['count'][t] - load[t])) + 1):
            x[r].append(c)
            load[t] += c
            rec(r, k + 1, x, load)
            load[t] -= c
            x[r].pop()

    rec(0, 0, [[] for _ in range(R)], np.zeros(pool['ntypes']))
    rows = []
    for x in allocs:
        w = n.copy()
        for r in range(R):
            w[r] -= sum(x[r])
        srvset = [np.zeros(0)]
        for r in range(R):
            if not pool['pools'][r]:
                blocks = [np.zeros(K[r])]
            else:
                blocks = [np.zeros(0)]
                for k in range(len(pool['pools'][r])):
                    comps = _compositions(int(x[r][k]), int(pool['len'][r][k]))
                    blocks = [np.concatenate([b, c]) for b in blocks for c in comps]
                blocks = [np.concatenate([b, np.zeros(K[r] - b.size)]) for b in blocks]
            srvset = [np.concatenate([s, b]) for s in srvset for b in blocks]
        if is_siro:
            bufset = [w.astype(float)]
        elif np.sum(w) == 0:
            bufset = [np.zeros(W)]
        else:
            vi = []
            for r in range(R):
                vi.extend([r + 1] * int(w[r]))
            perms = {tuple(p) for p in multiset_perms(vi)}
            bufset = [np.concatenate([np.zeros(W - len(p)), np.array(p, dtype=float)]) for p in sorted(perms)]
        for b in bufset:
            for s in srvset:
                rows.append(np.concatenate([b, s]))
    if not rows:
        return np.zeros((0, W + nsrv))
    space = np.unique(np.array(rows, dtype=float), axis=0)
    return space[::-1]


def _compositions(total, parts):
    """All vectors of PARTS non-negative integers summing to TOTAL."""
    if parts == 0:
        return [np.zeros(0)] if total == 0 else []
    if parts == 1:
        return [np.array([float(total)])]
    out = []
    for first in range(total, -1, -1):
        for rest in _compositions(total - first, parts - 1):
            out.append(np.concatenate([[float(first)], rest]))
    return out
