"""
Pair form of a multiserver MAP service for the CTMC state space.

Port of MATLAB solver_ctmc_mapc.m. JMT and both LDES engines sample a MAP (or
MMPP2) service at a station through ONE sampler per (station, class): every draw
starts in the phase the previous draw ENDED in, draws chained in service-start
order, idle periods included. With c > 1 servers the next start can happen while
earlier draws are still in progress, so the landing phase of a draw must be
known when it starts. With V = (-D0)^-1 D1, a draw started in h ends in j with
probability V(h,j). A busy server is therefore a PAIR (i,j): current phase i and
predetermined landing j, which moves i->k at D0(i,k)V(k,j)/V(i,j) and completes
at D1(i,j)/V(i,j). The per-class memory variable keeps its meaning of carried
phase h in 0..p-1, but it is the landing of the most recently STARTED draw: a
start from h enters pair (h,j) w.p. V(h,j) and sets h := j, while phase moves
and completions leave it alone. For a renewal MAP the chain reduces in law to
PH/c; single-server stations are left untouched.

The pair law is stored as the lifted MAP D0p (conditioned moves) and
D1p((i,j),(j,j')) = D1(i,j)/V(i,j)*V(j,j'), equivalent in law to the original.
The bookkeeping read by the state handlers is sn.ctmcmapc[(ist, r)] with keys
p, pairs (T x 2, 0-based), V and done.
"""

import numpy as np

from ....constants import ProcessType
from ....lang.base import SchedStrategy, NodeType


def _sched_of(sn, ist):
    sched = sn.sched.get(ist) if isinstance(sn.sched, dict) else sn.sched[ist]
    return sched


def ctmc_mapc(sn):
    """Rewrite the MAP service of every multiserver FCFS station into pair form (in place)."""
    M = int(sn.nstations)
    R = int(sn.nclasses)
    mapc = getattr(sn, 'ctmcmapc', None)
    if not mapc:
        mapc = {}
    old_phasessz = np.array(sn.phasessz, dtype=int, copy=True)
    old_phaseshift = np.array(sn.phaseshift, dtype=int, copy=True)
    changed = []
    for ist in range(M):
        ind = int(sn.stationToNode[ist])
        nt = sn.nodetype[ind]
        if nt == NodeType.SOURCE or _sched_of(sn, ist) != SchedStrategy.FCFS:
            continue
        c = float(np.asarray(sn.nservers).ravel()[ist])
        if not np.isfinite(c) or c <= 1:
            continue
        for r in range(R):
            if (ist, r) in mapc or sn.procid[ist, r] not in (ProcessType.MAP, ProcessType.MMPP2):
                continue
            D0 = np.atleast_2d(np.asarray(sn.proc[ist][r][0], dtype=float))
            D1 = np.atleast_2d(np.asarray(sn.proc[ist][r][1], dtype=float))
            p = D0.shape[0]
            V = np.linalg.solve(-D0, D1)
            V[np.abs(V) < 1e-14] = 0.0
            pairs = np.array([(i, j) for i in range(p) for j in range(p) if V[i, j] > 0], dtype=int)
            T = pairs.shape[0]
            tp = -np.ones((p, p), dtype=int)
            for t in range(T):
                tp[pairs[t, 0], pairs[t, 1]] = t
            done = np.zeros(T)
            D0p = np.zeros((T, T))
            D1p = np.zeros((T, T))
            for t in range(T):
                i, j = pairs[t]
                done[t] = D1[i, j] / V[i, j]
                D0p[t, t] = D0[i, i]
                for k in range(p):
                    if k != i and D0[i, k] != 0 and tp[k, j] >= 0:
                        D0p[t, tp[k, j]] = D0[i, k] * V[k, j] / V[i, j]
                for jn in range(p):
                    if tp[j, jn] >= 0:
                        D1p[t, tp[j, jn]] = done[t] * V[j, jn]
            from ...mam.map_analysis import map_pie
            sn.proc[ist][r] = [D0p, D1p]
            sn.pie[ist][r] = np.asarray(map_pie(D0p, D1p), dtype=float).ravel()
            sn.mu[ist][r] = -np.diag(D0p)
            sn.phi[ist][r] = done / (-np.diag(D0p))
            sn.phases[ist, r] = T
            sn.phasessz[ist, r] = T
            mapc[(ist, r)] = {'p': p, 'pairs': pairs, 'V': V, 'done': done}
            changed.append((ist, r))
    if not changed:
        return sn
    sn.ctmcmapc = mapc
    for ist in sorted(set(i for i, _ in changed)):
        cum = 0
        for r in range(R):
            sn.phaseshift[ist, r] = cum
            cum += int(sn.phasessz[ist, r])
        if sn.phaseshift.shape[1] > R:  # (M, R+1) layout: trailing column is the total
            sn.phaseshift[ist, R] = cum
        _rebuild_state(sn, ist, old_phasessz, old_phaseshift, mapc)
    return sn


def _rebuild_state(sn, ist, old_phasessz, old_phaseshift, mapc):
    """Map the declared initial row of station ist to the pair layout.

    A job in MAP phase i is placed in the first pair (i,j); the buffer and the local
    variables, the carried phase included, are unchanged. The row is [buf | srv | vars]
    with one trailing memory variable per MAP class; a row of any other width is left
    alone, which only costs the initial-state pruning.
    """
    if getattr(sn, 'state', None) is None or getattr(sn, 'stationToStateful', None) is None:
        return
    isf = int(np.asarray(sn.stationToStateful).ravel()[ist])
    try:
        old = np.atleast_2d(np.asarray(sn.state[isf], dtype=float))
    except (IndexError, KeyError, TypeError):
        return
    R = int(sn.nclasses)
    Kold = int(np.sum(old_phasessz[ist, :]))
    nmem = sum(1 for r in range(R) if sn.procid[ist, r] in (ProcessType.MAP, ProcessType.MMPP2)
               and int(old_phasessz[ist, r]) > 1)
    nbuf = old.shape[1] - Kold - nmem
    if nbuf < 0:
        return
    rows = []
    for row in old:
        buf = row[:nbuf]
        srv_old = row[nbuf:nbuf + Kold]
        var = row[nbuf + Kold:]
        srv = np.zeros(int(np.sum(sn.phasessz[ist, :])))
        for r in range(R):
            blk = srv_old[old_phaseshift[ist, r]:old_phaseshift[ist, r] + old_phasessz[ist, r]]
            base = int(sn.phaseshift[ist, r])
            if (ist, r) not in mapc:
                srv[base:base + int(sn.phasessz[ist, r])] = blk
                continue
            pairs = mapc[(ist, r)]['pairs']
            for i in np.flatnonzero(blk > 0):
                t = int(np.flatnonzero(pairs[:, 0] == i)[0])
                srv[base + t] += blk[i]
        rows.append(np.concatenate([buf, srv, var]))
    new = np.array(rows)
    sn.state[isf] = new[0] if np.ndim(sn.state[isf]) == 1 else new
