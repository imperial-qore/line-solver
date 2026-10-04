"""
Synchronous call DAG carrying reference-task customers into a layer.

Native Python port of ``matlab/src/api/lqn/lqn_ref_routes.m``. Indices are the
0-based ones of the native LayeredNetworkStruct; a position absent from a group
is -1 where MATLAB writes 0.

Resolves, for the caller set of one layer, the synchronous call graph along which
reference (REF) task customers descend to those callers, and the mean number of
times each entry and each call is invoked per REF cycle. SolverLN gives every
caller of a layer its own client chain; the caller set of one group, ``members``,
is what the 'refpath' interlock method turns into ONE chain.

Nodes are ENTRIES. Visits are computed topologically, v(u) = sum over parents of
v(p)*w(a)*callmean, and routes are COUNTED in the same pass, never enumerated.
Only SYNC calls are followed: an ASYNC call terminates blocking, and forwarding
is already flattened into pseudo-SYNC arcs before the layers are built.

Each group carries:
    reftask       task index of the REF task at the root
    headIsCaller  true when the REF task is a caller of this layer
    members       callers of this layer on the DAG, descent order
    entries       every DAG entry, topologically ordered, root entries first
    etask         lqn.parent of each entry
    ismember      true where that entry's task is a caller of this layer
    vEntry        mean invocations of each entry per REF cycle
    actweight     list, per entry, of a (2, k) array [aidx; executions per invocation]
    calls         (n, 5) rows [cidx, fromPos, toPos, aidx, vCall], vCall per REF cycle
    prefixPos     positions in entries forming the prefix, topological order
    prefixTerm    true where that prefix position is a first caller
    npaths        distinct REF-to-caller routes, counted
    poolmin       min of lqn.maxmult over the DAG tasks, a DIAGNOSTIC only

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from types import SimpleNamespace
from typing import List, Optional, Sequence, Tuple

import numpy as np

from ...constants import CallType


def lqn_ref_routes(lqn, callers: Sequence[int], maxpaths: Optional[float] = None,
                   serverSet: Optional[Sequence[int]] = None) -> Tuple[List[SimpleNamespace], str]:
    """
    Synchronous call DAG carrying reference-task customers into a layer.

    Args:
        lqn: LayeredNetworkStruct (0-based indices).
        callers: Task indices that call the layer's server.
        maxpaths: Refuse the layer above this many reference-task routes into it. Default 32.
        serverSet: Server elements of the layer. A prefix node whose task is one of them is
            recursion and refuses the layer.

    Returns:
        (R, why): R is a list with one group per reference task reaching CALLERS; WHY is
        non-empty when the layer must fall back to another interlock method, R being empty.
    """
    R: List[SimpleNamespace] = []
    why = ''
    if maxpaths is None:
        maxpaths = 32
    server_set = [] if serverSet is None else [int(s) for s in np.atleast_1d(serverSet)]
    ncalls = int(getattr(lqn, 'ncalls', 0) or 0)
    callers = [int(c) for c in np.atleast_1d(callers)] if callers is not None else []
    if ncalls == 0 or len(callers) == 0:
        return R, why

    callers = sorted(set(callers))
    is_caller = np.zeros(lqn.nidx, dtype=bool)
    is_caller[callers] = True

    succ = _sync_successors(lqn)

    reftasks = []
    for t in range(lqn.ntasks):
        tidx = lqn.tshift + t
        if _isref(lqn, tidx):
            reftasks.append(tidx)
    if not reftasks:
        return R, why

    # A caller reachable from two REF tasks is two INDEPENDENT customer pools; merging them would invent a
    # correlation. The whole LAYER falls back, since a refused caller may lie on another group's path.
    nref_of = np.zeros(lqn.nidx, dtype=int)
    for r in reftasks:
        seen = _reachable_entries(lqn, succ, r)
        mem = set(int(_parent(lqn, e)) for e in np.flatnonzero(seen))
        mem = set(m for m in mem if m >= 0 and is_caller[m])
        if is_caller[r]:
            mem.add(r)
        for m in mem:
            nref_of[m] += 1
    overloaded = np.flatnonzero(nref_of > 1)
    if overloaded.size > 0:
        o = int(overloaded[0])
        why = ("task '%s' is reachable from %d reference tasks, whose customer pools are independent"
               % (_name_of(lqn, o), int(nref_of[o])))
        return R, why

    for r in reftasks:
        grp, gwhy = _build_group(lqn, succ, r, is_caller, maxpaths, server_set)
        if gwhy:
            return [], gwhy
        if grp is not None:
            R.append(grp)
    return R, why


def _isref(lqn, idx: int) -> bool:
    isref = np.asarray(lqn.isref).ravel()
    return idx < isref.size and bool(isref[idx])


def _parent(lqn, idx: int) -> int:
    par = np.asarray(lqn.parent).ravel()
    return int(par[idx])


def _entries_of(lqn, tidx: int) -> List[int]:
    eo = lqn.entriesof.get(tidx, []) if isinstance(lqn.entriesof, dict) else lqn.entriesof[tidx]
    return [int(e) for e in np.atleast_1d(eo)] if eo is not None else []


def _acts_of(lqn, eidx: int) -> List[int]:
    ao = lqn.actsof.get(eidx, []) if isinstance(lqn.actsof, dict) else lqn.actsof[eidx]
    return [int(a) for a in np.atleast_1d(ao)] if ao is not None else []


def _sync_successors(lqn):
    """Per calling ENTRY, rows [cidx, called entry, calling activity, callmean] of its SYNC calls.

    The calling entry is NOT lqn.parent of the activity (that is its TASK), so actsof is inverted over the
    ENTRY range instead.
    """
    succ = [[] for _ in range(lqn.nidx)]
    entry_of_act = -np.ones(lqn.nidx, dtype=int)
    for e in range(lqn.nentries):
        eidx = lqn.eshift + e
        for a in _acts_of(lqn, eidx):
            entry_of_act[a] = eidx
    calltype = np.asarray(lqn.calltype).ravel()
    callpair = np.asarray(lqn.callpair)
    for cidx in range(lqn.ncalls):
        if int(calltype[cidx]) != int(CallType.SYNC):
            continue
        aidx = int(callpair[cidx, 0])
        eidx_to = int(callpair[cidx, 1])
        if aidx < 0 or eidx_to < 0:
            continue
        eidx_from = int(entry_of_act[aidx])
        if eidx_from < 0:
            continue
        succ[eidx_from].append((cidx, eidx_to, aidx, _call_mean_of(lqn, cidx)))
    return succ


def _call_mean_of(lqn, cidx: int) -> float:
    """Mean number of invocations carried by call CIDX."""
    m = float('nan')
    cpm = getattr(lqn, 'callproc_mean', None)
    if cpm is not None and np.size(cpm) > cidx:
        m = float(np.asarray(cpm).ravel()[cidx])
    if not np.isfinite(m):
        cp = getattr(lqn, 'callproc', None)
        proc = None
        if isinstance(cp, dict):
            proc = cp.get(cidx)
        elif cp is not None and len(cp) > cidx:
            proc = cp[cidx]
        if proc is not None:
            if hasattr(proc, 'getMean'):
                m = float(proc.getMean())
            elif hasattr(proc, 'mean'):
                m = float(proc.mean)
    if not np.isfinite(m):
        callpair = getattr(lqn, 'callpair', None)
        if callpair is not None and np.asarray(callpair).shape[0] > cidx and np.asarray(callpair).shape[1] > 2:
            m = float(np.asarray(callpair)[cidx, 2])
    if not np.isfinite(m):
        m = 1.0
    return m


def _name_of(lqn, idx: int) -> str:
    hn = getattr(lqn, 'hashnames', None)
    if hn is not None and len(hn) > idx and hn[idx] is not None and str(hn[idx]) != '':
        return str(hn[idx])
    return '#%d' % idx


def _reachable_entries(lqn, succ, tidx: int) -> np.ndarray:
    """Every entry reachable from task TIDX over SYNC calls, without pruning."""
    seen = np.zeros(lqn.nidx, dtype=bool)
    stack = list(_entries_of(lqn, tidx))
    while stack:
        eidx = stack.pop()
        if seen[eidx]:
            continue
        seen[eidx] = True
        for row in succ[eidx]:
            stack.append(row[1])
    return seen


def _build_group(lqn, succ, r: int, is_caller, maxpaths, server_set):
    """One group: the synchronous DAG below REF task R, its per-cycle visits, and the prefix above the callers."""
    roots = _entries_of(lqn, r)
    if not roots:
        return None, ''

    # Depth-first sweep with an on-stack marker. A back edge is REFUSED: a cycle makes v(u) a geometric series
    # the chain would have to express as a self-loop through the layer's own server.
    WHITE, GREY, BLACK = 0, 1, 2
    color = np.zeros(lqn.nidx, dtype=int)
    post = []
    for e0 in roots:
        if color[e0] != WHITE:
            continue
        stack = [[e0, 0]]
        while stack:
            u, ci = stack[-1]
            if ci == 0:
                color[u] = GREY
            kids = [row[1] for row in succ[u]]
            if ci < len(kids):
                stack[-1][1] = ci + 1
                v = kids[ci]
                if color[v] == GREY:
                    return None, "the synchronous call graph below '%s' is cyclic" % _name_of(lqn, v)
                elif color[v] == WHITE:
                    stack.append([v, 0])
            else:
                color[u] = BLACK
                post.append(u)
                stack.pop()

    entries = post[::-1]  # topological order
    if not entries:
        return None, ''
    pos = -np.ones(lqn.nidx, dtype=int)
    for i, e in enumerate(entries):
        pos[e] = i
    etask = [_parent(lqn, e) for e in entries]
    ismem = np.array([t >= 0 and bool(is_caller[t]) for t in etask], dtype=bool)

    if not ismem.any():
        return None, ''  # this REF reaches none of the callers

    actweight = []
    calls = []
    for i, u in enumerate(entries):
        aw, awwhy = _entry_act_weights(lqn, u)
        if awwhy:
            return None, awwhy
        actweight.append(aw)
        for (cidx, v, aidx, cmean) in succ[u]:
            if pos[v] < 0:
                continue
            w = 0.0
            if aw.shape[1] > 0:
                hit = np.flatnonzero(aw[0, :] == aidx)
                if hit.size > 0:
                    w = float(aw[1, hit[0]])
            calls.append([cidx, i, int(pos[v]), aidx, w * cmean])
    calls = np.array(calls, dtype=float).reshape(-1, 5)

    # Topological visits: the REF task selects among its entries with equal probability.
    n = len(entries)
    vEntry = np.zeros(n)
    rootpos = [int(pos[e]) for e in roots if pos[e] >= 0]
    for rp in rootpos:
        vEntry[rp] = 1.0 / len(roots)
    for i in range(n):
        for k in np.flatnonzero(calls[:, 1] == i):
            j = int(calls[k, 2])
            vEntry[j] = vEntry[j] + vEntry[i] * calls[k, 4]
    vCall = np.array([vEntry[int(calls[k, 1])] * calls[k, 4] for k in range(calls.shape[0])])
    if calls.shape[0] > 0:
        calls[:, 4] = vCall

    # The prefix: entries reachable from a root without passing THROUGH a caller.
    inPrefix = np.zeros(n, dtype=bool)
    prefTerm = np.zeros(n, dtype=bool)
    npath = np.zeros(n)
    for rp in rootpos:
        inPrefix[rp] = True
        npath[rp] = 1
    for i in range(n):
        if not inPrefix[i]:
            continue
        if ismem[i]:
            prefTerm[i] = True
            continue  # do not descend past a caller
        # A hop whose task is a server of this layer would put the customer in two places at once.
        if server_set and etask[i] in server_set:
            return None, ("task '%s' is both an intermediate on the reference path and a server of this layer"
                          % _name_of(lqn, etask[i]))
        for k in np.flatnonzero(calls[:, 1] == i):
            j = int(calls[k, 2])
            inPrefix[j] = True
            npath[j] = npath[j] + npath[i]

    np_routes = float(np.sum(npath[prefTerm]))
    if np_routes > maxpaths:
        return None, ('the reference path into this layer carries %g distinct routes, '
                      'above config.interlock_maxpaths = %g' % (np_routes, maxpaths))

    members = []
    for t, m in zip(etask, ismem):
        if m and t not in members:
            members.append(t)

    g = SimpleNamespace(
        reftask=r,
        headIsCaller=bool(is_caller[r]),
        members=members,
        entries=list(entries),
        etask=list(etask),
        ismember=ismem,
        vEntry=vEntry,
        actweight=actweight,
        calls=calls,
        prefixPos=[int(i) for i in np.flatnonzero(inPrefix)],
        prefixTerm=prefTerm[inPrefix],
        npaths=np_routes,
        poolmin=_pool_min(lqn, etask),
    )
    return g, ''


def _pool_min(lqn, tasks) -> float:
    """Smallest thread pool on the path. DIAGNOSTIC ONLY: the chain population is never capped by it."""
    p = float('inf')
    mm = getattr(lqn, 'maxmult', None)
    if mm is None:
        return p
    mm = np.asarray(mm).ravel()
    for t in sorted(set(int(t) for t in tasks if t >= 0)):
        if mm.size > t:
            p = min(p, float(mm[t]))
    return p


def _entry_act_weights(lqn, eidx: int):
    """(2, k) array [aidx; mean executions per invocation of entry EIDX], propagated over lqn.graph."""
    acts = _acts_of(lqn, eidx)
    if not acts:
        return np.zeros((2, 0)), ''
    nodeset = [eidx] + acts
    pos = {u: i for i, u in enumerate(nodeset)}
    m = len(nodeset)
    graph = lqn.graph
    A = np.zeros((m, m))
    for i, u in enumerate(nodeset):
        row = np.asarray(graph[u, :]).ravel()
        for v in np.flatnonzero(row):
            if int(v) in pos:
                A[i, pos[int(v)]] = float(row[v])

    # Topological propagation with an explicit cycle test: an activity-graph loop makes the executions a
    # geometric series the chain cannot express, so the layer falls back rather than truncating the count.
    remaining = np.sum(A > 0, axis=0).astype(int)
    remaining[0] = 0  # the entry is the source
    w = np.zeros(m)
    w[0] = 1.0
    queue = [int(i) for i in np.flatnonzero(remaining == 0)]
    done = np.zeros(m, dtype=bool)
    ndone = 0
    while queue:
        i = queue.pop(0)
        if done[i]:
            continue
        done[i] = True
        ndone += 1
        for j in np.flatnonzero(A[i, :] > 0):
            w[j] = w[j] + w[i] * A[i, j]
            remaining[j] -= 1
            if remaining[j] <= 0 and not done[j]:
                queue.append(int(j))
    if ndone < m:
        return np.zeros((2, 0)), "the activity graph of entry '%s' contains a loop" % _name_of(lqn, eidx)
    aw = np.vstack([np.array(acts, dtype=float), w[1:]])
    return aw, ''
