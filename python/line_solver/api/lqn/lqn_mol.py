"""Method of Layers on the SRVN decomposition of an entry-only LQN.

Port of ``matlab/src/api/lqn/lqn_mol.m``. A compact, self-contained
reimplementation of the layered fixed point ``SolverLN`` runs, restricted to
LQNs in which every entry binds exactly one activity and there are no activity
precedences. It decomposes the model the way ``lqns --srvn-layering`` does --
one submodel per processor and one per called task -- and sweeps them in the two
phases of Rolia-Sevcik's Method of Layers: all software (task) submodels, then
all hardware (processor) ones.

Every submodel is a closed multiclass queueing network with ONE station and one
class per client task, so it is solved by :func:`pfqn_qdamva` rather than by
building a ``Network``. The surrogate client delay of ``SolverLN`` collapses into
the think-time vector ``Z`` of that call.

WHERE THIS DIFFERS FROM ``SolverLN``'s ``srvn.cs``: a submodel here carries one
class per client TASK with visit-weighted demands, where ``srvn.cs`` carries one
class per activity and encodes the call multiplicities as routing. On an
entry-only model the two agree on the structure and differ only in the
aggregation, so the throughputs and processor utilizations track closely while
entry response times spread more.

SCOPE. Entry-only models. Activity graphs (fork/join, OR-branches, loops, second
phases, forwarding), asynchronous calls, caches, setup tasks, admission
constraints, replication and open arrivals are REFUSED, not approximated, and
named when they are.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import math
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ...constants import SchedStrategy
from ..pfqn.qdamva import pfqn_qdamva

__all__ = ['lqn_mol']

# GlobalConstants.FineTol, the reference's own "effectively zero".
_FINE_TOL = 1e-8

_CALL_KIND = {1: 'synchronous', 2: 'asynchronous', 3: 'forwarding'}

# THE PYTHON lsn IS 0-BASED, unlike MATLAB's and the C++ port's, which number
# from 1 and leave slot 0 unused: here hosts are 0..nhosts-1, tasks
# tshift..tshift+ntasks-1, and so on, with nidx the COUNT rather than the
# largest index. Every loop below is written in that convention; transcribing
# the reference's 1-based ranges instead shifts every element by one and reads
# an activity where an entry was meant, which is what the assertions catch
# first and the numbers would not.


def _sched_of(lsn, idx: int):
    """The scheduling strategy of element `idx`, as a SchedStrategy."""
    s = lsn.sched.get(idx)
    if isinstance(s, SchedStrategy):
        return s
    try:
        return SchedStrategy(int(s))
    except (TypeError, ValueError):
        return s


def _sched_name(s: Any) -> str:
    return getattr(s, 'name', str(s))


def _is_sched(s: Any, name: str) -> bool:
    return _sched_name(s).upper() == name


def _vec(a, n: int) -> np.ndarray:
    """A (1 x m) or nested field flattened to a plain length-m vector."""
    if a is None:
        return np.zeros(n)
    return np.asarray(a, dtype=float).ravel()


def _mean_of(d) -> float:
    """The mean of a hostdem/think entry, which is a float or a Distribution."""
    if d is None:
        return 0.0
    if isinstance(d, (int, float, np.floating)):
        return float(d)
    for attr in ('getMean', 'get_mean'):
        if hasattr(d, attr):
            return float(getattr(d, attr)())
    return 0.0


def _mol_assert(lsn) -> None:
    """Refuse every feature this decomposition does not represent, NAMING the
    element, rather than returning a number that quietly ignores it."""
    hn = lsn.hashnames
    for eidx in range(lsn.eshift, lsn.eshift + lsn.nentries):
        acts = lsn.actsof.get(eidx, [])
        if len(acts) != 1:
            raise ValueError(
                'lqn_mol: entry %s binds %d activities. lqn_mol solves entry-only models; '
                'use SolverLN for an activity graph.' % (hn[eidx], len(acts)))
    graph = np.asarray(lsn.graph)
    ph = _vec(getattr(lsn, 'actphase', None), lsn.nacts)
    for a in range(lsn.nacts):
        aidx = lsn.ashift + a
        # A successor INSIDE the activity band is a precedence; an edge to an
        # entry or a task is the ordinary binding every entry-only model has.
        row = graph[aidx, lsn.ashift:lsn.ashift + lsn.nacts]
        if np.any(row != 0):
            raise ValueError(
                'lqn_mol: activity %s has an activity precedence. lqn_mol solves entry-only '
                'models; use SolverLN for an activity graph.' % hn[aidx])
        if a < ph.size and int(ph[a]) != 1:
            raise ValueError('lqn_mol: activity %s is in phase %d. lqn_mol supports phase 1 only.'
                             % (hn[aidx], int(ph[a])))
    for cidx in range(lsn.ncalls):
        ct = int(lsn.calltype[cidx])
        if ct != 1:
            name = lsn.callhashnames[cidx] if getattr(lsn, 'callhashnames', None) is not None \
                else 'call %d' % (cidx + 1)
            raise ValueError('lqn_mol: call %s is %s. lqn_mol supports synchronous calls only.'
                             % (name, _CALL_KIND.get(ct, 'none')))
    nelem = lsn.tshift + lsn.ntasks
    cache = _vec(getattr(lsn, 'iscache', None), nelem)
    setup = _vec(getattr(lsn, 'hassetup', None), nelem)
    repl = _vec(getattr(lsn, 'repl', None), nelem)
    lincon = getattr(lsn, 'lincon', None) or {}
    for idx in range(nelem):
        if idx < cache.size and cache[idx]:
            raise ValueError('lqn_mol: %s is a cache task, which lqn_mol does not model.' % hn[idx])
        if idx < setup.size and setup[idx]:
            raise ValueError('lqn_mol: %s has a setup time, which lqn_mol does not model.' % hn[idx])
        if idx < repl.size and repl[idx] != 1.0:
            raise ValueError('lqn_mol: %s is replicated %d times, which lqn_mol does not model.'
                             % (hn[idx], int(repl[idx])))
        if lincon.get(idx) is not None:
            raise ValueError('lqn_mol: %s carries an admission constraint, which lqn_mol does not '
                             'model.' % hn[idx])
        s = _sched_of(lsn, idx)
        if idx < lsn.nhosts:
            ok = any(_is_sched(s, n) for n in ('PS', 'FCFS', 'INF'))
        else:
            ok = any(_is_sched(s, n) for n in ('PS', 'FCFS', 'INF', 'REF'))
        if not ok:
            raise ValueError('lqn_mol: %s is scheduled %s, which lqn_mol does not model.'
                             % (hn[idx], _sched_name(s)))
    if getattr(lsn, 'callgroups', None):
        raise ValueError('lqn_mol: this model uses routed call groups, which lqn_mol does not '
                         'model.')
    arr = getattr(lsn, 'arrival', None) or {}
    for eidx in range(lsn.eshift, lsn.eshift + lsn.nentries):
        if arr.get(eidx) is not None:
            raise ValueError('lqn_mol: entry %s has an open arrival. lqn_mol solves closed models '
                             'only.' % hn[eidx])


def _mol_mu(N: np.ndarray, c: float) -> np.ndarray:
    """The queue-dependent rate multiplier row of a c-server station over a
    population of ``sum(N)``.

    ``pfqn_lldfun`` SKIPS a constant row, so a single server must come back as a
    row of ones and not as a scalar 1, or the multiserver term is never applied.
    """
    smax = max(2, int(math.ceil(max(0.0, float(np.sum(N))))))
    if not math.isfinite(c) or c <= 1.0:
        return np.ones((1, smax))
    return np.minimum(np.arange(1, smax + 1, dtype=float), c).reshape(1, smax)


def lqn_mol(lsn, options: Optional[Dict[str, Any]] = None
            ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]:
    """Method of Layers on an entry-only, closed, synchronous LQN.

    Args:
        lsn: LayeredNetworkStruct, from ``LayeredNetwork.getStruct()``.
        options: optional dict with ``iter_max`` (200), ``iter_tol`` (1e-6) and
            ``relax_factor`` (0.5).

    Returns:
        ``(QN, UN, RN, TN, info)``, four ``(nidx+1,)`` vectors in the column
        convention ``SolverLN`` and LQNS report, so they line up with
        ``LN(model).get_avg_table()`` cell for cell:

        =======  ================  =====================  ==============  ============
        index    QN (QLen)         UN (Util)              RN (RespT)      TN (Tput)
        =======  ================  =====================  ==============  ============
        host     NaN               processor utilization  NaN             NaN
        task     sum of entry T*S  sum of entry proc util NaN             cycle rate
        entry    T*S               processor utilization  response time   throughput
        act      as its entry      as its entry           as its entry    as its entry
        =======  ================  =====================  ==============  ============
    """
    options = dict(options or {})
    iter_max = int(options.get('iter_max', 200))
    iter_tol = float(options.get('iter_tol', 1e-6))
    om = float(options.get('relax_factor', 0.5))

    _mol_assert(lsn)

    nidx = lsn.nidx
    eidxs = list(range(lsn.eshift, lsn.eshift + lsn.nentries))
    tidxs = list(range(lsn.tshift, lsn.tshift + lsn.ntasks))
    hidxs = list(range(lsn.nhosts))
    ncalls = lsn.ncalls
    parent = _vec(lsn.parent, nidx)
    isref = _vec(lsn.isref, nidx) != 0

    # ---- static per-entry data ------------------------------------------------
    actof = np.zeros(nidx, dtype=int)
    dem = np.zeros(nidx)
    taskof = np.zeros(nidx, dtype=int)
    hostof = np.zeros(nidx, dtype=int)
    for tidx in tidxs:
        hostof[tidx] = int(parent[tidx])
    for eidx in eidxs:
        actof[eidx] = int(lsn.actsof[eidx][0])
        d = _mean_of(lsn.hostdem.get(actof[eidx]))
        dem[eidx] = 0.0 if math.isnan(d) else d
        taskof[eidx] = int(parent[eidx])

    # ---- static per-call data -------------------------------------------------
    callsrc = np.zeros(ncalls, dtype=int)
    calldst = np.zeros(ncalls, dtype=int)
    cally = np.zeros(ncalls)
    entryOfAct = np.zeros(nidx, dtype=int)
    for eidx in eidxs:
        entryOfAct[actof[eidx]] = eidx
    callsFrom: List[List[int]] = [[] for _ in range(nidx)]
    callsTo: List[List[int]] = [[] for _ in range(nidx)]
    for cidx in range(ncalls):
        callsrc[cidx] = entryOfAct[int(lsn.callpair[cidx, 0])]
        calldst[cidx] = int(lsn.callpair[cidx, 1])
        y = float(lsn.callpair[cidx, 2])
        cally[cidx] = 0.0 if math.isnan(y) else y
        callsFrom[callsrc[cidx]].append(cidx)
        callsTo[calldst[cidx]].append(cidx)

    # ---- populations, from maxmult (mult is wrong for INF tasks) --------------
    npop = np.ones(nidx)
    maxmult = _vec(getattr(lsn, 'maxmult', None), nidx)
    for idx in range(lsn.tshift + lsn.ntasks):
        m = float(maxmult[idx]) if idx < maxmult.size else 1.0
        if not math.isfinite(m) or m < 1.0:
            m = 1.0
        npop[idx] = m

    # ---- layer sets -----------------------------------------------------------
    # One hardware layer per populated host, one software layer per called
    # non-reference task, as buildLayers.m draws them.
    hostLayers = [h for h in hidxs if lsn.tasksof.get(h)]
    isCalled = np.zeros(nidx, dtype=bool)
    for cidx in range(ncalls):
        isCalled[taskof[calldst[cidx]]] = True
    taskLayers = [t for t in tidxs if not isref[t] and isCalled[t]]

    # ---- fixed-point state ----------------------------------------------------
    residt = dem.copy()
    servt = np.zeros(nidx)
    callservt = np.zeros(ncalls)
    thinkt = np.zeros(nidx)
    share = np.zeros(nidx)
    Xtask = np.zeros(nidx)
    Xentry = np.zeros(nidx)
    busyth = np.zeros(nidx)
    zref = np.zeros(nidx)
    for tidx in tidxs:
        # lqn_ref_thinktime: a reference task's declared think time, and zero for
        # every other task and for a negative or non-finite one.
        z = 0.0
        if isref[tidx]:
            z = _mean_of(lsn.think.get(tidx))
            if not math.isfinite(z) or z < 0.0:
                z = 0.0
        zref[tidx] = z
        thinkt[tidx] = z
        es = lsn.entriesof.get(tidx, [])
        if es:
            share[es] = 1.0 / len(es)
    # Seed servt bottom-up over the call graph so a callee is priced before its
    # caller; a cycle just leaves the residual demand seeded at 0.
    servt[eidxs] = dem[eidxs]
    for _ in range(max(1, lsn.nentries)):
        for eidx in eidxs:
            s = residt[eidx]
            for cidx in callsFrom[eidx]:
                s += cally[cidx] * servt[calldst[cidx]]
            servt[eidx] = s
    for cidx in range(ncalls):
        callservt[cidx] = servt[calldst[cidx]]

    def cycle_outside(tidx: int, excl: int) -> float:
        """Time a thread of ``tidx`` spends away from ``excl`` in one cycle: its
        think time, its own processor residence, and its blocking at every
        callee other than ``excl``."""
        z = thinkt[tidx]
        for eidx in lsn.entriesof.get(tidx, []):
            w = share[eidx] * residt[eidx]
            for cidx in callsFrom[eidx]:
                if taskof[calldst[cidx]] != excl:
                    w += share[eidx] * cally[cidx] * callservt[cidx]
            z += w
        return z

    def host_servers(hidx: int) -> float:
        """LINE scales a station utilization into [0,1] whatever its
        multiplicity, and reports busy SERVERS at an infinite server."""
        return 1.0 if _is_sched(_sched_of(lsn, hidx), 'INF') else npop[hidx]

    def throughputs() -> None:
        """Reference tasks set the pace; every other rate follows from the call
        rates, so the entries are visited in call-graph order until stable."""
        for t in tidxs:
            if isref[t]:
                cyc = thinkt[t]
                for eidx in lsn.entriesof.get(t, []):
                    cyc += share[eidx] * servt[eidx]
                Xtask[t] = npop[t] / cyc if cyc > _FINE_TOL else 0.0
                for eidx in lsn.entriesof.get(t, []):
                    Xentry[eidx] = Xtask[t] * share[eidx]
            else:
                Xtask[t] = 0.0
                for eidx in lsn.entriesof.get(t, []):
                    Xentry[eidx] = 0.0
        for _ in range(max(1, lsn.ntasks)):
            for eidx in eidxs:
                if isref[taskof[eidx]]:
                    continue
                x = 0.0
                for cidx in callsTo[eidx]:
                    x += cally[cidx] * Xentry[callsrc[cidx]]
                Xentry[eidx] = x
            for t in tidxs:
                if isref[t]:
                    continue
                es = lsn.entriesof.get(t, [])
                Xtask[t] = float(np.sum(Xentry[es])) if es else 0.0
                if Xtask[t] > _FINE_TOL:
                    share[es] = Xentry[es] / Xtask[t]
        # Mean busy threads, by Little's law over the entries the task serves.
        # This is the occupancy the think-time closure needs, and it is exact
        # given the throughputs -- unlike the layer AMVA's own U, which is X*L*g
        # with g a reciprocal rate multiplier and not a server count.
        for t in tidxs:
            u = 0.0
            for eidx in lsn.entriesof.get(t, []):
                u += Xentry[eidx] * servt[eidx]
            busyth[t] = u

    resid = float('inf')
    it = 0
    while it < iter_max:
        it += 1
        servt_prev = servt.copy()
        thinkt_prev = thinkt.copy()

        # ---- phase 1: software layers (thread contention at each called task)
        for tidx in taskLayers:
            # The task is the station, its caller tasks the classes.
            callers = sorted({taskof[callsrc[cidx]]
                              for eidx in lsn.entriesof.get(tidx, [])
                              for cidx in callsTo[eidx]})
            K = len(callers)
            if K == 0:
                continue
            gcl = np.ones(K)
            if not _is_sched(_sched_of(lsn, tidx), 'INF'):
                # An infinite-thread task never queues for a thread.
                L = np.zeros((1, K))
                N = np.zeros(K)
                Z = np.zeros(K)
                for k, ctask in enumerate(callers):
                    N[k] = npop[ctask]
                    d = 0.0
                    for eidx in lsn.entriesof.get(ctask, []):
                        for cidx in callsFrom[eidx]:
                            if taskof[calldst[cidx]] == tidx:
                                d += share[eidx] * cally[cidx] * servt[calldst[cidx]]
                    L[0, k] = d
                    Z[k] = cycle_outside(ctask, tidx)
                # The AMVA U output is X*L*g, where g is the RECIPROCAL RATE
                # MULTIPLIER at the current congestion, not 1/c -- it is not a
                # busy-server count, so nothing here reads it. Occupancy comes
                # from Little's law in throughputs().
                _, _, _, _, R = pfqn_qdamva(L, N, Z, _mol_mu(N, npop[tidx]))
                for k in range(K):
                    if L[0, k] > _FINE_TOL:
                        gcl[k] = R[0, k] / L[0, k]
            for k, ctask in enumerate(callers):
                for eidx in lsn.entriesof.get(ctask, []):
                    for cidx in callsFrom[eidx]:
                        if taskof[calldst[cidx]] != tidx:
                            continue
                        newv = gcl[k] * servt[calldst[cidx]]
                        callservt[cidx] = om * newv + (1.0 - om) * callservt[cidx]

        # ---- phase 2: hardware layers (processor contention at each host) -----
        for hidx in hostLayers:
            # The processor is the station, its tasks the classes.
            tsks = lsn.tasksof.get(hidx, [])
            K = len(tsks)
            f = np.ones(K)
            if not _is_sched(_sched_of(lsn, hidx), 'INF'):
                # A delay processor never queues.
                L = np.zeros((1, K))
                N = np.zeros(K)
                Z = np.zeros(K)
                for k, tidx in enumerate(tsks):
                    N[k] = npop[tidx]
                    d = 0.0
                    z = thinkt[tidx]
                    for eidx in lsn.entriesof.get(tidx, []):
                        d += share[eidx] * dem[eidx]
                        for cidx in callsFrom[eidx]:
                            z += share[eidx] * cally[cidx] * callservt[cidx]
                    L[0, k] = d
                    Z[k] = z
                _, _, _, _, R = pfqn_qdamva(L, N, Z, _mol_mu(N, npop[hidx]))
                for k in range(K):
                    if L[0, k] > _FINE_TOL:
                        f[k] = R[0, k] / L[0, k]
            for k, tidx in enumerate(tsks):
                for eidx in lsn.entriesof.get(tidx, []):
                    residt[eidx] = f[k] * dem[eidx]

        # ---- recompose entry service times ------------------------------------
        for eidx in eidxs:
            s = residt[eidx]
            for cidx in callsFrom[eidx]:
                s += cally[cidx] * callservt[cidx]
            servt[eidx] = om * s + (1.0 - om) * servt[eidx]

        # ---- throughputs, entry shares, think-time closure ---------------------
        throughputs()
        for tidx in tidxs:
            if isref[tidx]:
                thinkt[tidx] = zref[tidx]
                continue
            if Xtask[tidx] <= _FINE_TOL:
                continue
            # Idle time of a thread per cycle. updateThinkTimes splits this into
            # an INF arm (njobs - util) and a finite arm (njobs*abs(1-util)) only
            # because LINE reports busy SERVERS at an infinite server and a busy
            # FRACTION at a finite one; carrying the count in both cases makes
            # the two arms the same expression.
            newz = max(0.0, abs(npop[tidx] - busyth[tidx]) / Xtask[tidx] - zref[tidx])
            thinkt[tidx] = om * newz + (1.0 - om) * thinkt[tidx]

        # Both halves of the state must settle: servt alone can sit still for an
        # iteration while the think times are still moving.
        resid = max(
            float(np.max(np.abs(servt[eidxs] - servt_prev[eidxs])
                         / np.maximum(1.0, np.abs(servt[eidxs])))) if eidxs else 0.0,
            float(np.max(np.abs(thinkt[tidxs] - thinkt_prev[tidxs])
                         / np.maximum(1.0, np.abs(thinkt[tidxs])))) if tidxs else 0.0)
        if resid < iter_tol:
            break
    throughputs()

    # ---- assemble the reported vectors ----------------------------------------
    QN = np.full(nidx, np.nan)
    UN = np.full(nidx, np.nan)
    RN = np.full(nidx, np.nan)
    TN = np.full(nidx, np.nan)
    for eidx in eidxs:
        hidx = hostof[taskof[eidx]]
        procutil = Xentry[eidx] * dem[eidx] / host_servers(hidx)
        QN[eidx] = Xentry[eidx] * servt[eidx]
        UN[eidx] = procutil
        RN[eidx] = servt[eidx]
        TN[eidx] = Xentry[eidx]
        aidx = actof[eidx]
        QN[aidx] = QN[eidx]
        UN[aidx] = UN[eidx]
        RN[aidx] = RN[eidx]
        TN[aidx] = TN[eidx]
    for tidx in tidxs:
        es = lsn.entriesof.get(tidx, [])
        QN[tidx] = float(np.sum(QN[es])) if es else 0.0
        UN[tidx] = float(np.sum(UN[es])) if es else 0.0
        RN[tidx] = np.nan
        TN[tidx] = Xtask[tidx]
    for hidx in hidxs:
        u = 0.0
        for tidx in lsn.tasksof.get(hidx, []):
            u += UN[tidx]
        QN[hidx] = np.nan
        UN[hidx] = u
        RN[hidx] = np.nan
        # No throughput is defined at a processor, as in LQNS.
        TN[hidx] = np.nan

    info = {
        'iter': it, 'resid': resid, 'servt': servt, 'residt': residt,
        'callservt': callservt, 'thinkt': thinkt, 'share': share,
        'hostLayers': hostLayers, 'taskLayers': taskLayers,
    }
    return QN, UN, RN, TN, info
