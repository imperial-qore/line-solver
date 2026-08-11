"""
Dedicated Next-Reaction Method (NRM) SSA solver.

Native Python port of the JAR `Solver_ssa_nrm` / `Solver_ssa_nrm_space`
handlers and their analyzers (`Solver_ssa_analyzer_nrm`,
`Solver_ssa_analyzer_nrm_space`).

The direct method (`solver_ssa_nrm`) runs Gibson & Bruck's original
Next-Reaction Method using an indexed min-heap of absolute firing times and
computes the steady-state performance indices directly during the run. It
supports INF, EXT, PS, DPS, GPS, PSPRIO, DPSPRIO, GPSPRIO, FCFS and LCFS
scheduling. FCFS/LCFS stations are simulated with an ordered per-node buffer
so that the service rate is proportional to the jobs actually in service.

Port from: java/src/main/java/jline/solvers/ssa/handlers/Solver_ssa_nrm.java
"""

import numpy as np
from collections import deque
from typing import List, Optional, Tuple

from ....constants import GlobalConstants, PollingType
# Import the enums from the same module the NetworkStruct uses, so that the
# SchedStrategy / NodeType members compare equal to sn.sched / sn.nodetype.
from ...sn import SchedStrategy, NodeType, sn_region_members
from ...state import toMarginalAggr
from ...state.polling import polling_info, polling_next, polling_budget
from .handler import SolverSSAOptions, SolverSSAReturn
from ....lang.base import ReplacementStrategy


def _np_get(npm, key, default=None):
    """Read field KEY from a cache nodeparam, which may be a dict or an object."""
    if npm is None:
        return default
    if isinstance(npm, dict):
        return npm.get(key, default)
    return getattr(npm, key, default)


# =============================================================================
# Indexed min-heap (mirror of jline.util.IndexedMinHeap)
# =============================================================================
class IndexedMinHeap:
    """
    A min-heap over a fixed set of integer keys 0..n-1, each associated with a
    double priority. Supports O(log n) update of any key by index and O(1) peek
    at the minimum. Keys are reaction indices, priorities are absolute firing
    times tau_k.
    """

    def __init__(self, n: int):
        self.n = n
        self.heap = list(range(n))   # heap[pos] = reaction index at position pos
        self.pos = list(range(n))    # pos[k]    = position of reaction k in heap
        self.key = [float('inf')] * n  # key[k]   = priority of reaction k

    def peek_min(self) -> int:
        """Return the index of the reaction with the smallest key. O(1)."""
        return self.heap[0]

    def peek_min_key(self) -> float:
        """Return the smallest key value. O(1)."""
        return self.key[self.heap[0]]

    def update(self, k: int, new_key: float) -> None:
        """Update the priority of reaction k and restore heap order. O(log n)."""
        self.key[k] = new_key
        self._sift_up(self.pos[k])
        self._sift_down(self.pos[k])

    def get_key(self, k: int) -> float:
        """Read the current priority of reaction k. O(1)."""
        return self.key[k]

    def build_heap(self) -> None:
        """Establish heap order after bulk-writing key[]. O(n)."""
        for p in range(self.n // 2 - 1, -1, -1):
            self._sift_down(p)

    def _swap(self, i: int, j: int) -> None:
        self.heap[i], self.heap[j] = self.heap[j], self.heap[i]
        self.pos[self.heap[i]] = i
        self.pos[self.heap[j]] = j

    def _sift_up(self, p: int) -> None:
        while p > 0:
            parent = (p - 1) // 2
            if self.key[self.heap[parent]] > self.key[self.heap[p]]:
                self._swap(parent, p)
                p = parent
            else:
                break

    def _sift_down(self, p: int) -> None:
        while True:
            smallest = p
            l = 2 * p + 1
            r = 2 * p + 2
            if l < self.n and self.key[self.heap[l]] < self.key[self.heap[smallest]]:
                smallest = l
            if r < self.n and self.key[self.heap[r]] < self.key[self.heap[smallest]]:
                smallest = r
            if smallest != p:
                self._swap(p, smallest)
                p = smallest
            else:
                break


# =============================================================================
# Phase slot map
#
# The state vector counts jobs per (node, class, PHASE) rather than per
# (node, class), so that phase-type service is represented exactly instead of
# being collapsed onto its mean rate. Phase counts differ per (station, class)
# via sn.phasessz, hence the explicit offset map rather than arithmetic on R.
#
# The layout is chosen so that a single-phase model is bit-for-bit the old one:
# with nph == 1 everywhere, ph_off[ind, r] == ind*R + r and therefore
# slot(ind, r, 0) == ind*R + r, exactly the flat class index the engine used
# before. Every exponential model must reproduce its previous results, which is
# the self-check for this generalization.
# =============================================================================
class _Smap:
    """
    Reverse map from a state index to its (node, class, phase), plus the
    forward offsets.

    Every consumer that used to decode a state index arithmetically
    (pos // R, pos % R) must go through this instead: with unequal phase counts
    the flat arithmetic no longer identifies the node.
    """

    __slots__ = ('node', 'cls', 'phase', 'ph_off', 'nph', 'R', 'NS', 'expanded')

    def __init__(self, sn, I: int, R: int):
        nph = np.ones((I, R), dtype=int)
        for ind in range(I):
            if sn.isstation[ind]:
                ist = int(sn.nodeToStation[ind])
                for r in range(R):
                    nph[ind, r] = max(1, int(sn.phasessz[ist, r]))
        ph_off = np.zeros((I, R), dtype=int)
        NS = 0
        for ind in range(I):
            for r in range(R):
                ph_off[ind, r] = NS
                NS += int(nph[ind, r])
        self.nph = nph
        self.ph_off = ph_off
        self.R = R
        self.NS = NS
        self.expanded = NS > I * R      # False for a purely exponential model
        self.node = np.zeros(NS, dtype=int)
        self.cls = np.zeros(NS, dtype=int)
        self.phase = np.zeros(NS, dtype=int)
        for ind in range(I):
            for r in range(R):
                for kk in range(int(nph[ind, r])):
                    self.node[ph_off[ind, r] + kk] = ind
                    self.cls[ph_off[ind, r] + kk] = r
                    self.phase[ph_off[ind, r] + kk] = kk

    def node_slots(self, ind: int) -> slice:
        """
        Every state slot of node ind. The slots of a node are contiguous by
        construction (the offsets run node-major, then class, then phase), so
        the whole node is one slice.
        """
        r_last = self.R - 1
        return slice(int(self.ph_off[ind, 0]),
                     int(self.ph_off[ind, r_last] + self.nph[ind, r_last]))

    def class_slots(self, ind: int, r: int) -> slice:
        """The state slots holding the phases of class r at node ind."""
        return slice(int(self.ph_off[ind, r]),
                     int(self.ph_off[ind, r] + self.nph[ind, r]))


def _class_counts(X, smap: _Smap, ind: int) -> np.ndarray:
    """
    Per-class populations at node ind, summing each class over its phases. The
    scheduling rate laws are class-level: they are unchanged by phase
    expansion, and only the per-phase share (see _kir_frac) is layered on top.
    """
    R = smap.R
    npop = np.zeros(R)
    for r in range(R):
        npop[r] = float(np.sum(X[smap.class_slots(ind, r)]))
    return npop


def _class_pop(X, smap: _Smap, ind: int, r: int) -> float:
    """Population of class r at node ind, summed over its phases."""
    return float(np.sum(X[smap.class_slots(ind, r)]))


def _kir_frac(X, slot: int, smap: _Smap, ind: int, r: int) -> float:
    """
    Share of its class that the job population in one phase represents: kir/nir.

    The class-level rate law is split across the class's phases in this ratio,
    which is exactly how State.afterEventStation writes every phase-aware case
    (e.g. DPS uses (kir/nir) * [class share]). For a single-phase class this is
    1 whenever the class is present, so an exponential model is unaffected.
    """
    nir = float(np.sum(X[smap.class_slots(ind, r)]))
    if nir <= 0:
        return 0.0
    return float(X[slot]) / nir


def _entry_probs(sn, jnd: int, s: int, nphjs: int) -> np.ndarray:
    """
    Entry-phase distribution of a class-s job arriving at node jnd: pie of its
    service process there. A non-station node, or a station whose process is
    absent (a disabled class), has a single phase entered with probability 1.
    """
    pentry = np.zeros(nphjs)
    if nphjs <= 1:
        pentry[0] = 1.0
        return pentry
    ist = int(sn.nodeToStation[jnd])
    p = np.atleast_1d(np.asarray(sn.pie[ist][s], dtype=float)).ravel()
    if p.size == 0 or np.all(np.isnan(p)) or float(np.nansum(p)) <= 0:
        # no entry distribution declared: enter the first phase
        pentry[0] = 1.0
        return pentry
    take = min(nphjs, p.size)
    pentry[:take] = p[:take]
    return pentry / float(np.sum(pentry))


# ordered-buffer policies share the departure rate law and differ only in promotion rule / arrival rule (LCFSPR); see _pick_from_buffer and _is_preemptive.
_BUFFERED_SCHED = (SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO,
                   SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT,
                   SchedStrategy.LCFSPR)

# Of those, the policies whose state keeps the buffer as per-class counts rather
# than as an ordered list of class ids (State.fromMarginalAndRunning).
_COUNT_BUFFERED_SCHED = (SchedStrategy.SIRO, SchedStrategy.SEPT, SchedStrategy.LEPT)

# preemptive arrival displaces the incumbent into the buffer; LCFSPI is deliberately absent (unsupported by SolverSSA/SolverCTMC).
_PREEMPTIVE_SCHED = (SchedStrategy.LCFSPR,)

# PAS/OI keep the full ordered list (no server/buffer split), oldest-first; sn.sched normalizes OI to PAS so one code path covers both.
_LIST_SCHED = (SchedStrategy.PAS,)

# list-based stations seed their initial buffer the same way as count-based ones, despite differing rate laws.
_INIT_BUFFERED_SCHED = _BUFFERED_SCHED + _LIST_SCHED

# Scheduling policies supported by the NRM solver.
_ALLOWED_SCHED = (
    SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS,
    SchedStrategy.LPS, SchedStrategy.DPS, SchedStrategy.GPS,
    SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO,
    SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO,
    SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT,
    SchedStrategy.LCFSPR, SchedStrategy.PAS, SchedStrategy.POLLING,
)


def _nrm_eligible(sn) -> bool:
    """True when the NRM engine can run this model, matching the explicit
    ``method='nrm'`` gate: every station scheduling in ``_ALLOWED_SCHED`` and
    every per-feature guard (FCR, routing, balking/reneging, phase-type service,
    cache). The NRM simulates open and closed models alike; it only lacks
    Fork/Join node handling, so those are excluded. Used to prefer the NRM on the
    ``default`` and ``parallel`` dispatch paths so they match the engine's real
    capability rather than always running the serial engine.
    """
    def _sched(ist):
        return sn.sched.get(ist) if isinstance(sn.sched, dict) else sn.sched[ist]
    if not all(_sched(ist) in _ALLOWED_SCHED for ist in range(sn.nstations)):
        return False
    from ...sn import sn_has_fork_join
    if sn_has_fork_join(sn):
        return False
    return (_fcr_nrm_ok(sn) and _routing_nrm_ok(sn) and _impatience_nrm_ok(sn)
            and _phase_nrm_ok(sn) and _cache_nrm_ok(sn))


# =============================================================================
# Polling controller helpers
# =============================================================================
# A polling station has a single server that cycles through the class buffers it
# serves. The NRM carries the controller [mode, pos, swk, ctr] in the auxiliary
# per-node buffer of the polling node (bufs[ind]); mode 0 = PARKED, 1 = SERVING
# pos, 2 = SWITCHING towards pos, swk = 1-based phase of the switchover PH, ctr =
# visit budget. This mirrors the MATLAB solver_ssa_nrm controller exactly. Class
# indices are 0-based here (the buffer stores class ids as r+1 elsewhere, but the
# controller stores the 0-based pos, matching the polling_* helpers).

def _draw_from_dist(w):
    """Index drawn from the (unnormalized, nonnegative) weight vector w. Returns
    0 when the total weight is not positive."""
    w = np.asarray(w, dtype=float)
    tot = float(w.sum())
    if tot <= 0.0:
        return 0
    u = np.random.random() * tot
    c = 0.0
    for i in range(w.size):
        c += float(w[i])
        if u < c:
            return i
    return int(w.size) - 1


def _cache_access(sn, ind, read_class, contents):
    """Simulate one cache READ at cache node IND by a class-READ_CLASS job over the
    cache state CONTENTS (totalCacheCapacity content slots followed, when a
    retrieval system is present, by a per-item retrieval-occupancy bitmap).
    Returns (out_class, new_contents, category): out_class is the class the job
    leaves in -- out_class == -1 means the request was absorbed as a delayed hit
    and produces nothing -- and category is 1 hit, 2 miss/retrieval-complete, 3
    delayed-hit, 4 begin-retrieval. A faithful sample-path port of
    State.afterEventCache (READ, is_simulation): non-retrieval hit/miss with all
    replacement policies, plus the retrieval (delayed-hit) system."""
    from ...state.after_event_cache import _cpos, _get_ac
    npm = sn.nodeparam[ind]
    m = np.atleast_1d(_np_get(npm, 'itemcap', [0])).astype(int)
    h = int(len(m))
    ac = _np_get(npm, 'accost', None)
    repl = _np_get(npm, 'replacestrat', None)
    repl_name = repl.name if hasattr(repl, 'name') else str(repl)
    if repl_name == 'HLRU':
        repl_name = 'LRU'
    tcc = int(_np_get(npm, 'total_cache_capacity', 0) or 0)
    if tcc <= 0:
        tcc = int(np.sum(m))
    qadm = float(_np_get(npm, 'qlru', 1.0))
    pread = _np_get(npm, 'pread', None)
    p = np.atleast_1d(np.asarray(pread[read_class], dtype=float))
    hitcl = np.atleast_1d(_np_get(npm, 'hitclass', [])).astype(int)
    misscl = np.atleast_1d(_np_get(npm, 'missclass', [])).astype(int)
    _rci_raw = _np_get(npm, 'retrieval_class_indices', None)
    rci = np.array(list(_rci_raw), dtype=int) if _rci_raw is not None and len(_rci_raw) else np.array([], dtype=int)
    rci = rci[rci >= 0]
    retr_classes = _np_get(npm, 'retrieval_classes', None)
    rsc = _np_get(npm, 'retrieval_system_capacity', 0)
    has_retr = bool(np.any(np.atleast_1d(rsc).astype(float) > 0)) if rsc is not None else False
    is_from_retrieval = int(read_class) in set(int(x) for x in rci.tolist())
    var = list(contents)        # read-only original
    varp = list(contents)       # working copy (returned)

    # default access-cost chain (accost=None): miss->list1 (reject 0), hit in list j promotes to j+1, hit in last list stays.
    if ac is None and h > 0:
        Rmat = np.zeros((h + 1, h + 1))
        for jj in range(h):
            Rmat[jj, jj + 1] = 1.0
        Rmat[h, h] = 1.0
        def _acv(row, col):
            return Rmat[row, col]
    else:
        def _acv(row, col):
            return _get_ac(ac, read_class, k, row, col)

    k = _draw_from_dist(p)      # 0-indexed item; the cached value is k+1
    miss_row = np.array([_acv(0, col) for col in range(h + 1)])
    l = _draw_from_dist(miss_row)   # 0 => reject, 1..h => list l-1
    posk = -1
    for col in range(min(tcc, len(var))):
        if var[col] == k + 1:
            posk = col
            break
    if is_from_retrieval:
        posk = -1   # a returning retrieval always COMPLETES its own miss

    total_upper = int(np.sum(m))
    if posk != -1:
        # ===================== CACHE HIT =====================
        out_class = int(hitcl[read_class])
        category = 1
        if posk < total_upper - int(m[h - 1]):
            # hit in list i < h: promote toward the last list
            repl_hit = 'LRU' if repl_name == 'QLRU' else repl_name
            mcum = np.cumsum(m)
            i_list = int(np.searchsorted(mcum, posk, side='right'))
            j_pos = posk - (int(np.sum(m[:i_list])) if i_list > 0 else 0)
            hit_row = np.array([_acv(1 + i_list, 1 + jj)
                                for jj in range(i_list, h)], dtype=float)
            inew = i_list + _draw_from_dist(hit_row)
            if repl_hit == 'FIFO':
                if inew != i_list:
                    varp[_cpos(m, i_list, j_pos)] = var[_cpos(m, inew, int(m[inew]) - 1)]
                    if m[inew] > 1:
                        hn = _cpos(m, inew, 0); tn = _cpos(m, inew, int(m[inew]) - 1)
                        varp[hn + 1:tn + 1] = var[hn:tn]
                    varp[_cpos(m, inew, 0)] = k + 1
            elif repl_hit == 'RR':
                r_pos = int(np.random.randint(int(m[inew])))
                varp[_cpos(m, i_list, j_pos)] = var[_cpos(m, inew, r_pos)]
                varp[_cpos(m, inew, r_pos)] = k + 1
            else:  # LRU, SFIFO
                if j_pos > 0:
                    hi = _cpos(m, i_list, 0)
                    varp[hi + 1:hi + j_pos + 1] = var[hi:hi + j_pos]
                varp[_cpos(m, i_list, 0)] = var[_cpos(m, inew, int(m[inew]) - 1)]
                if m[inew] > 1:
                    hn = _cpos(m, inew, 0); tn = _cpos(m, inew, int(m[inew]) - 1)
                    varp[hn + 1:tn + 1] = var[hn:tn]
                varp[_cpos(m, inew, 0)] = k + 1
        else:
            # hit in the last list h
            j_pos = posk - (int(np.sum(m[:h - 1])) if h > 1 else 0)
            if repl_name in ('LRU', 'HLRU', 'QLRU'):
                hh = _cpos(m, h - 1, 0)
                moved = var[hh + j_pos]
                if j_pos > 0:
                    varp[hh + 1:hh + j_pos + 1] = var[hh:hh + j_pos]
                varp[hh] = moved
        return out_class, varp, category

    # ===================== CACHE MISS / retrieval =====================
    if has_retr and not is_from_retrieval:
        rclass = -1
        if retr_classes is not None:
            rc_arr = np.atleast_2d(np.asarray(retr_classes))
            if k < rc_arr.shape[0] and read_class < rc_arr.shape[1]:
                rclass = int(rc_arr[k, read_class])
        if rclass != -1:
            in_retrieval = (tcc + k < len(var)) and var[tcc + k] != 0
            if in_retrieval:
                # DELAYED HIT: served by the in-flight retrieval and absorbed.
                return -1, varp, 3
            # BEGIN retrieval: switch to the item's retrieval class, mark item.
            varp[tcc + k] = 1
            return rclass, varp, 4

    # COMPLETE the miss (retrieval return, or a plain miss): clear the bit and
    # admit item k+1 per the replacement policy.
    if is_from_retrieval and (tcc + k < len(varp)):
        varp[tcc + k] = 0
    out_class = int(misscl[read_class])
    category = 2
    listidx = l - 1
    if listidx >= 0:
        head = _cpos(m, listidx, 0)
        tail = _cpos(m, listidx, int(m[listidx]) - 1)
        if repl_name in ('FIFO', 'LRU', 'SFIFO'):
            if m[listidx] > 1:
                varp[head + 1:tail + 1] = var[head:tail]
            varp[head] = k + 1
        elif repl_name == 'RR':
            r_pos = int(np.random.randint(int(m[listidx])))
            varp[_cpos(m, listidx, r_pos)] = k + 1
        elif repl_name == 'QLRU':
            if np.random.random() <= qadm:
                if m[listidx] > 1:
                    varp[head + 1:tail + 1] = var[head:tail]
                varp[head] = k + 1
    return out_class, varp, category


def _poll_serve_gate(ctrl, r):
    """1 when the polling controller ctrl is SERVING class r, else 0. This is the
    single-server gate that turns a class-r service departure on only while the
    server attends class r."""
    if ctrl is not None and len(ctrl) >= 2 and int(ctrl[0]) == 1 and int(ctrl[1]) == r:
        return 1.0
    return 0.0


def _poll_sw_rate(ctrl, pinfo):
    """Total leaving rate of the switchover phase the controller currently
    occupies, i.e. -D0(swk,swk) of the switchover PH into buffer pos; 0 unless the
    server is SWITCHING. The competition between advancing to another phase and
    absorbing is resolved at firing time by the run loop."""
    if ctrl is None or len(ctrl) < 3 or int(ctrl[0]) != 2:
        return 0.0
    pos = int(ctrl[1])
    swk = int(ctrl[2])
    D0 = pinfo['sw_d0'][pos]
    if D0 is None:
        return 0.0
    return float(-D0[swk - 1, swk - 1])


def _poll_land_ctrl(pinfo, q, mode, budget):
    """Controller row [mode, pos, swk, ctr] the server lands in after polling_next
    resolves (q, mode, budget): SERVING q with the visit budget, SWITCHING into q
    with the entry phase drawn from the switchover PH (1-based), or PARKED at the
    canonical q. Mirrors State.pollingLand specialized to the single-server
    polling station the NRM carries (exponential service, no in-service phase)."""
    if mode == 1:
        return [1, int(q), 0, int(budget)]
    if mode == 2:
        swk = _draw_from_dist(pinfo['sw_pie'][q]) + 1
        return [2, int(q), int(swk), 0]
    return [0, int(q), 0, 0]  # parked


def _sched_of(sn, ist: int):
    """Scheduling strategy at station index ist."""
    return sn.sched[ist]


def _is_buffered(ind: int, sn) -> bool:
    """True for stations whose waiting jobs are held in an ordered buffer."""
    if not sn.isstation[ind]:
        return False
    ist = int(sn.nodeToStation[ind])
    return _sched_of(sn, ist) in _BUFFERED_SCHED


def _is_list_sched(ind: int, sn) -> bool:
    """
    True for stations whose buffer holds the FULL ordered job list rather than
    only the waiting jobs.
    """
    if not sn.isstation[ind]:
        return False
    return _sched_of(sn, int(sn.nodeToStation[ind])) in _LIST_SCHED


def _is_preemptive(ind: int, sn) -> bool:
    """
    True for the preempt-resume policies, whose arrivals displace an incumbent
    instead of queueing behind it.
    """
    if not sn.isstation[ind]:
        return False
    return _sched_of(sn, int(sn.nodeToStation[ind])) in _PREEMPTIVE_SCHED


def _is_retrial_station(ind: int, sn) -> bool:
    """
    True for stations with a retrial orbit: their freed servers are not filled
    by promotion, only by a successful RETRY.

    sn.retrialType == 0 is ambiguous (it is also the "no retrial" default), so
    the test reads retrialProc, whose entry is None exactly when no retrial is
    configured.
    """
    if not sn.isstation[ind]:
        return False
    proc = getattr(sn, 'retrialProc', None)
    if proc is None:
        return False
    ist = int(sn.nodeToStation[ind])
    if ist < 0 or ist >= len(proc):
        return False
    return any(p is not None for p in proc[ist])


def _pick_preempted(nvec, buf, jnd: int, arr_class: int, R: int, smap: _Smap) -> int:
    """
    Class (1-based) of the incumbent displaced by an arrival of class arr_class
    (0-based) at node jnd, drawn in proportion to the servers' class
    occupancies, as State.afterEventStation weights its preemption branches by
    si_preempt/sum(space_srv). Returns 0 when no server holds a job.

    nvec already counts the arrival, so it is discounted here to recover the
    pre-arrival in-service composition (in-service = population minus buffer
    occupancy).
    """
    insvc = _class_counts(nvec, smap, jnd)
    for r in range(R):
        insvc[r] -= sum(1 for cl in buf if cl == r + 1)
        if r == arr_class:
            insvc[r] -= 1.0        # discount the job that just arrived
    insvc[insvc < 0] = 0.0
    tot = float(np.sum(insvc))
    if tot <= 0:
        return 0
    u = np.random.random() * tot
    acc = 0.0
    last = 0
    for r in range(R):
        if insvc[r] > 0:
            last = r + 1
            acc += insvc[r]
            if u < acc:
                return r + 1
    return last


def _pick_from_buffer(buf, sn, ist: int) -> int:
    """
    Index of the waiting job that the discipline at station ist promotes into
    service. The buffer is ordered newest-first / oldest-last, matching the
    convention of State.afterEventStation's space_buf.
    """
    sched = _sched_of(sn, ist)
    n = len(buf)
    if sched == SchedStrategy.FCFS:
        return n - 1                       # oldest
    if sched in (SchedStrategy.LCFS, SchedStrategy.LCFSPR):
        return 0                           # newest / most recently preempted
    if sched == SchedStrategy.SIRO:
        # uniform draw over waiting jobs reproduces State.afterEventStation's (nir-sir)/(ni-sum(sir)) promotion probability exactly.
        return int(np.floor(np.random.random() * n))
    arr = list(buf)
    if sched == SchedStrategy.HOL:
        # Highest priority (lowest classprio value); FCFS within the group, so
        # the oldest = the last matching position.
        best = min(_class_prio(sn, c - 1) for c in arr)
        for i in range(n - 1, -1, -1):
            if _class_prio(sn, arr[i] - 1) == best:
                return i
    if sched in (SchedStrategy.SEPT, SchedStrategy.LEPT):
        # Shortest (SEPT) or longest (LEPT) expected processing time, i.e. the
        # highest or lowest service rate. Oldest first within a class.
        rates_row = np.asarray(sn.rates[ist], dtype=float).reshape(-1)
        vals = [rates_row[c - 1] for c in arr]
        finite = [v for v in vals if np.isfinite(v)]
        if not finite:
            return n - 1
        best = max(finite) if sched == SchedStrategy.SEPT else min(finite)
        for i in range(n - 1, -1, -1):
            if np.isfinite(vals[i]) and vals[i] == best:
                return i
    raise ValueError(f"_pick_from_buffer: unsupported buffered policy {sched}")


def _has_service(sn, ist: int, k: int) -> bool:
    """Whether class k has an enabled service process at station ist."""
    return not np.isnan(sn.rates[ist, k])


def _class_prio(sn, r: int) -> int:
    if sn.classprio is None:
        return 0
    cp = np.asarray(sn.classprio).flatten()
    return int(cp[r]) if r < len(cp) else 0


# =============================================================================
# PS-family sharing factors
#
# Each returns the multiplier applied to the class-r service rate, i.e. the
# fraction of total service capacity class r receives in population state
# NPOP. All mirror the corresponding case of State.afterEventStation
# specialized to exponential (single-phase) service, where the phase
# population kir equals the class population nir.
# =============================================================================
_WEIGHTED_SCHED = (SchedStrategy.DPS, SchedStrategy.GPS,
                   SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO)


def _sched_weights(sn, ist: int, R: int) -> np.ndarray:
    """Normalized per-class DPS/GPS weights of station ist, from sn.schedparam."""
    w = np.asarray(sn.schedparam, dtype=float)[ist, :R].copy()
    tot = float(np.sum(w))
    if tot <= 0:
        raise ValueError(
            f"Station {ist} has weighted (DPS/GPS) scheduling with non-positive total weight.")
    return w / tot


def _build_sched_weights(sn, R: int):
    """
    Normalized weights per station, None for policies that do not read them.
    State.afterEventStation rejects multi-server DPS/GPS, so fail on those
    rather than silently simulating a different station.
    """
    wnorm = [None] * sn.nstations
    for ist in range(sn.nstations):
        s = _sched_of(sn, ist)
        if s in _WEIGHTED_SCHED:
            wnorm[ist] = _sched_weights(sn, ist, R)
            if sn.nservers[ist] > 1:
                raise ValueError(f"Multi-server {s} stations are not supported yet.")
    return wnorm


def _dpsshare(w, npop, r: int) -> float:
    """DPS: class r receives w_r*n_r/(w.n) of the single server."""
    den = float(np.dot(w, npop))
    if den <= 0:
        return 0.0
    return float(w[r] * npop[r] / den)


def _gpsshare(w, npop, r: int) -> float:
    """GPS: w_r/(w.c) with c_s = 1{n_s>0}: weights split across active classes."""
    if npop[r] <= 0:
        return 0.0
    den = float(np.dot(w, (np.asarray(npop) > 0).astype(float)))
    if den <= 0:
        return 0.0
    return float(w[r] / den)


def _is_urgent(sn, npop, r: int) -> bool:
    """True when class r is in the most urgent non-empty priority group.
    LINE orders priorities with lower value = more urgent."""
    occupied = [s for s in range(len(npop)) if npop[s] > 0]
    if not occupied:
        return False
    return _class_prio(sn, r) == min(_class_prio(sn, s) for s in occupied)


def _prio_group(sn, npop, r: int) -> np.ndarray:
    """Population vector restricted to the priority group of class r."""
    act = np.zeros(len(npop), dtype=float)
    pr = _class_prio(sn, r)
    for s in range(len(npop)):
        if _class_prio(sn, s) == pr:
            act[s] = npop[s]
    return act


def _prio_pop(sn, npop, r: int, c: float) -> float:
    """Population the lld factor is evaluated at: the full station population
    below capacity, the priority-group population above it."""
    ni = float(np.sum(npop))
    if ni <= c or not _is_urgent(sn, npop, r):
        return ni
    return float(np.sum(_prio_group(sn, npop, r)))


def _prio_vec(sn, npop, r: int, c: float) -> np.ndarray:
    """Population vector the cd factor is evaluated at for DPSPRIO/GPSPRIO.
    PSPRIO instead uses the full vector in both branches; that asymmetry is
    inherited from State.afterEventStation and reproduced here."""
    ni = float(np.sum(npop))
    if ni <= c or not _is_urgent(sn, npop, r):
        return np.asarray(npop, dtype=float)
    return _prio_group(sn, npop, r)


def _psprioshare(sn, npop, r: int, c: float) -> float:
    """PSPRIO: PS below capacity; above it only the most urgent non-empty group
    shares the servers and every other class is frozen."""
    ni = float(np.sum(npop))
    if ni <= 0:
        return 0.0
    if ni <= c:
        return float((npop[r] / ni) * min(ni, c))
    if not _is_urgent(sn, npop, r):
        return 0.0
    niprio = float(np.sum(_prio_group(sn, npop, r)))
    return float((npop[r] / niprio) * min(niprio, c))


def _dpsprioshare(sn, w, npop, r: int, c: float) -> float:
    """DPSPRIO: DPS below capacity, DPS restricted to the urgent group above."""
    ni = float(np.sum(npop))
    if ni <= 0:
        return 0.0
    if ni <= c:
        return _dpsshare(w, npop, r)
    if not _is_urgent(sn, npop, r):
        return 0.0
    return _dpsshare(w, _prio_group(sn, npop, r), r)


def _gpsprioshare(sn, w, npop, r: int, c: float) -> float:
    """GPSPRIO: GPS below capacity, GPS restricted to the urgent group above."""
    ni = float(np.sum(npop))
    if ni <= 0:
        return 0.0
    if ni <= c:
        return _gpsshare(w, npop, r)
    if not _is_urgent(sn, npop, r):
        return 0.0
    return _gpsshare(w, _prio_group(sn, npop, r), r)


def _is_ps_family(sched) -> bool:
    """The processor-sharing family, whose utilization is the share of service
    capacity a class receives rather than a server-occupancy count."""
    return sched in (SchedStrategy.PS, SchedStrategy.LPS, SchedStrategy.DPS,
                     SchedStrategy.GPS, SchedStrategy.PSPRIO,
                     SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO)


def _ps_family_util(sn, w, npop, k: int, servers: float, sched) -> float:
    """Class-k utilization of a PS-family station: the share of service capacity
    class k receives, divided by the server count. The lld/cd scalings are
    excluded on purpose, as they rescale work rather than server occupancy."""
    if sched in (SchedStrategy.PS, SchedStrategy.LPS):
        total = float(np.sum(npop))
        if total <= 0:
            return 0.0
        return float((npop[k] / total) * min(servers, total) / servers)
    if sched == SchedStrategy.DPS:
        return _dpsshare(w, npop, k) / servers
    if sched == SchedStrategy.GPS:
        return _gpsshare(w, npop, k) / servers
    if sched == SchedStrategy.PSPRIO:
        return _psprioshare(sn, npop, k, servers) / servers
    if sched == SchedStrategy.DPSPRIO:
        return _dpsprioshare(sn, w, npop, k, servers) / servers
    if sched == SchedStrategy.GPSPRIO:
        return _gpsprioshare(sn, w, npop, k, servers) / servers
    return 0.0


# =============================================================================
# G-network signals
#
# A signal class never joins the station it reaches: it acts on the jobs already
# there and is annihilated (api.state.signal_removal / MATLAB
# State.afterEventStationSignal). That makes it an arrival-side effect exactly
# like balking, so the departure rate is unchanged and only the arrival outcome
# differs. The reference enumerates every victim subset with its probability
# because it builds a generator; a simulator instead draws the batch size and
# the victims, which is equivalent and avoids the enumeration.
# =============================================================================
class _Sig:
    __slots__ = ('on', 'sn', 'issignal', 'nonsignal')

    def __init__(self, sn):
        self.on = False
        iss = getattr(sn, 'issignal', None)
        if iss is None or np.size(iss) == 0 or not np.any(np.asarray(iss)):
            return
        self.on = True
        self.sn = sn
        self.issignal = np.asarray(iss).ravel().astype(bool)
        self.nonsignal = [r for r in range(len(self.issignal)) if not self.issignal[r]]


def _sig_is_signal_arrival(sig, dest_pos, R, smap: _Smap) -> bool:
    """True when the state slot dest_pos is a signal class at a station."""
    s = int(smap.cls[dest_pos])
    jnd = int(smap.node[dest_pos])
    return bool(sig.issignal[s]) and bool(sig.sn.isstation[jnd])


def _sig_remove_job(nvec, smap: _Smap, jnd: int, r: int) -> None:
    """
    Remove one class-r job from node jnd, drawing which phase it was in from
    the per-phase populations. Jobs of a class in the same phase are
    exchangeable, so the draw is proportional to the phase occupancies; for a
    single-phase class this is the plain decrement of the class's only slot.
    """
    sl = smap.class_slots(jnd, r)
    counts = np.asarray(nvec[sl], dtype=float)
    tot = float(np.sum(counts))
    if tot <= 0:
        return
    u = np.random.random() * tot
    acc = 0.0
    for kk in range(counts.size):
        acc += counts[kk]
        if counts[kk] > 0 and u < acc:
            nvec[sl.start + kk] -= 1.0
            return
    nvec[sl.start + int(np.flatnonzero(counts > 0)[-1])] -= 1.0


def _sig_apply(sig, nvec, buffers, dest_pos, R, mi, smap: _Smap):
    """Apply the arrival of a signal class at a station: pick victims, remove."""
    from ...state.signal_removal import is_catastrophe_signal, signal_batch_pmf
    from ....constants import RemovalPolicy
    sn = sig.sn
    jnd = int(smap.node[dest_pos])
    cls = int(smap.cls[dest_pos])
    ist = int(sn.nodeToStation[jnd])

    # CATASTROPHE empties the station of every job, ignoring the batch-size
    # distribution: a catastrophe removes all jobs by definition.
    if is_catastrophe_signal(sn, cls):
        nvec[smap.node_slots(jnd)] = 0.0
        buffers[jnd].clear()
        return

    # eligible victim classes: forJobClass restricts to that class, else every non-signal class (Gelenbe negative customer); sn.signaltarget is 0-based (-1=none).
    tgt = -1
    st = getattr(sn, 'signaltarget', None)
    if st is not None and cls < len(st):
        tgt = int(st[cls])
    tgtclasses = [tgt] if tgt >= 0 else list(sig.nonsignal)
    tgtclasses = [r for r in tgtclasses if _class_pop(nvec, smap, jnd, r) > 0]
    ntot = float(sum(_class_pop(nvec, smap, jnd, r) for r in tgtclasses))
    if not tgtclasses or ntot <= 0:
        return   # no victim: the signal simply vanishes

    # batch size drawn from the reference pmf, clipped at the eligible population.
    kvals, kprobs = signal_batch_pmf(sn, cls, int(ntot))
    u = np.random.random()
    acc = 0.0
    k = int(kvals[-1])
    for kv, kp in zip(kvals, kprobs):
        acc += kp
        if acc >= u:
            k = int(kv)
            break

    policy = RemovalPolicy.RANDOM
    pol = getattr(sn, 'signalrempolicy', None)
    if pol is not None and cls < len(pol) and pol[cls] is not None:
        policy = pol[cls]

    for _step in range(k):
        if not _sig_remove_one(sn, nvec, buffers, jnd, ist, tgtclasses,
                               policy, R, mi, smap):
            break   # already drained


def _sig_remove_one(sn, nvec, buffers, jnd, ist, tgtclasses, policy, R, mi,
                    smap: _Smap) -> bool:
    """Remove one victim under the signal's removal policy."""
    from ....constants import RemovalPolicy
    buf = buffers[jnd]
    wait_idx = [i for i in range(len(buf)) if (buf[i] - 1) in tgtclasses]
    nwait = len(wait_idx)
    nsrv = 0.0
    for r in tgtclasses:
        nsrv += max(0.0, _class_pop(nvec, smap, jnd, r)
                    - sum(1 for cl in buf if cl == r + 1))
    if nwait == 0 and nsrv <= 0:
        return False

    # FCFS/LCFS rank the ordered (newest-first) buffer by age; count-only buffers (SIRO/SEPT/LEPT) have no age and degenerate to a uniform draw.
    is_ordered = _sched_of(sn, ist) in (SchedStrategy.FCFS, SchedStrategy.HOL,
                                        SchedStrategy.LCFS)
    age_ordered = is_ordered and policy in (RemovalPolicy.FCFS, RemovalPolicy.LCFS)
    if age_ordered and nwait > 0:
        pick = wait_idx[-1] if policy == RemovalPolicy.FCFS else wait_idx[0]
        victim = buf[pick] - 1
        del buffers[jnd][pick]
        _sig_remove_job(nvec, smap, jnd, victim)
        return True

    # RANDOM draws uniformly over waiting and in-service alike; FCFS/LCFS drain
    # the waiting line before reaching into the servers.
    if policy == RemovalPolicy.RANDOM:
        total = nwait + nsrv
    else:
        total = nwait if nwait > 0 else nsrv
    if nwait > 0 and (policy != RemovalPolicy.RANDOM
                      or np.random.random() * total < nwait):
        pick = wait_idx[int(np.random.random() * nwait)]
        victim = buf[pick] - 1
        del buffers[jnd][pick]
        _sig_remove_job(nvec, smap, jnd, victim)
        return True

    # an in-service victim, uniform over the eligible in-service jobs
    target = np.random.random() * nsrv
    acc = 0.0
    for r in tgtclasses:
        cnt = max(0.0, _class_pop(nvec, smap, jnd, r)
                  - sum(1 for cl in buf if cl == r + 1))
        acc += cnt
        if cnt > 0 and target < acc:
            _sig_remove_job(nvec, smap, jnd, r)
            # freed server pulls the head of line, i.e. the leaving buffer job (in-service = population minus buffer occupancy).
            total_new = float(np.sum(nvec[smap.node_slots(jnd)]))
            if len(buffers[jnd]) > max(0.0, total_new - mi[jnd]):
                buffers[jnd].pop()   # head of line: the oldest waiting job
            return True
    return False


# =============================================================================
# Balking and reneging
# =============================================================================
class _Balk:
    """
    Per (station,class) balking threshold table, for the QUEUE_LENGTH strategy.

    A balked job is lost: it has left its source but never joins the
    destination, so the departure RATE is unchanged and only the arrival
    outcome differs (State.afterEventStation scales the admitted branches by
    1-balkProb and adds a balked branch of probability balkProb that leaves the
    destination state untouched). Only QUEUE_LENGTH is a pure function of the
    state vector; EXPECTED_WAIT / COMBINED depend on the mean wait and are
    rejected by the analyzers.
    """

    __slots__ = ('on', 'strategy', 'thresholds', 'isstation', 'nodeToStation')

    def __init__(self, sn):
        self.on = False
        bs = getattr(sn, 'balkingStrategy', None)
        if bs is None or np.size(bs) == 0:
            return
        from ....lang.base import BalkingStrategy
        if not np.any(np.asarray(bs) == int(BalkingStrategy.QUEUE_LENGTH)):
            return
        self.on = True
        self.strategy = np.asarray(bs)
        self.thresholds = sn.balkingThresholds
        self.isstation = sn.isstation
        self.nodeToStation = sn.nodeToStation


def _balk_draw(balk, nvec, dest_pos, R, smap: _Smap) -> bool:
    """
    True if the job routed to state slot dest_pos balks. The threshold table is
    scanned in order and the FIRST interval containing the pre-arrival total
    station population wins, matching State.afterEventStation.
    """
    from ....lang.base import BalkingStrategy
    jnd = int(smap.node[dest_pos])
    s = int(smap.cls[dest_pos])
    if not balk.isstation[jnd]:
        return False
    ist = int(balk.nodeToStation[jnd])
    if ist < 0 or balk.strategy[ist, s] != int(BalkingStrategy.QUEUE_LENGTH):
        return False
    qlen = float(np.sum(nvec[smap.node_slots(jnd)]))   # pre-arrival population
    th = balk.thresholds[ist][s] if balk.thresholds is not None else None
    if not th:
        return False
    balk_prob = 0.0
    for t in th:
        if qlen >= t[0] and qlen <= t[1]:
            balk_prob = float(t[2])
            break
    return balk_prob > 0.0 and np.random.random() < balk_prob


def _is_physical_capacity(sn, ist: int, job_class: int) -> bool:
    """True when station ist declares a PHYSICAL drop rule (anything other than
    WaitingQueue) for job_class. Mirrors State.isPhysicalCapacity: only a
    physical-drop station loses a refused arrival; a WaitingQueue region parks
    it instead (handled by the FCR gate, not here)."""
    droprule = getattr(sn, 'droprule', None)
    if droprule is None:
        return False
    dr = np.asarray(droprule)
    if ist >= dr.shape[0] or (dr.ndim > 1 and job_class >= dr.shape[1]):
        return False
    from ...sn.network_struct import DropStrategy as _Drop
    val = int(dr[ist, job_class]) if dr.ndim > 1 else int(dr[ist])
    return val != int(_Drop.WAITQ)


def _capacity_loss(sn, nvec, dest_pos, R, smap: _Smap) -> bool:
    """True if an OPEN-class job routed to state slot dest_pos is lost because
    its destination is a physically finite-capacity station already at cap.
    Mirrors State.afterEventStation's hasRoom + arrivalIsLost gate, which the
    NRM reaction network otherwise omits -- letting a capped queue (e.g. an
    open PAS/OI M/M/K/N station) overflow far past sn.cap under simulation.

    A refused CLOSED job must BLOCK, not vanish from the conserved population,
    so _arrival_is_lost returns False for it and this gate does not drop it.
    Mirrors MATLAB/Java capacityLoss in solver_ssa_nrm."""
    from ...state.after_event_station import _arrival_is_lost
    jnd = int(smap.node[dest_pos])
    dst_c = int(smap.cls[dest_pos])
    nodeToStation = np.asarray(sn.nodeToStation).flatten()
    if jnd >= nodeToStation.size:
        return False
    ist = int(nodeToStation[jnd])
    if ist < 0:
        return False
    if not _is_physical_capacity(sn, ist, dst_c) or not _arrival_is_lost(sn, ist, dst_c):
        return False
    cap = np.asarray(sn.cap).flatten() if getattr(sn, 'cap', None) is not None else None
    total = sum(_class_pop(nvec, smap, jnd, r) for r in range(R))   # pre-arrival
    if cap is not None and ist < cap.size and np.isfinite(cap[ist]) and total >= cap[ist]:
        return True
    classcap = getattr(sn, 'classcap', None)
    if classcap is not None:
        cc = np.asarray(classcap)
        if ist < cc.shape[0] and dst_c < cc.shape[1]:
            lim = cc[ist, dst_c]
            if lim > 0 and _class_pop(nvec, smap, jnd, dst_c) >= lim:
                return True
    return False


# =============================================================================
# Finite capacity regions (DROP rule)
#
# A region constrains an aggregate of the per-class populations of its member
# stations, which is a linear function of the NRM state vector, so admission is
# a multiplicative 0/1 gate on the routing draw. The DROP rule censors the
# refused transition, and censoring an exponential transition is exactly what
# zeroing its share of the propensity does. WAITQ instead parks refused jobs in
# a per-region FIFO, which is extra state the reaction network does not carry,
# so those models are routed to the serial engine and never reach here.
# =============================================================================
class _Fcr:
    """Per-region member nodes and admission caps; mirrors the FCR precompute
    of the serial engine (solver_ssa) and of MATLAB solver_ssa_nrm."""

    __slots__ = ('on', 'member_node', 'class_cap', 'global_cap', 'mem_cap',
                 'sz', 'A', 'b', 'F', 'waitq', 'any_waitq')

    def __init__(self, sn):
        self.on = False
        self.any_waitq = False
        F = int(getattr(sn, 'nregions', 0) or 0)
        if F == 0:
            return
        self.on = True
        self.F = F
        K = sn.nclasses
        M = sn.nstations
        # per (region,class) DROP/WAITQ rule compared BY NAME; DROP destroys, WAITQ parks in the region FIFO.
        self.waitq = np.zeros((F, K), dtype=bool)
        rule = getattr(sn, 'regionrule', None)
        if rule is not None and np.size(rule) > 0:
            from ....lang.base import DropStrategy as _BaseDrop
            drop_id = float(int(_BaseDrop.DROP))
            arr = np.asarray(rule, dtype=float)
            for f in range(F):
                for r in range(K):
                    if arr[f, r] != drop_id:
                        self.waitq[f, r] = True
                        self.any_waitq = True
        self.member_node = np.zeros((F, sn.nnodes), dtype=bool)
        self.class_cap = []
        self.global_cap = np.full(F, np.inf)
        self.mem_cap = np.full(F, np.inf)
        self.sz = []
        self.A = [None] * F
        self.b = [None] * F
        for f in range(F):
            Rmat = np.asarray(sn.region[f], dtype=float)      # M x (K+1)
            # region membership read from sn.regionmembers, not derived from region[f]; see _kb/04-networkstruct.md Python native network.py section.
            memvec = -np.ones(M)
            rmm = getattr(sn, 'regionmaxmem', None)
            if rmm is not None and len(rmm) > f and rmm[f] is not None and len(np.asarray(rmm[f]).ravel()):
                memvec = np.asarray(rmm[f], dtype=float).ravel()
            mask = sn_region_members(sn, f, Rmat, memvec)
            members = np.flatnonzero(mask)
            ccap = np.full(K, np.inf)
            for r in range(K):
                cv = Rmat[members, r]
                cv = cv[cv != -1]
                if cv.size:
                    ccap[r] = cv.min()
            self.class_cap.append(ccap)
            gv = Rmat[members, K]
            gv = gv[gv != -1]
            if gv.size:
                self.global_cap[f] = gv.min()
            if rmm is not None and len(rmm) > f and rmm[f] is not None and len(np.asarray(rmm[f]).ravel()):
                mv = np.asarray(rmm[f], dtype=float).ravel()[members]
                mv = mv[mv != -1]
                if mv.size:
                    self.mem_cap[f] = mv.min()
            self.sz.append(np.asarray(sn.regionsz, dtype=float)[f, :])
            rlc = getattr(sn, 'regionlincon', None)
            if rlc is not None and len(rlc) > f and rlc[f] is not None:
                Af, bf = rlc[f][0], rlc[f][1]
                if Af is not None and np.asarray(Af).size:
                    self.A[f] = np.asarray(Af, dtype=float)
                    self.b[f] = np.asarray(bf, dtype=float).ravel()
            for ist in members:
                self.member_node[f, int(sn.stationToNode[ist])] = True


def _fcr_violates(xn, ccap, gcap, memcap, sz, A, b) -> bool:
    """True if per-class population vector xn breaks any admission constraint."""
    if np.any(xn > ccap) or float(np.sum(xn)) > gcap or float(np.dot(xn, sz)) > memcap:
        return True
    if A is not None:
        return bool(np.any(A.dot(xn) > b))
    return False


def _fcr_region_pop(nvec, member_row, R, smap: _Smap):
    """Per-class population of a region, read off the NRM state vector."""
    x = np.zeros(R)
    for jnd in np.flatnonzero(member_row):
        x += _class_counts(nvec, smap, int(jnd))
    return x


def _fcr_admits(fcr, nvec, src_node, src_class, dst_node, dst_class, R,
                smap: _Smap) -> bool:
    """
    True if a class-dst_class job may enter dst_node, having just left src_node
    as src_class. Only regions containing the destination can refuse the move; a
    move whose source is in the same region frees a slot first, so the departure
    is accounted for before the arrival is tested.
    """
    return _fcr_refusing_region(fcr, nvec, src_node, src_class, dst_node,
                                dst_class, R, smap) < 0


def _fcr_refusing_region(fcr, nvec, src_node, src_class, dst_node, dst_class, R,
                         smap: _Smap) -> int:
    """
    Index of the FIRST region that refuses a class-dst_class job entering
    dst_node, having just left src_node as src_class; -1 if every region admits
    it. Same admission test as :func:`_fcr_admits`, but it names the refusing
    region so the caller can consult that region's DROP/WAITQ rule. Mirrors
    MATLAB fcrRefusingRegion. A src_node < 0 means the mover has no live source
    in the state (a WAITQ release, whose job already left its source when it was
    parked), so no source slot is freed.
    """
    if not fcr.on:
        return -1
    for f in range(fcr.F):
        if not fcr.member_node[f, dst_node]:
            continue   # this region does not constrain the destination
        x = _fcr_region_pop(nvec, fcr.member_node[f], R, smap)
        if src_node >= 0 and fcr.member_node[f, src_node]:
            x[src_class] -= 1.0
        x[dst_class] += 1.0
        if _fcr_violates(x, fcr.class_cap[f], fcr.global_cap[f], fcr.mem_cap[f],
                         fcr.sz[f], fcr.A[f], fcr.b[f]):
            return f
    return -1


def _draw_from_dist(p) -> int:
    """Index drawn from the (unnormalized, nonnegative) weight vector p."""
    tot = float(np.sum(p))
    if tot <= 0.0:
        return 0
    u = np.random.random()
    c = 0.0
    for i in range(len(p)):
        c += p[i] / tot
        if c > u:
            return i
    return len(p) - 1


def _fcr_release_cascade(fcr, nvec, buffers, fcr_buf, mi, R, sn, smap: _Smap,
                         svcph=None, bufph_node=None):
    """
    Strict-FIFO head-of-line release of parked WAITQ tokens: admit each region's
    FIFO head while the admission constraints permit, applying the arrival to
    the destination station (entry-phase slot plus buffer join). Mirrors MATLAB
    fcrReleaseCascade. A token encodes (dst_node, dst_class) as dst_node*R +
    dst_class; the entry phase is drawn at release, as a routed arrival draws it.
    Loops until a full pass frees nothing, so a release that frees capacity
    elsewhere cascades. Returns (number admitted, svc_changed).
    """
    released = 0
    svc_changed = False
    progress = True
    while progress:
        progress = False
        for f in range(len(fcr_buf)):
            if not fcr_buf[f]:
                continue
            tok = fcr_buf[f][0]
            dst_node = tok // R
            dst_class = tok % R
            # The parked job already left its source, so admission is tested with
            # the source term absent (src_node = -1).
            if _fcr_refusing_region(fcr, nvec, -1, dst_class, dst_node, dst_class, R, smap) >= 0:
                continue   # head-of-line: this FIFO stays blocked
            if bufph_node is not None and bufph_node[dst_node]:
                # buffered-PH destination lands in the class total slot; entry into service/phase is decided in _apply_arrival_buffer.
                nvec[int(smap.ph_off[dst_node, dst_class])] += 1.0
            else:
                pentry = _entry_probs(sn, dst_node, dst_class, int(smap.nph[dst_node, dst_class]))
                ke = _draw_from_dist(pentry)
                dslot = int(smap.ph_off[dst_node, dst_class]) + ke
                nvec[dslot] += 1.0
            if _apply_arrival_buffer(dst_node, dst_class, nvec, buffers, mi, R, sn, smap,
                                     svcph, bufph_node):
                svc_changed = True
            del fcr_buf[f][0]
            released += 1
            progress = True
    return released, svc_changed


def _impatience_nrm_ok(sn) -> bool:
    """
    True when the model's balking and reneging are both forms the NRM can
    evaluate off the state vector.

    QUEUE_LENGTH balking is a pure function of the population, so the NRM draws
    it at firing time; EXPECTED_WAIT and COMBINED depend on the mean waiting
    time. Reneging abandons at the aggregate rate (waiting count)*mu, which is
    only correct when patience is memoryless; phase-type patience would need
    each waiting job's remaining phase. The serial engine rejects the same two
    combinations outright.
    """
    from ....lang.base import BalkingStrategy, ImpatienceType
    from ....constants import ProcessType
    bs = getattr(sn, 'balkingStrategy', None)
    if bs is not None and np.size(bs) > 0:
        bad = np.asarray(bs)
        if np.any((bad != 0) & (bad != int(BalkingStrategy.QUEUE_LENGTH))):
            return False
    ic = getattr(sn, 'impatienceClass', None)
    it = getattr(sn, 'impatienceType', None)
    if ic is not None and it is not None and np.size(ic) > 0:
        ic = np.asarray(ic)
        it = np.asarray(it)
        # impatienceType==0 is ambiguous and cannot be tested as `!= EXP` in python; see _kb/07-cross-language-parity.md impatienceType note.
        exp_id = int(ProcessType.EXP.value if hasattr(ProcessType.EXP, 'value')
                     else ProcessType.EXP)
        bad = (ic == int(ImpatienceType.RENEGING)) & (it != 0) & (it != exp_id)
        if np.any(bad):
            return False
    return True


# NRM phase-expansion-exact stations/EXT exclusion; see _kb/06-solver-catalog.md SSA Phase-type service section.
_PHASE_EXACT_SCHED = (
    SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.LPS,
    SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.PSPRIO,
    SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO,
    SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO,
    SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT,
)

# The non-preemptive buffered family, whose phase-type service is expanded via
# svcph. LCFSPR is intentionally excluded (see _PHASE_EXACT_SCHED).
_BUFPH_SCHED = (SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO,
                SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT)


def _phase_nrm_ok(sn) -> bool:
    """
    True when every non-exponential service sits at a station whose rate law the
    NRM expands exactly.

    Phase expansion splits the class-level share across a class's phases in the
    ratio kir/nir, which needs only the per-phase populations -- true of the
    INF/PS family, where every job present is in service. The non-preemptive
    buffered family (FCFS/LCFS/SIRO/HOL/SEPT/LEPT) is also expanded, via the
    auxiliary in-service phase multiset svcph (see _BUFPH_SCHED). Preemptive
    LCFSPR (would need the preempted job's phase remembered) and POLLING (a
    controller carries a single in-service job) are excluded and still need the
    serial engine.
    """
    from ....constants import ProcessType
    pid = getattr(sn, 'procid', None)
    if pid is None or np.size(pid) == 0:
        return True
    pid = np.atleast_2d(np.asarray(pid, dtype=object))
    # ProcessType is a plain Enum here with EXP=1 (not MATLAB's 0); compare via .value, never a bare int cast.
    exp_id = int(ProcessType.EXP.value)
    dis_id = int(ProcessType.DISABLED.value)
    imm_id = int(ProcessType.IMMEDIATE.value)
    for ist in range(min(pid.shape[0], int(sn.nstations))):
        for r in range(min(pid.shape[1], int(sn.nclasses))):
            v = pid[ist, r]
            if v is None or (isinstance(v, float) and np.isnan(v)):
                continue
            vid = int(v.value) if hasattr(v, 'value') else int(v)
            # An immediate class carries no (D0,D1) and is simulated as a fast
            # exponential (see _refresh_phase_fields), so it is single-phase.
            if vid in (exp_id, dis_id, imm_id):
                continue
            if _sched_of(sn, ist) not in _PHASE_EXACT_SCHED:
                return False
    return True


def _cache_nrm_ok(sn) -> bool:
    """True for every Cache node. The NRM models a cache access as an immediate
    state-dependent class switch (read -> hit/miss/retrieval) at the cache node,
    applying the same replacement logic as State.afterEventCache to the cache
    contents carried alongside the buffers, INCLUDING the retrieval (delayed-hit)
    system: a miss for an item not yet being fetched begins a retrieval (the job
    is routed to the fetch queue and returns to complete the miss), and a
    concurrent request for an item already being fetched is absorbed as a delayed
    hit -- matching the serial engine's sample-path semantics."""
    return True


def _routing_nrm_ok(sn) -> bool:
    """
    True when every routing strategy in the model is one the NRM resolves at
    firing time off the state vector. JSQ and SQ both select from the
    candidate queue lengths, so both are native. RROBIN/WRROBIN need a rotation
    pointer, which no rate reads, so the NRM carries it as auxiliary state
    alongside the buffers rather than as a reaction-network dimension (see _Rr).
    RL needs an external policy and still requires the serial engine.
    """
    from ....constants import RoutingStrategy
    if getattr(sn, 'routing', None) is None:
        return True
    statedep = {int(RoutingStrategy.RL.value)}
    kch = int(RoutingStrategy.SQ.value)
    arr = np.asarray(sn.routing, dtype=object)
    for ind in range(arr.shape[0]):
        for r in range(arr.shape[1]):
            v = arr[ind, r]
            if v is None:
                continue
            v = int(v.value) if hasattr(v, 'value') else int(v)
            if v in statedep:
                return False
    return True


def _fcr_nrm_ok(sn) -> bool:
    """
    True: finite capacity regions are supported under both rules. DROP is
    reproduced by censoring the refused transition; WAITQ parks the refused job
    in a per-region FIFO carried explicitly by the NRM engine
    (:func:`_fcr_release_cascade`) and admits it head-of-line as capacity frees.
    No region rule forces a serial fallback. Kept as a predicate so the dispatch
    structure mirrors MATLAB solver_ssa_analyzer.fcrNrmOK.
    """
    return True


def _oi_pos_rates(mu_fun, G, c):
    """
    Per-position (Delta_mu, ejected class) of the ordered list c, skipping the
    positions that receive no service. Shared by the rate law and the departure
    rewrite so both split the completion the same way.
    """
    from ...state.after_event_station import passAndSwap
    out = []
    c0 = [int(x) - 1 for x in c]           # 0-based classes for svcRateFun
    mu_prev = 0.0                          # mu of the empty prefix is 0
    for p in range(len(c)):
        mu_cur = float(mu_fun(c0[:p + 1]))
        ratep = mu_cur - mu_prev
        mu_prev = mu_cur
        if ratep <= 0:
            continue                       # position p receives no service
        _cnew, dep_class, _chain = passAndSwap(list(c), p, G)
        out.append((p, ratep, int(dep_class)))
    return out


def _oi_param(sn, ind):
    """(svcRateFun, swapGraph) of an order-independent / pass-and-swap node."""
    npj = sn.nodeparam[ind] if (getattr(sn, 'nodeparam', None) is not None
                                and ind in sn.nodeparam) else None
    mu_fun = npj.get('svcRateFun') if isinstance(npj, dict) else None
    G = npj.get('swapGraph') if isinstance(npj, dict) else None
    if mu_fun is None:
        raise ValueError('PAS/OI station has no service rate function mu(c); '
                         'set it via set_service(lambda c: ...).')
    return mu_fun, G


def _pas_in_svc(sn, ind: int, c, r: int) -> int:
    """
    Number of class-r (0-based) jobs in service at a PAS/OI station holding the
    ordered list c: the positions whose marginal rate increment Delta_mu is
    positive. Mirrors the sir that the PAS branch of State.toMarginal reports,
    which is what solver_ctmc_analyzer divides by the server count.

    The class counted is the one AT the position, not the one pass-and-swap
    would eject, so no swap graph enters here. Counting positions rather than
    rate shares is what makes a job served by several server types count once
    instead of 1/rate.
    """
    if len(c) == 0:
        return 0
    mu_fun, _G = _oi_param(sn, ind)
    c0 = [int(x) - 1 for x in c]           # 0-based classes for svcRateFun
    n = 0
    mu_prev = 0.0
    for p in range(len(c)):
        mu_cur = float(mu_fun(c0[:p + 1]))
        if mu_cur - mu_prev > 0 and c0[p] == r:
            n += 1
        mu_prev = mu_cur
    return n


def _oirate(mu_fun, G, c, r: int) -> float:
    """
    Aggregate class-r (0-based) departure rate of an order-independent station
    holding the ordered list c. Mirrors the DEP branch of
    after_event_station_pas: every position contributes its own service token at
    Delta_mu, and pass-and-swap decides which class actually leaves.
    """
    if len(c) == 0:
        return 0.0
    return float(sum(rp for (_p, rp, dep) in _oi_pos_rates(mu_fun, G, c)
                     if dep == r + 1))


def _oi_depart(sn, ind: int, buf, r: int):
    """
    Apply the pass-and-swap rewrite for a class-r (0-based) departure at OI
    station ind. The completing position is drawn among those whose
    pass-and-swap ejects class r, weighted by that position's own Delta_mu --
    the same split after_event_station_pas enumerates.
    """
    from ...state.after_event_station import passAndSwap
    mu_fun, G = _oi_param(sn, ind)
    c = list(buf)
    cand = [(p, rp) for (p, rp, dep) in _oi_pos_rates(mu_fun, G, c) if dep == r + 1]
    if not cand:
        return deque(c)   # this class cannot depart from the current list
    u = np.random.random() * sum(rp for (_p, rp) in cand)
    acc = 0.0
    pick = cand[-1][0]
    for (p, rp) in cand:
        acc += rp
        if u < acc:
            pick = p
            break
    cnew, _dep, _chain = passAndSwap(c, pick, G)
    return deque(int(x) for x in cnew)


# =============================================================================
# Round-robin routing pointers
#
# The pointer that RROBIN/WRROBIN walk is a per-(node,class) local variable, not
# a population, and no rate depends on it: it only decides where a departure
# goes. In a generator that makes it a genuine extra state dimension, but a
# simulator carries it as auxiliary state alongside the buffers, which is what
# happens here -- so it must NOT enter the reaction network.
# after_event_router advances the pointer on the departure and the routing
# closure then reads state_AFTER, so the destination used is the one the pointer
# lands on: advance first, then select.
# =============================================================================
class _Rr:
    """
    Per-(node,class) round-robin cycles and pointers, seeded from the initial
    state. RROBIN's state slot holds a destination NODE INDEX, WRROBIN's a
    POSITION in the weighted cycle (destinations repeated by weight), matching
    _update_rrobin_pointer. Positions are 0-based here.
    """

    def __init__(self, sn, state):
        from ....constants import RoutingStrategy
        self.on = False
        self.isrr = {}
        self.cycle = {}
        self.pos = {}
        if getattr(sn, 'routing', None) is None:
            return
        rr_val = int(RoutingStrategy.RROBIN.value)
        wrr_val = int(RoutingStrategy.WRROBIN.value)
        R = sn.nclasses
        for ind in range(sn.nnodes):
            for r in range(R):
                v = sn.routing[ind, r]
                if v is None:
                    continue
                v = int(v.value) if hasattr(v, 'value') else int(v)
                if v not in (rr_val, wrr_val):
                    continue
                is_wrr = (v == wrr_val)
                cyc = self._cycle_of(sn, ind, r, is_wrr)
                if cyc is None or len(cyc) == 0:
                    continue
                self.on = True
                self.isrr[(ind, r)] = True
                self.cycle[(ind, r)] = list(cyc)
                self.pos[(ind, r)] = self._seed(sn, state, ind, r, cyc, is_wrr)

    @staticmethod
    def _cycle_of(sn, ind, r, is_wrr):
        # WRROBIN cycle comes from nodeparam; plain RROBIN derives it from the connection matrix instead.
        if is_wrr:
            from ...state.routing_pointer import wrr_weighted_outlinks
            wol = wrr_weighted_outlinks(sn, ind, r)
            if wol:
                return wol
        nparam = (sn.nodeparam[ind] if getattr(sn, 'nodeparam', None) is not None
                  and ind in sn.nodeparam else None)
        if isinstance(nparam, dict):
            entry = nparam.get(r)
            if isinstance(entry, dict):
                ol = entry.get('outlinks')
                if ol is not None and len(ol) > 0:
                    return list(ol)
            ol = nparam.get('outlinks')
            if ol is not None and r < len(ol):
                cand = np.atleast_1d(ol[r])
                if len(cand) > 0:
                    return [int(x) for x in cand]
        cm = getattr(sn, 'connmatrix', None)
        if cm is not None:
            cm = np.asarray(cm)
            if ind < cm.shape[0]:
                return [int(x) for x in np.flatnonzero(cm[ind, :] > 0)]
        return None

    @staticmethod
    def _seed(sn, state, ind, r, cyc, is_wrr):
        # pointer seeded from the initial state slot when available (full state layout only); native aggregate state starts the cycle at 0.
        R = sn.nclasses
        nvars = getattr(sn, 'nvars', None)
        if nvars is None or not sn.isstateful[ind]:
            return 0
        try:
            st = state[int(sn.nodeToStateful[ind])]
            if st is None:
                return 0
            st = np.atleast_2d(np.asarray(st, dtype=float))
            slot = int(np.sum(np.asarray(nvars[ind]).flatten()[:R + r + 1])) - 1
            if slot < 0 or slot >= st.shape[1]:
                return 0
            v = st[0, slot]
        except (TypeError, KeyError, IndexError):
            return 0
        if is_wrr:
            # a POSITION in the weighted cycle, stored 1-based by
            # _update_rrobin_pointer
            iv = int(round(v))
            return iv - 1 if 1 <= iv <= len(cyc) else 0
        # a destination NODE INDEX
        for j, d in enumerate(cyc):
            if int(d) == int(round(v)):
                return j
        return 0


def _rr_next(rr: _Rr, ind: int, r: int) -> int:
    """Advance the pointer cyclically and return the destination it lands on."""
    cyc = rr.cycle[(ind, r)]
    p = rr.pos[(ind, r)]
    p = 0 if p >= len(cyc) - 1 else p + 1
    rr.pos[(ind, r)] = p
    return int(cyc[p])


# =============================================================================
# Common setup: stoichiometry, initial state, buffers, propensities
# =============================================================================
def _build_nrm_problem(sn):
    """
    Build the stoichiometry matrix S, reaction-to-(node,class) maps, the initial
    aggregate state nvec0 and per-node buffers, the per-station rates / server
    counts, the propensity functions a[j](X, bufs) and their dependency sets D.

    Mirrors the shared front-end of Solver_ssa_nrm / Solver_ssa_nrm_space.
    """
    R = sn.nclasses
    I = sn.nnodes
    epstol = GlobalConstants.Zero

    # The full (MATLAB/JAR) state layout encodes the buffer order and the phase
    # of every job explicitly and needs phasessz/phaseshift to decode; the native
    # Python NetworkStruct instead stores the aggregate per-class marginal, which
    # is exactly R wide (Network._refresh_state writes one entry per class). The
    # layout is therefore decided by the ROW WIDTH, not by phasessz being
    # populated: refresh_struct sets phasessz unconditionally (all-ones when no
    # class is phase-type), so testing it classified every phase-type native
    # model as full-state and handed a per-class marginal to toMarginalAggr,
    # which read it as phase counts and collapsed the whole population onto
    # class 0. A full-state row whose width happens to equal R has all-ones
    # phasessz and no buffer section, so the two decodings coincide there and
    # reading it as native is exact.
    def _is_full_state(state_row) -> bool:
        return int(np.atleast_2d(np.asarray(state_row, dtype=float)).shape[1]) != R
    # mu/phi/pie/phasessz derived from sn.proc (not natively populated by getStruct) and restored afterwards so the mutation does not leak.
    from ..ctmc.handler import _refresh_phase_fields
    _refresh_phase_fields(sn, immediate_as_rate=True)

    smap = _Smap(sn, I, R)
    nph, ph_off, NS = smap.nph, smap.ph_off, smap.NS

    # buffered phase-type service: see _kb/06-solver-catalog.md SSA svcph auxiliary-structure section.
    bufph_class = np.zeros((I, R), dtype=bool)
    for ind in range(I):
        if sn.isstation[ind] and _sched_of(sn, int(sn.nodeToStation[ind])) in _BUFPH_SCHED:
            for r in range(R):
                if int(nph[ind, r]) > 1:
                    bufph_class[ind, r] = True
    bufph_node = np.any(bufph_class, axis=1)
    maxnph = int(np.max(nph)) if nph.size else 1

    # Cache nodes modeled as an immediate state-dependent class switch; see _kb/06-solver-catalog.md SSA Cache nodes section.
    is_cache_node = np.zeros(I, dtype=bool)
    is_cache_read = np.zeros((I, R), dtype=bool)
    cache_retr_dest = np.full((I, R), -1, dtype=int)
    for ind in range(I):
        if sn.nodetype[ind] == NodeType.CACHE:
            is_cache_node[ind] = True
            npm = sn.nodeparam[ind]
            pread_c = _np_get(npm, 'pread', None)
            hitcl_c = np.atleast_1d(_np_get(npm, 'hitclass', []))
            _rci_raw = _np_get(npm, 'retrieval_class_indices', None)
            rci_set = set(int(x) for x in list(_rci_raw)) if _rci_raw is not None and len(_rci_raw) else set()
            rci_set = set(x for x in rci_set if x >= 0)
            for r in range(R):
                pr = pread_c[r] if (pread_c is not None and r < len(pread_c)) else None
                hc = int(hitcl_c[r]) if r < len(hitcl_c) else -1
                if pr is not None:
                    pr_a = np.atleast_1d(np.asarray(pr, dtype=float))
                    # hitclass/missclass are 0-indexed class indices, -1 = none.
                    if pr_a.size > 0 and not np.any(np.isnan(pr_a)) and (hc >= 0 or r in rci_set):
                        is_cache_read[ind, r] = True
            # BEGUN-retrieval destination slot resolved once from rtnodes (the retrieval class's cache row).
            for rc in rci_set:
                row = np.asarray(sn.rtnodes[ind * R + rc, :]).ravel()
                nz = np.nonzero(row > 0)[0]
                if nz.size > 0:
                    dslot = int(nz[0])
                    jnd = dslot // R
                    s = dslot % R
                    cache_retr_dest[ind, rc] = int(ph_off[jnd, s])

    # polling controllers track only the controller state, not per-job phase, so phase-type service at a polling station is rejected rather than mis-modeled.
    poll_info = {}
    is_poll = np.zeros(I, dtype=bool)
    poll_on = False
    for ind in range(I):
        if not sn.isstation[ind]:
            continue
        ist = int(sn.nodeToStation[ind])
        if _sched_of(sn, ist) != SchedStrategy.POLLING:
            continue
        pinfo = polling_info(sn, ind)
        poll_info[ind] = pinfo
        is_poll[ind] = True
        poll_on = True
        for r in range(R):
            if int(nph[ind, r]) > 1:
                raise ValueError(
                    'NRM polling supports exponential service only; station %d '
                    'class %d has phase-type service. Use method="serial".'
                    % (ist, r))

    # --- Stoichiometry & reaction mapping (self-loops included) --------------
    Srows_l: List[np.ndarray] = []
    fromIdx: List[int] = []
    fromIR: List[Tuple[int, int]] = []
    toIdx: List[List[int]] = []
    probIR: List[List[float]] = []
    dep_phase: List[int] = []       # service phase each departure reaction consumes
    is_buf_svc: List[bool] = []     # departure reaction of a buffered-PH class (reads svcph)
    is_cache_rx: List[bool] = []    # cache-access reaction (read -> hit/miss at a cache node)
    cache_hit_slot: List[int] = []  # nvec slot the hit-class job is produced into
    cache_miss_slot: List[int] = [] # nvec slot the miss-class job is produced into

    rtnodes = (np.asarray(sn.rtnodes) if sn.rtnodes is not None
               else np.zeros((I * R, I * R)))
    isslc = np.asarray(sn.isslc).flatten() if sn.isslc is not None else np.zeros(R)

    # departure reactions per (node,class,phase) fire at mu(k)*phi(k); the destination draw carries routing-prob * entry-phase-prob.
    k = 0
    for ind in range(I):
        for r in range(R):
            for kk in range(int(nph[ind, r])):
                fromIR.append((ind, r))
                # cache access: consume the read-class job, resolve hit/miss and contents update at firing via _cache_access; no static routing for the read class.
                if is_cache_read[ind, r]:
                    npm = sn.nodeparam[ind]
                    hitcl = np.atleast_1d(_np_get(npm, 'hitclass', []))
                    misscl = np.atleast_1d(_np_get(npm, 'missclass', []))
                    fromIdx.append(int(ph_off[ind, r]))
                    is_buf_svc.append(False)
                    dep_phase.append(kk)
                    toIdx.append([])
                    probIR.append([])
                    row = np.zeros(NS)
                    row[int(ph_off[ind, r])] = -1.0
                    Srows_l.append(row)
                    is_cache_rx.append(True)
                    cache_hit_slot.append(int(ph_off[ind, int(hitcl[r])]))
                    cache_miss_slot.append(int(ph_off[ind, int(misscl[r])]))
                    k += 1
                    continue
                is_cache_rx.append(False)
                cache_hit_slot.append(0)
                cache_miss_slot.append(0)
                # buffered-PH source departure: only in-service jobs carry a phase (tracked in svcph); nvec's total slot is decremented regardless of which phase completed.
                if bufph_class[ind, r]:
                    fromIdx.append(int(ph_off[ind, r]))
                    is_buf_svc.append(True)
                else:
                    fromIdx.append(int(ph_off[ind, r]) + kk)
                    is_buf_svc.append(False)
                dep_phase.append(kk)
                toIdx.append([])
                probIR.append([])
                row = np.zeros(NS)
                if isslc[r]:
                    row[fromIdx[k]] = -np.inf
                else:
                    row[fromIdx[k]] = -1.0
                    for jnd in range(I):
                        for s in range(R):
                            p = rtnodes[ind * R + r, jnd * R + s]
                            if p <= 0:
                                continue
                            if bufph_class[jnd, s]:
                                # buffered-PH destination arrival collapses to the total slot; service entry and phase are decided at firing from server occupancy and pie.
                                dslot = int(ph_off[jnd, s])
                                toIdx[k].append(dslot)
                                probIR[k].append(float(p))
                                row[dslot] += p
                                continue
                            pentry = _entry_probs(sn, jnd, s, int(nph[jnd, s]))
                            for ke in range(int(nph[jnd, s])):
                                if pentry[ke] <= 0:
                                    continue
                                dslot = int(ph_off[jnd, s]) + ke
                                toIdx[k].append(dslot)
                                probIR[k].append(float(p * pentry[ke]))
                                row[dslot] += p * pentry[ke]
                Srows_l.append(row)
                k += 1
    n_dep_only = k          # departure reactions occupy 0..n_dep_only-1

    # phase-transition reactions move a job between its own service phases (D0 off-diagonal); they never leave the node.
    phase_rx = []           # (reaction index, from phase, to phase, rate)
    for ind in range(I):
        if not sn.isstation[ind]:
            continue
        ist = int(sn.nodeToStation[ind])
        for r in range(R):
            if int(nph[ind, r]) <= 1:
                continue
            proc_ir = (sn.proc[ist][r] if (sn.proc is not None and ist < len(sn.proc)
                                           and r < len(sn.proc[ist])) else None)
            # _refresh_phase_fields normalizes every enabled process to
            # [D0, D1]; anything else carries no D0 off-diagonal to expand.
            if not isinstance(proc_ir, (list, tuple)) or len(proc_ir) < 2:
                continue
            D0 = np.atleast_2d(np.asarray(proc_ir[0], dtype=float))
            if D0.shape[0] < int(nph[ind, r]) or np.any(np.isnan(D0)):
                continue
            for ka in range(int(nph[ind, r])):
                for kb in range(int(nph[ind, r])):
                    if ka == kb or D0[ka, kb] <= 0:
                        continue
                    fromIR.append((ind, r))
                    dep_phase.append(0)
                    toIdx.append([])
                    probIR.append([])
                    row = np.zeros(NS)
                    if bufph_class[ind, r]:
                        # buffered-PH phase transition changes svcph only (nvec, the class total, is untouched), so its stoichiometry column is all zeros.
                        fromIdx.append(int(ph_off[ind, r]))
                    else:
                        fromIdx.append(int(ph_off[ind, r]) + ka)
                        row[int(ph_off[ind, r]) + ka] = -1.0
                        row[int(ph_off[ind, r]) + kb] = 1.0
                    Srows_l.append(row)
                    phase_rx.append((k, ka, kb, float(D0[ka, kb])))
                    k += 1

    Srows = (np.asarray(Srows_l) if Srows_l else np.zeros((0, NS)))
    num = k                 # departure + phase-transition reactions
    is_phase = np.zeros(num, dtype=bool)
    phase_rate = np.zeros(num)
    phase_from = np.zeros(num, dtype=int)   # source phase of a phase-transition reaction
    phase_to = np.zeros(num, dtype=int)     # target phase of a phase-transition reaction
    for (rx, ka, kb, rate) in phase_rx:
        is_phase[rx] = True
        phase_rate[rx] = rate
        phase_from[rx] = ka
        phase_to[rx] = kb

    # reneging: waiting jobs abandon at rate impatienceMu; modeled as an extra all-zero-column reaction (a bare -1 at the source slot), excluded from throughput; source slot is the class's one phase since only INF/PS admit non-exponential service.
    n_dep = k               # departure + phase-transition reactions
    renege_rx = []          # (reaction index, node, class, mu)
    ic = getattr(sn, 'impatienceClass', None)
    if ic is not None and np.size(ic) > 0:
        from ....lang.base import ImpatienceType
        ic = np.asarray(ic)
        mu_all = np.asarray(sn.impatienceMu, dtype=float)
        extra = []
        for ist in range(sn.nstations):
            ind = int(sn.stationToNode[ist])
            for r in range(R):
                if ic[ist, r] == int(ImpatienceType.RENEGING) and mu_all[ist, r] > 0:
                    row = np.zeros(NS)
                    slot = int(ph_off[ind, r])
                    row[slot] = -1.0   # abandons and leaves the system
                    extra.append(row)
                    fromIdx.append(slot)
                    fromIR.append((ind, r))
                    toIdx.append([])
                    probIR.append([])
                    renege_rx.append((k, ind, r, float(mu_all[ist, r])))
                    k += 1
        if extra:
            Srows = np.vstack([Srows, np.asarray(extra)])

    # retrial: orbiting job retries at retrialMu, succeeding only if a server is free; a no-op otherwise, so the reaction has an all-zero stoichiometry column and its dependency set is supplied explicitly.
    retry_rx = []           # (reaction index, node, class, mu)
    rproc = getattr(sn, 'retrialProc', None)
    if rproc is not None:
        rmu = np.asarray(getattr(sn, 'retrialMu', np.zeros((sn.nstations, R))), dtype=float)
        extra = []
        for ist in range(sn.nstations):
            if ist >= len(rproc) or not any(p is not None for p in rproc[ist]):
                continue
            ind = int(sn.stationToNode[ist])
            for r in range(R):
                if rmu[ist, r] > 0:
                    extra.append(np.zeros(NS))   # a retry moves no job between nodes
                    fromIdx.append(int(ph_off[ind, r]))
                    fromIR.append((ind, r))
                    toIdx.append([])
                    probIR.append([])
                    retry_rx.append((k, ind, r, float(rmu[ist, r])))
                    k += 1
        if extra:
            Srows = np.vstack([Srows, np.asarray(extra)])

    # polling switchover is appended as a zero-stoichiometry reaction per node with a timed switchover; immediate switchovers fold into the enclosing event.
    poll_sw_node_of = {}     # reaction index -> polling node
    if poll_on:
        extra = []
        for ind in range(I):
            pinfo = poll_info.get(ind)
            if pinfo is None or not bool(np.any(pinfo['has_sw'])):
                continue
            extra.append(np.zeros(NS))   # a switchover moves no job
            fromIdx.append(int(ph_off[ind, 0]))   # sentinel slot; never read as a class
            fromIR.append((ind, 0))
            toIdx.append([])
            probIR.append([])
            poll_sw_node_of[k] = ind
            k += 1
        if extra:
            Srows = np.vstack([Srows, np.asarray(extra)])

    # pad phase flags/rates over the appended renege/retry reactions so every per-reaction array shares one index.
    is_phase = np.concatenate([is_phase, np.zeros(k - num, dtype=bool)])
    phase_rate = np.concatenate([phase_rate, np.zeros(k - num)])
    phase_from = np.concatenate([phase_from, np.zeros(k - num, dtype=int)])
    phase_to = np.concatenate([phase_to, np.zeros(k - num, dtype=int)])
    dep_phase = dep_phase + [0] * (k - num)
    is_buf_svc = np.concatenate([np.asarray(is_buf_svc, dtype=bool),
                                 np.zeros(k - len(is_buf_svc), dtype=bool)])
    is_cache_rx = np.concatenate([np.asarray(is_cache_rx, dtype=bool),
                                  np.zeros(k - len(is_cache_rx), dtype=bool)])
    cache_hit_slot = np.concatenate([np.asarray(cache_hit_slot, dtype=int),
                                     np.zeros(k - len(cache_hit_slot), dtype=int)])
    cache_miss_slot = np.concatenate([np.asarray(cache_miss_slot, dtype=int),
                                      np.zeros(k - len(cache_miss_slot), dtype=int)])

    is_renege = np.zeros(k, dtype=bool)
    renege_mu = np.zeros(k)
    for (rx, _ind, _r, mu) in renege_rx:
        is_renege[rx] = True
        renege_mu[rx] = mu
    is_retry = np.zeros(k, dtype=bool)
    retry_mu = np.zeros(k)
    for (rx, _ind, _r, mu) in retry_rx:
        is_retry[rx] = True
        retry_mu[rx] = mu
    is_poll_sw = np.zeros(k, dtype=bool)
    poll_sw_node = np.full(k, -1, dtype=int)
    for rx, ind in poll_sw_node_of.items():
        is_poll_sw[rx] = True
        poll_sw_node[rx] = ind

    S = Srows.T  # NS x numReactions

    # --- Initial state vector and per-node buffers ---------------------------
    nvec0 = np.zeros(NS)
    buffers0: List[deque] = [deque() for _ in range(I)]
    # a cache node's buffers slot carries cache CONTENTS (not a queueing buffer); any valid ordered placement is a correct warm start.
    for ind in range(I):
        if is_cache_node[ind]:
            npm = sn.nodeparam[ind]
            tcc = int(_np_get(npm, 'total_cache_capacity', 0) or 0)
            if tcc <= 0:
                tcc = int(np.sum(np.atleast_1d(_np_get(npm, 'itemcap', [0]))))
            # with a retrieval system the contents are followed by a per-item occupancy bitmap (State.spaceCache), started empty.
            rsc = _np_get(npm, 'retrieval_system_capacity', 0)
            if rsc is not None and np.any(np.atleast_1d(np.asarray(rsc, dtype=float)) > 0):
                nitems = int(_np_get(npm, 'nitems', 0) or 0)
                buffers0[ind] = list(range(1, tcc + 1)) + [0] * nitems
            else:
                buffers0[ind] = list(range(1, tcc + 1))
    # svcph0[ind][r,k]: in-service phase multiset of each buffered-PH node, filled after the waiting buffer is known; see _kb/06-solver-catalog.md SSA svcph section.
    svcph0: List[Optional[np.ndarray]] = [
        (np.zeros((R, maxnph)) if bufph_node[ind] else None) for ind in range(I)]
    for ind in range(I):
        if sn.isstateful[ind]:
            isf = int(sn.nodeToStateful[ind])
            state_i = sn.state[isf]
            if state_i is None:
                raise RuntimeError("State matrix for stateful node %d is null" % ind)
            state_i = np.atleast_2d(np.asarray(state_i, dtype=float))

            full_state = _is_full_state(state_i)
            if full_state:
                _, nir = toMarginalAggr(sn, ind, state_i)
            else:
                # native representation: the row already holds per-class counts
                nir = state_i[0, :R]
            nir = np.atleast_1d(np.asarray(nir).flatten())

            for r in range(R):
                v = nir[r] if r < nir.shape[0] else 0.0
                if np.isinf(v):
                    if sn.nodetype[ind] == NodeType.SOURCE:
                        v = 1.0
                    else:
                        raise RuntimeError("Infinite population error.")
                # class population spread across phases via the entry distribution pie (the phase a job starts service in); buffered-PH keeps its total in slot 0, detail lives in svcph0.
                if int(nph[ind, r]) <= 1 or bufph_class[ind, r]:
                    nvec0[int(ph_off[ind, r])] = v
                else:
                    pe = _entry_probs(sn, ind, r, int(nph[ind, r]))
                    left = v
                    for ke in range(int(nph[ind, r])):
                        if ke == int(nph[ind, r]) - 1:
                            take = left
                        else:
                            take = min(left, float(round(v * pe[ke])))
                        nvec0[int(ph_off[ind, r]) + ke] = take
                        left -= take

            # buffers populated from the raw state vector where the full layout carries it; native representation starts buffers empty (standard init places jobs at reference stations).
            ist = int(sn.nodeToStation[ind])
            # Only stations have FCFS/LCFS scheduling; skip non-station stateful
            # nodes such as RROBIN dispatchers/Routers and Caches (ist == -1).
            if ist >= 0 and not full_state and _sched_of(sn, ist) in _LIST_SCHED:
                # a list station's rate law needs len(buf)==total; the native aggregate state seeds it class-ascending (one of several equivalent uniqueperms rows of the same ergodic chain, so the choice only affects the initial transient).
                for r2 in range(R):
                    buffers0[ind].extend(
                        [r2 + 1] * int(round(_class_pop(nvec0, smap, ind, r2))))
            if ist >= 0 and full_state and _sched_of(sn, ist) in _INIT_BUFFERED_SCHED:
                sumK = int(np.sum(np.asarray(sn.phasessz[ist]).flatten()))
                sumNvars = int(np.sum(np.asarray(sn.nvars[ind]).flatten())) if sn.nvars is not None else 0
                bufCols = state_i.shape[1] - sumK - sumNvars
                if _sched_of(sn, ist) in _LIST_SCHED:
                    # PAS/OI stores the full ordered list left-aligned, oldest first, already in NRM order (copied verbatim).
                    for pos in range(bufCols):
                        classId = int(round(state_i[0, pos]))
                        if 1 <= classId <= R:
                            buffers0[ind].append(classId)
                elif _sched_of(sn, ist) in _COUNT_BUFFERED_SCHED:
                    # SIRO/SEPT/LEPT's unordered per-class counts are expanded into the NRM's ordered list; order within it is immaterial for these disciplines.
                    for r2 in range(min(R, bufCols)):
                        buffers0[ind].extend([r2 + 1] * int(round(state_i[0, r2])))
                else:
                    # FCFS/HOL/LCFS keep an ordered list of class ids
                    for pos in range(bufCols):
                        classId = int(round(state_i[0, pos]))
                        if 1 <= classId <= R:
                            buffers0[ind].append(classId)  # addLast
                        # classId == 0 means empty position, skip

            # in-service phase multiset seeded from class-total-minus-waiting, phases drawn from pie; only in-service jobs get a phase.
            if bufph_node[ind]:
                for r in range(R):
                    waiting_r = sum(1 for cl in buffers0[ind] if cl == r + 1)
                    tot_r = float(nir[r]) if r < nir.shape[0] else 0.0
                    insvc_r = max(0.0, tot_r - waiting_r)
                    if int(nph[ind, r]) <= 1:
                        svcph0[ind][r, 0] = insvc_r
                    else:
                        pe = _entry_probs(sn, ind, r, int(nph[ind, r]))
                        left = insvc_r
                        for ke in range(int(nph[ind, r])):
                            if ke == int(nph[ind, r]) - 1:
                                take = left
                            else:
                                take = min(left, float(round(insvc_r * pe[ke])))
                            svcph0[ind][r, ke] = take
                            left -= take

    # --- Per-node rates and server counts ------------------------------------
    mi = np.zeros(I)
    rates = np.zeros((I, R))
    for ind in range(I):
        if sn.isstation[ind]:
            ist = int(sn.nodeToStation[ind])
            for r in range(R):
                muir = sn.rates[ist, r]
                if not np.isnan(muir):
                    rates[ind, r] = muir
            mi[ind] = sn.nservers[ist]
        else:
            for r in range(R):
                rates[ind, r] = GlobalConstants.Immediate
            mi[ind] = GlobalConstants.MaxInt
        if np.isinf(mi[ind]):
            mi[ind] = GlobalConstants.MaxInt

    # limited load-dependent scaling multiplies the aggregate service rate by ntot; default 1 (inert) without load dependence.
    lldMat = getattr(sn, 'lldscaling', None)
    if lldMat is not None and not isinstance(lldMat, np.ndarray):
        lldMat = np.asarray(lldMat)
    if lldMat is not None and (lldMat.ndim < 2 or lldMat.size == 0):
        lldMat = None
    lldlimit = lldMat.shape[1] if lldMat is not None else 0

    def _lldfac(ist, ntot):
        if lldMat is None or ntot < 1:
            return 1.0
        idx = int(min(int(round(ntot)), lldlimit)) - 1
        if 0 <= ist < lldMat.shape[0] and 0 <= idx < lldMat.shape[1]:
            return float(lldMat[ist, idx])
        return 1.0

    # class-dependent scaling maps the per-class population vector to a 1xR rate-scaling vector, evaluated per firing.
    cdCell = getattr(sn, 'cdscaling', None)
    # joint-dependent (non-product-form eta_i) scaling; folded multiplicatively
    # into the same firing-rate factor as class dependence.
    jdCell = getattr(sn, 'jdscaling', None)

    def _one_fac(cell, ist, nvec, r):
        if cell is None:
            return 1.0
        fn = (cell.get(ist, None) if isinstance(cell, dict)
              else (cell[ist] if ist < len(cell) else None))
        if not callable(fn):
            return 1.0
        v = np.atleast_1d(np.asarray(
            fn(np.asarray(nvec, dtype=float).flatten()), dtype=float)).flatten()
        return float(v[min(int(r), v.size - 1)])

    def _cdfac(ist, nvec, r):
        # class-r rate-scaling factor: product of the class- and joint-dependence
        # handles' class-r components, 1 when neither is declared.
        return _one_fac(cdCell, ist, nvec, r) * _one_fac(jdCell, ist, nvec, r)

    # Normalized per-class weights of the DPS/GPS-family stations.
    wnorm = _build_sched_weights(sn, R)

    # service-process event rate: mu(k)*phi(k) for a departure, D0(k,k') for a phase change; reduces to the plain exponential rate for a single-phase class.
    rate_of = np.zeros(len(fromIdx))
    for j in range(len(fromIdx)):
        ind, r = fromIR[j]
        if j < num and is_phase[j]:
            rate_of[j] = phase_rate[j]
        elif j < num and sn.isstation[ind] and int(nph[ind, r]) > 1:
            ist = int(sn.nodeToStation[ind])
            kk = dep_phase[j]
            rate_of[j] = (float(np.atleast_1d(sn.mu[ist][r])[kk])
                          * float(np.atleast_1d(sn.phi[ist][r])[kk]))
        else:
            rate_of[j] = rates[ind, r]

    # --- Propensity functions a[j](X, bufs) ----------------------------------
    # Every share helper below is a CLASS-level law and is reused unchanged on
    # the class counts; phase expansion enters solely as the kir/nir factor,
    # which splits that class share across the class's phases exactly as
    # State.afterEventStation does.
    def make_prop(jj):
        ind, rr = fromIR[jj]
        fidx = fromIdx[jj]
        base = rate_of[jj]
        kfrac = (lambda X: _kir_frac(X, fidx, smap, ind, rr))
        if not sn.isstation[ind]:
            return lambda X, bufs, svc: (base * kfrac(X)
                                    * min(1.0, _class_pop(X, smap, ind, rr)))
        ist = int(sn.nodeToStation[ind])
        sched = _sched_of(sn, ist)

        # buffered-PH departure/phase-transition rates read the in-service multiset svc, not nvec; see _kb/06-solver-catalog.md SSA svcph section.
        if bufph_class[ind, rr]:
            kk_ph = int(phase_from[jj]) if is_phase[jj] else int(dep_phase[jj])
            def f_bufph(X, bufs, svc):
                npop = _class_counts(X, smap, ind)
                return (base * float(svc[ind][rr, kk_ph])
                        * _lldfac(ist, float(np.sum(npop))) * _cdfac(ist, npop, rr))
            return f_bufph

        # a phase-transition reaction is throttled by the same scheduling rate law as a departure (PS-family jobs advance phases only at their shared rate).
        if sched == SchedStrategy.EXT:
            # a Source's arrival rate is state-independent; applying the kir/nir share would latch it to zero once the fictitious token count crosses 0.
            return lambda X, bufs, svc: base
        if sched == SchedStrategy.INF:
            return lambda X, bufs, svc: (base * kfrac(X) * _class_pop(X, smap, ind, rr)
                                    * _cdfac(ist, _class_counts(X, smap, ind), rr))
        if sched in (SchedStrategy.PS, SchedStrategy.LPS):
            # LPS shares the PS rate law in State.afterEventStation: the
            # sharing limit is the server count, so min(ni,c) already covers it.
            if R == 1:
                def f_ps1(X, bufs, svc):
                    nir = _class_pop(X, smap, ind, rr)
                    return (base * kfrac(X) * min(mi[ind], nir) * _lldfac(ist, nir)
                            * _cdfac(ist, _class_counts(X, smap, ind), rr))
                return f_ps1
            def f_ps(X, bufs, svc):
                npop = _class_counts(X, smap, ind)
                total = epstol + float(np.sum(npop))
                return (base * kfrac(X) * (npop[rr] / total) * min(mi[ind], total)
                        * _lldfac(ist, total) * _cdfac(ist, npop, rr))
            return f_ps
        if sched == SchedStrategy.DPS:
            wrow = wnorm[ist]
            def f_dps(X, bufs, svc):
                npop = _class_counts(X, smap, ind)
                return (base * kfrac(X) * _dpsshare(wrow, npop, rr)
                        * _lldfac(ist, float(np.sum(npop))) * _cdfac(ist, npop, rr))
            return f_dps
        if sched == SchedStrategy.GPS:
            wrow = wnorm[ist]
            def f_gps(X, bufs, svc):
                npop = _class_counts(X, smap, ind)
                return (base * kfrac(X) * _gpsshare(wrow, npop, rr)
                        * _lldfac(ist, float(np.sum(npop))) * _cdfac(ist, npop, rr))
            return f_gps
        if sched in _LIST_SCHED:
            # PAS/OI class-r departure rate = summed Delta_mu over positions whose pass-and-swap ejects a class-r job; see after_event_station_pas.
            mu_fun, G = _oi_param(sn, ind)
            return lambda X, bufs, svc: _oirate(mu_fun, G, bufs[ind], rr)
        if sched in _BUFFERED_SCHED:
            def f_fcfs(X, bufs, svc):
                waiting = sum(1 for cl in bufs[ind] if cl == rr + 1)
                npop = _class_counts(X, smap, ind)
                in_service = npop[rr] - waiting
                if in_service <= 0:
                    return 0.0
                total = float(np.sum(npop))
                return (base * kfrac(X) * in_service * _lldfac(ist, total)
                        * _cdfac(ist, npop, rr))
            return f_fcfs
        if sched == SchedStrategy.POLLING:
            # A polling station has a single server that serves exactly one job,
            # of the class its controller currently attends. The class-r service
            # departure therefore fires only while the controller is SERVING
            # class r, at the plain service rate of the one job in service --
            # never scaled by the class population, since the other class-r jobs
            # wait in the buffer for the server to come back. No kir/nir share
            # is applied: service at a polling station is exponential (phase-type
            # is rejected in _build_nrm_problem), so the class is single-phase.
            def f_poll(X, bufs, svc):
                npop = _class_counts(X, smap, ind)
                return (base * _poll_serve_gate(bufs[ind], rr)
                        * _lldfac(ist, float(np.sum(npop)))
                        * _cdfac(ist, npop, rr))
            return f_poll
        if sched == SchedStrategy.PSPRIO:
            def f_psprio(X, bufs, svc):
                npop = _class_counts(X, smap, ind)
                return (base * kfrac(X) * _psprioshare(sn, npop, rr, mi[ind])
                        * _lldfac(ist, _prio_pop(sn, npop, rr, mi[ind]))
                        * _cdfac(ist, npop, rr))
            return f_psprio
        if sched == SchedStrategy.DPSPRIO:
            wrow = wnorm[ist]
            def f_dpsprio(X, bufs, svc):
                npop = _class_counts(X, smap, ind)
                return (base * kfrac(X) * _dpsprioshare(sn, wrow, npop, rr, mi[ind])
                        * _lldfac(ist, _prio_pop(sn, npop, rr, mi[ind]))
                        * _cdfac(ist, _prio_vec(sn, npop, rr, mi[ind]), rr))
            return f_dpsprio
        if sched == SchedStrategy.GPSPRIO:
            wrow = wnorm[ist]
            def f_gpsprio(X, bufs, svc):
                npop = _class_counts(X, smap, ind)
                return (base * kfrac(X) * _gpsprioshare(sn, wrow, npop, rr, mi[ind])
                        * _lldfac(ist, _prio_pop(sn, npop, rr, mi[ind]))
                        * _cdfac(ist, _prio_vec(sn, npop, rr, mi[ind]), rr))
            return f_gpsprio
        return lambda X, bufs, svc: base * kfrac(X) * _class_pop(X, smap, ind, rr)

    a = [make_prop(j) for j in range(num)]

    # --- Reneging propensities ----------------------------------------------
    # Only the jobs actually waiting can abandon, so the rate is the class-r
    # buffer occupancy, the same quantity the FCFS-family rate law relies on.
    def make_renege(j):
        ind, rr = fromIR[j]
        mu = renege_mu[j]
        return lambda X, bufs, svc: mu * sum(1 for cl in bufs[ind] if cl == rr + 1)

    # Retrial propensities: only jobs actually in orbit retry, and only a free
    # server admits them. The orbit is the buffer, and in-service is population
    # minus orbit occupancy -- the "total > mi" test the other policies use is
    # INVALID here, because a departure does not promote, so the orbit can be
    # full while servers sit idle.
    def make_retry(j):
        ind, rr = fromIR[j]
        mu = retry_mu[j]

        def f_retry(X, bufs, svc):
            orbit = sum(1 for cl in bufs[ind] if cl == rr + 1)
            if orbit <= 0:
                return 0.0
            in_service = float(np.sum(X[smap.node_slots(ind)])) - len(bufs[ind])
            return mu * orbit if in_service < mi[ind] else 0.0
        return f_retry

    # Polling switchover propensity: the total leaving rate of the current
    # switchover phase, 0 unless the controller is SWITCHING. The choice between
    # advancing to another phase and absorbing (arriving at the target buffer) is
    # resolved at firing time in the run loop, exactly as a routed departure
    # resolves its destination after it fires.
    def make_poll_sw(j):
        ind = int(poll_sw_node[j])
        pinfo = poll_info[ind]
        return lambda X, bufs, svc: _poll_sw_rate(bufs[ind], pinfo)

    for j in range(n_dep, len(fromIdx)):
        if is_poll_sw[j]:
            a.append(make_poll_sw(j))
        elif is_retry[j]:
            a.append(make_retry(j))
        else:
            a.append(make_renege(j))

    # FCR under DROP destroys the refused job at firing (never scales the propensity, which would wrongly keep it at the source); see _kb/06-solver-catalog.md SSA FCR update note.
    fcr = _Fcr(sn)

    # --- Propensity dependencies D[k] ----------------------------------------
    D: List[List[int]] = [list() for _ in range(S.shape[1])]
    for kk in range(S.shape[1]):
        if is_retry[kk]:
            # a retry has an all-zero stoichiometry column, so its dependency set (every reaction at the affected node) is supplied explicitly rather than derived.
            node_k = fromIR[kk][0]
            D[kk] = [j for j in range(len(fromIdx))
                     if int(smap.node[fromIdx[j]]) == node_k]
            continue
        J = [i for i in range(S.shape[0]) if S[i, kk] != 0.0]
        vecd = []
        for pos in J:
            # affected slots decoded through the phase-expansion slot map (never arithmetically), collecting every phase-slot of each affected node.
            ind = int(smap.node[pos])
            sl = smap.node_slots(ind)
            vecd.extend(range(sl.start, sl.stop))
        if vecd:
            vecd_unique = list(dict.fromkeys(vecd))
            vecs = []
            seen = set()
            for u in vecd_unique:
                for kk2 in range(S.shape[1]):
                    # no fcr.on propensity widening: regions no longer gate rates (see the FCR note above); admission resolves at firing time on the drawn destination.
                    touches = S[u, kk2] < 0.0
                    if touches and kk2 not in seen:
                        seen.add(kk2)
                        vecs.append(kk2)
            D[kk] = vecs

    # a retry's all-zero stoichiometry column excludes it from the sign-based dependency derivation, so it must be refreshed by hand whenever its node changes.
    for j in np.flatnonzero(is_retry):
        indj = fromIR[j][0]
        slots = smap.node_slots(indj)
        for kk in range(S.shape[1]):
            if np.any(S[slots, kk] != 0.0) or fromIR[kk][0] == indj:
                if j not in D[kk]:
                    D[kk].append(int(j))

    # Remove self-loop (-inf) markings now that D accounts for them
    S[np.isinf(S)] = 0.0

    return (S, D, a, fromIdx, fromIR, nvec0, buffers0, mi, fcr, _Balk(sn),
            is_renege, _Sig(sn), is_retry, _Rr(sn, sn.state), smap, n_dep_only,
            is_phase, is_poll, poll_info, is_poll_sw, poll_sw_node,
            svcph0, bufph_node, is_buf_svc, dep_phase, phase_from, phase_to,
            is_cache_rx, cache_hit_slot, cache_miss_slot, is_cache_node, is_cache_read,
            cache_retr_dest)


def _build_routing(S):
    """Per-reaction routing for inverse-CDF destination sampling.

    Mirrors MATLAB next_reaction_method_direct: P = S with negative entries
    incremented by 1, so each firing moves exactly one job to a single
    stochastically chosen destination (keeping the marginal state integer).
    Without this a fractional routing column would be applied in full, draining
    immediate pass-through nodes (e.g. ClassSwitch) below zero and biasing
    station queue lengths.
    """
    num_reactions = S.shape[1]
    rows = S.shape[0]
    dest_row = [None] * num_reactions
    cdf = [None] * num_reactions
    nnzP = [0] * num_reactions
    for k in range(num_reactions):
        dest = []
        pr = []
        for row in range(rows):
            v = S[row, k]
            p = v + 1.0 if v < 0.0 else v  # P = S; P(P<0) += 1
            if p > 0.0:
                dest.append(row)
                pr.append(p)
        nnzP[k] = len(dest)
        dest_row[k] = np.asarray(dest, dtype=int)
        cdf[k] = np.cumsum(pr) if pr else np.zeros(0)
    return nnzP, dest_row, cdf


def _sq_param(sn, ind, r):
    """d of a SQ-routed (node, class), from sn.nodeparam."""
    np_ = getattr(sn, 'nodeparam', None)
    if not np_ or ind not in np_ or not isinstance(np_[ind], dict):
        return 2
    entry = np_[ind].get(r)
    if not isinstance(entry, dict):
        return 2
    k = entry.get('d', 2)     # sub_sq default when nodeparam carries no d
    return int(k) if k else 2


def _build_sq_k(sn, fromIR, routing):
    """
    Number of candidates each SQ-routed reaction draws, 0 when it does
    not route with SQ.

    Power-of-k choices samples k candidates uniformly WITH replacement and joins
    the one holding the smallest total population, ties broken by first
    occurrence in the sampled tuple. This is the sampled form of the marginal
    enumerated by sub_sq in MNetwork.refreshRoutingMatrix and of LDES's
    selectSQDestination; drawing directly is equivalent for a
    simulator and avoids enumerating the ndest^d tuples. The m>0 variant
    forces the previous pick as the last candidate, which needs a
    per-(node,class) memory the reaction network does not carry, so those models
    are routed to the serial engine and never reach here.
    """
    from ....constants import RoutingStrategy
    nnzP, dest_row, _ = routing
    kch_val = int(RoutingStrategy.SQ.value)
    out = [0] * len(nnzP)
    if getattr(sn, 'routing', None) is None:
        return out
    for k in range(len(nnzP)):
        if nnzP[k] <= 1:
            continue
        ind, r = fromIR[k]
        v = sn.routing[ind, r]
        v = int(v.value) if hasattr(v, 'value') else int(v)
        if v != kch_val:
            continue
        kval = _sq_param(sn, ind, r)
        out[k] = max(1, min(kval, len(dest_row[k])))
    return out


def _build_jsq_flags(sn, fromIR, routing):
    """Flag reactions whose source class routes with JSQ.

    Such reactions select the destination with the smallest target-node
    population at firing time instead of inverse-CDF sampling.
    """
    from ....constants import RoutingStrategy
    nnzP, _, _ = routing
    jsq_val = int(RoutingStrategy.JSQ.value)
    flags = [False] * len(nnzP)
    if getattr(sn, 'routing', None) is None:
        return flags
    for k in range(len(nnzP)):
        if nnzP[k] > 1:
            ind, r = fromIR[k]
            v = sn.routing[ind, r]
            v = int(v.value) if hasattr(v, 'value') else int(v)
            if v == jsq_val:
                flags[k] = True
    return flags


def _fire_reaction(kfire, nvec, S, routing, src_row, jsq_flags=None, R=1,
                   fcr=None, fromIR=None, sq_k=None, balk=None,
                   sig=None, buffers=None, mi=None, rr=None, smap=None,
                   fcr_buf=None, sn=None):
    """Apply one firing; sample a single destination when there are several.

    JSQ-routed reactions join the candidate whose destination node holds the
    smallest total population (each candidate evaluated on its own queue,
    never the routing node's; ties split uniformly). All other reactions use
    inverse-CDF sampling over the static routing probabilities.

    A finite capacity region does NOT filter the routing draw: routing picks
    the destination first and the region decides admission at the destination's
    entry afterwards, dropping the job on refusal.

    Returns the destination state-row index (-1 for self-loop / no move).
    """
    nnzP, dest_row, cdf = routing
    if nnzP[kfire] > 1:
        cand = dest_row[kfire]
        if jsq_flags is not None and jsq_flags[kfire]:
            npop = np.full(len(cand), np.inf)
            for x in range(len(cand)):
                jnd = int(smap.node[int(cand[x])])
                npop[x] = float(np.sum(nvec[smap.node_slots(jnd)]))
            amins = np.flatnonzero(npop == npop.min())
            sel = int(amins[int(np.random.random() * len(amins))]) if len(amins) > 1 else int(amins[0])
        elif rr is not None and rr.on and (fromIR[kfire] in rr.isrr):
            # round-robin: advance the pointer, then take the destination it lands on.
            src_node, src_class = fromIR[kfire]
            jnd = _rr_next(rr, src_node, src_class)
            # a phase-type destination needs its entry phase drawn from pentry among the pointer-fixed node's candidates (RUN-10: fixing phase 0 biases the service time).
            cd = cdf[kfire]
            matches = [x for x in range(len(cand))
                       if int(smap.node[int(cand[x])]) == jnd]
            if not matches:
                raise ValueError(
                    'Round-robin selected node %d, which is not a routing '
                    'destination of node %d.' % (jnd, src_node))
            if len(matches) == 1:
                sel = matches[0]
            else:
                w = np.array([cd[x] - (cd[x - 1] if x > 0 else 0.0)
                              for x in matches], dtype=float)
                tot = float(w.sum())
                if tot <= 0.0:
                    sel = matches[0]
                else:
                    u = np.random.random() * tot
                    acc = 0.0
                    sel = matches[-1]
                    for i, x in enumerate(matches):
                        acc += w[i]
                        if acc > u:
                            sel = x
                            break
        elif sq_k is not None and sq_k[kfire] > 0:
            # SQ (k-choices): draw k candidates uniformly with replacement, keep the least loaded, ties broken by first occurrence.
            npop = np.empty(len(cand))
            for x in range(len(cand)):
                jnd = int(smap.node[int(cand[x])])
                npop[x] = float(np.sum(nvec[smap.node_slots(jnd)]))
            sel = 0
            best = np.inf
            draws = min(sq_k[kfire], len(cand))
            for _t in range(draws):
                x = int(np.random.random() * len(cand))
                if npop[x] < best:
                    best = npop[x]
                    sel = x
        else:
            u = np.random.random()
            cd = cdf[kfire]
            sel = len(cd) - 1
            for x in range(len(cd)):
                if cd[x] > u:
                    sel = x
                    break
        # balking is decided on the pre-arrival population; a balked job is lost (source still releases it).
        dest = int(dest_row[kfire][sel])
        lost = balk is not None and balk.on and _balk_draw(balk, nvec, dest, R, smap)
        # an open arrival at a full physically-capped destination is lost like a balk; mirrors State.afterEventStation's finite-capacity gate.
        if not lost and _capacity_loss(sn, nvec, dest, R, smap):
            lost = True
        # region refusal is decided on the pre-arrival population: DROP loses the job, WAITQ parks it head-of-line in the region FIFO; the source still departs either way.
        if not lost and fcr is not None and fcr.on:
            dst_n = int(smap.node[dest])
            dst_c = int(smap.cls[dest])
            fref = _fcr_refusing_region(fcr, nvec, fromIR[kfire][0], fromIR[kfire][1], dst_n, dst_c, R, smap)
            if fref >= 0:
                lost = True
                if fcr_buf is not None and fcr.waitq[fref, dst_c]:
                    fcr_buf[fref].append(dst_n * R + dst_c)
        if lost:
            nvec[src_row] -= 1.0
            return -1
        if sig is not None and sig.on and _sig_is_signal_arrival(sig, dest, R, smap):
            # the signal is annihilated on arrival: it never joins the station
            nvec[src_row] -= 1.0
            _sig_apply(sig, nvec, buffers, dest, R, mi, smap)
            return -1
        nvec[src_row] -= 1.0
        nvec[dest] += 1.0
        return dest
    if nnzP[kfire] == 1:
        dest = int(dest_row[kfire][0])
        lost = balk is not None and balk.on and _balk_draw(balk, nvec, dest, R, smap)
        # An open arrival at a full physically-capped destination is lost, as a
        # balked one is (renege/retry carry no destination and never reach here).
        if not lost and _capacity_loss(sn, nvec, dest, R, smap):
            lost = True
        # region gate applies to single-destination departures too; renege/retry columns have no destination and never reach this check.
        if not lost and fcr is not None and fcr.on:
            dst_n = int(smap.node[dest])
            dst_c = int(smap.cls[dest])
            fref = _fcr_refusing_region(fcr, nvec, fromIR[kfire][0], fromIR[kfire][1], dst_n, dst_c, R, smap)
            if fref >= 0:
                lost = True
                if fcr_buf is not None and fcr.waitq[fref, dst_c]:
                    fcr_buf[fref].append(dst_n * R + dst_c)
        if lost:
            nvec[src_row] -= 1.0
            return -1
        if sig is not None and sig.on and _sig_is_signal_arrival(sig, dest, R, smap):
            nvec[src_row] -= 1.0
            _sig_apply(sig, nvec, buffers, dest, R, mi, smap)
            return -1
        nvec += S[:, kfire]
        return dest
    nvec += S[:, kfire]
    return -1


def _draw_entry_phase(sn, jnd, s, nphjs):
    """Sample the service phase a class-s job starts in at node jnd from its
    entry distribution pie. A single-phase class always enters phase 0."""
    if nphjs <= 1:
        return 0
    return _draw_from_dist(_entry_probs(sn, jnd, s, nphjs))


def _update_buffers(kfire, nvec, buffers, fromIR, destPos, mi, R, sn, smap: _Smap,
                    svcph=None, bufph_node=None, is_buf_svc=None, dep_phase=None):
    """Maintain the ordered per-node buffers after reaction kfire fires. At a
    buffered-PH node the same events also move jobs in and out of the in-service
    phase multiset svcph; the returned flag says whether svcph changed so the
    caller forces a full propensity refresh (svcph is not in the stoichiometry,
    so the static dependency set misses it)."""
    ind = fromIR[kfire][0]
    svc_changed = False

    # buffered-PH departure drops the completing job's in-service phase (dep_phase); promotion refills the freed server at a fresh entry phase.
    if bufph_node is not None and bufph_node[ind] and is_buf_svc is not None and is_buf_svc[kfire]:
        r = fromIR[kfire][1]
        svcph[ind][r, int(dep_phase[kfire])] -= 1.0
        svc_changed = True

    # an order-independent departure is a pass-and-swap rewrite (chain shift + slot removal), not a promotion.
    if _is_list_sched(ind, sn) and buffers[ind]:
        buffers[ind] = _oi_depart(sn, ind, buffers[ind], fromIR[kfire][1])
        return svc_changed

    # departure from a buffered station promotes the discipline's chosen waiting job; a retrial station is the exception (orbit re-enters only via RETRY, never on departure).
    if (_is_buffered(ind, sn) and buffers[ind]
            and not _is_retrial_station(ind, sn)):
        pos = _pick_from_buffer(buffers[ind], sn, int(sn.nodeToStation[ind]))
        promoted = int(buffers[ind][pos])   # class id (1-based)
        del buffers[ind][pos]
        if bufph_node is not None and bufph_node[ind]:
            # The promoted waiting job starts service now, entering a phase drawn
            # from its entry distribution pie (the same allocation the init uses).
            ke = _draw_entry_phase(sn, ind, promoted - 1, int(smap.nph[ind, promoted - 1]))
            svcph[ind][promoted - 1, ke] += 1.0
            svc_changed = True

    # Destination is the (state-row) position selected by the routing draw
    if destPos < 0:
        return svc_changed
    arr_changed = _apply_arrival_buffer(int(smap.node[destPos]), int(smap.cls[destPos]),
                                        nvec, buffers, mi, R, sn, smap, svcph, bufph_node)
    return svc_changed or arr_changed


def _apply_arrival_buffer(jnd, r, nvec, buffers, mi, R, sn, smap: _Smap,
                          svcph=None, bufph_node=None):
    """Join a just-arrived class-r job to the ordered buffer of destination node
    jnd, if that node is buffered. nvec already includes the arrival. Shared by
    :func:`_update_buffers` (routed arrivals) and :func:`_fcr_release_cascade`
    (WAITQ releases), so the two paths cannot drift. Mirrors MATLAB
    applyArrivalBuffer. At a buffered-PH destination a job that enters service
    (rather than waiting) is added to svcph at a pie-drawn entry phase; the
    returned flag says whether svcph changed."""
    if _is_list_sched(jnd, sn):
        # PAS/OI arrival joins the back of the ordered list; capacity is the station's own cap, not a server/buffer split.
        if len(buffers[jnd]) < sn.cap[int(sn.nodeToStation[jnd])]:
            buffers[jnd].append(r + 1)      # newest last
        return False

    # Arrival at a buffered destination node
    if _is_buffered(jnd, sn):
        total_at_dest = float(np.sum(nvec[smap.node_slots(jnd)]))
        entered_service = False
        if _is_retrial_station(jnd, sn):
            # a retrial station breaks the usual buffer invariant (orbit may be occupied while servers idle), so an arrival must check the servers directly.
            in_svc = (total_at_dest - 1) - len(buffers[jnd])
            if in_svc >= mi[jnd]:
                buffers[jnd].appendleft(r + 1)
            else:
                entered_service = True
        elif total_at_dest > mi[jnd]:
            if _is_preemptive(jnd, sn):
                # preempt-resume arrival seizes a server; the displaced incumbent joins the buffer, leaving the new job in service.
                c = _pick_preempted(nvec, buffers[jnd], jnd, r, R, smap)
                if c > 0:
                    buffers[jnd].appendleft(c)   # addFirst
                entered_service = True
            else:
                # All servers busy - arriving job joins back of buffer
                buffers[jnd].appendleft(r + 1)   # addFirst
        else:
            # A server is free: the job goes straight into service.
            entered_service = True
        if entered_service and bufph_node is not None and bufph_node[jnd]:
            ke = _draw_entry_phase(sn, jnd, r, int(smap.nph[jnd, r]))
            svcph[jnd][r, ke] += 1.0
            return True
    return False


# =============================================================================
# Direct method: Gibson & Bruck NRM with indexed min-heap
# =============================================================================
def solver_ssa_nrm(sn, options: Optional[SolverSSAOptions] = None) -> SolverSSAReturn:
    """Steady-state analysis via the Next-Reaction Method with direct metrics."""
    if options is None:
        options = SolverSSAOptions()
    if options.seed and options.seed > 0:
        np.random.seed(options.seed)

    # immediate feedback (sn.immfeed) is not modeled by the NRM; self-loops are treated as ordinary re-queueing class-switching, so an immfeed model warns (use method='serial').
    _immfeed = getattr(sn, 'immfeed', None)
    if _immfeed is not None and np.any(_immfeed):
        import warnings
        warnings.warn("SolverSSA(method='nrm') does not model immediate feedback "
                      "(immfeed); self-loops are treated as class-switching with "
                      "re-queueing. Use method='serial' for immediate feedback.")

    samples = options.samples
    # mu/phi/pie/phasessz/phaseshift/phases and sn.proc rewrites are snapshotted and restored so nothing leaks into a later solver sharing this sn.
    import copy as _copy
    _orig = {f: _copy.deepcopy(getattr(sn, f))
             for f in ('mu', 'phi', 'pie', 'proc', 'phasessz', 'phaseshift',
                       'phases')}
    try:
        return _solver_ssa_nrm_run(sn, options, samples)
    finally:
        for f, v in _orig.items():
            setattr(sn, f, v)


def _solver_ssa_nrm_run(sn, options, samples) -> SolverSSAReturn:
    """Body of :func:`solver_ssa_nrm`, run with sn's phase fields populated."""
    # a model with Transition nodes is an SPN, not a queueing network, and is routed to the dedicated SPN runner (shares NRM clocks, own firing/vanishing-marking logic).
    for ind in range(sn.nnodes):
        if sn.nodetype[ind] == NodeType.TRANSITION:
            return _solver_ssa_nrm_spn(sn, options, samples)

    (S, D, a, fromIdx, fromIR, nvec0, buffers0, mi, fcr, balk, is_renege, sig,
     is_retry, rr, smap, n_dep_only, is_phase, is_poll, poll_info, is_poll_sw,
     poll_sw_node, svcph0, bufph_node, is_buf_svc, dep_phase, phase_from,
     phase_to, is_cache_rx, cache_hit_slot, cache_miss_slot, is_cache_node,
     is_cache_read, cache_retr_dest) = _build_nrm_problem(sn)

    M = sn.nstations
    K = sn.nclasses
    R = sn.nclasses
    NK = np.asarray(sn.njobs).flatten()
    _isslc = np.asarray(sn.isslc).flatten() if sn.isslc is not None else np.zeros(R)

    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))
    CN = np.zeros(K)
    XN = np.zeros(K)

    nvec = nvec0.copy()
    buffers = [deque(b) for b in buffers0]
    # Working copy of the in-service phase multiset (buffered-PH nodes only).
    svcph = [(arr.copy() if arr is not None else None) for arr in svcph0]
    # each polling controller is seeded to the first tangible state reachable from a canonical position; mirrors State.pollingInit / MATLAB NRM seed.
    poll_on = bool(np.any(is_poll))
    if poll_on:
        for ind in range(sn.nnodes):
            if not is_poll[ind]:
                continue
            pinfo = poll_info[ind]
            nbuf = _class_counts(nvec, smap, ind)
            q0, mode0, budget0 = polling_next(pinfo, 0, nbuf, R, arrived=True)
            buffers[ind] = _poll_land_ctrl(pinfo, q0, mode0, budget0)
    # Per-region WAITQ FIFO of parked (dst_node, dst_class) tokens, encoded as
    # dst_node*R + dst_class. Empty and untouched unless a region uses WAITQ.
    fcr_buf = [deque() for _ in range(fcr.F)] if fcr.on else []
    servers = np.asarray(sn.nservers).flatten()
    # Same normalized weights the propensities use, so the PS-family
    # utilization accumulator applies identical sharing factors.
    wnorm = _build_sched_weights(sn, R)
    # Per (cache node, read class) hit/miss event counts, used to reconstruct the
    # measured hit probability (State.afterEventCache convention).
    cache_prod = np.zeros((sn.nnodes, R))   # per (cache node, PRODUCED class) count

    numReactions = S.shape[1]
    Ak = np.array([a[i](nvec, buffers, svcph) for i in range(numReactions)])

    pq = IndexedMinHeap(numReactions)
    for kk in range(numReactions):
        pq.key[kk] = (-np.log(np.random.random()) / Ak[kk]) if Ak[kk] > 0.0 else float('inf')
    pq.build_heap()

    routing = _build_routing(S)
    jsq_flags = _build_jsq_flags(sn, fromIR, routing)
    sq_k = _build_sq_k(sn, fromIR, routing)

    sim_time = 0.0
    total_time = 0.0
    n = 0
    while n < samples:
        kfire = pq.peek_min()
        tau_fire = pq.peek_min_key()
        if np.isinf(tau_fire):
            break  # absorbing state

        dt = tau_fire - sim_time
        total_time += dt

        for ist in range(M):
            ind = int(sn.stationToNode[ist])
            sched = _sched_of(sn, ist)
            for kk in range(K):
                # nvec counts jobs per phase now, so the class population is the
                # sum over that class's phases.
                current_pop = _class_pop(nvec, smap, ind, kk)
                QN[ist, kk] += current_pop * dt

                # throughput sums the departure reactions of all a class's phases; phase-change and appended renege/retry reactions are excluded.
                dep_rate = 0.0
                for idx in range(n_dep_only):
                    if (fromIR[idx][0] == ind and fromIR[idx][1] == kk
                            and not is_phase[idx]):
                        dep_rate += Ak[idx]
                TN[ist, kk] += dep_rate * dt

                if sched in (SchedStrategy.INF, SchedStrategy.EXT):
                    UN[ist, kk] += current_pop * dt
                elif _is_ps_family(sched):
                    if _has_service(sn, ist, kk):
                        UN[ist, kk] += _ps_family_util(
                            sn, wnorm[ist], _class_counts(nvec, smap, ind), kk,
                            servers[ist], sched) * dt
                elif sched in _LIST_SCHED:
                    # PAS/OI utilization counts positions whose Delta_mu increment is positive, so a job served by several server types still counts once.
                    UN[ist, kk] += (_pas_in_svc(sn, ind, buffers[ind], kk)
                                    / servers[ist]) * dt
                elif sched == SchedStrategy.POLLING:
                    # polling class-kk utilization = fraction of time the controller is SERVING class kk.
                    ctrl = buffers[ind]
                    if len(ctrl) >= 2 and int(ctrl[0]) == 1 and int(ctrl[1]) == kk:
                        UN[ist, kk] += dt / servers[ist]
                elif sched in _BUFFERED_SCHED:
                    if _has_service(sn, ist, kk):
                        waiting = sum(1 for cl in buffers[ind] if cl == kk + 1)
                        in_service = current_pop - waiting
                        nsrv = servers[ist]
                        UN[ist, kk] += (in_service / nsrv) * dt

        sim_time = tau_fire
        cache_changed = False
        if is_cache_rx[kfire]:
            # cache access: read-class job draws an item, hit/miss decided against buffers[cacheNode] and replacement policy applied; mirrors State.afterEventCache.
            cn, rdc = fromIR[kfire]
            out_class, new_contents, cache_cat = _cache_access(sn, cn, rdc, buffers[cn])
            buffers[cn] = new_contents
            nvec[fromIdx[kfire]] -= 1.0            # consume the read-class job
            dest_pos = None
            if out_class >= 0:
                if cache_cat == 4:
                    # BEGIN retrieval routes the job to the fetch queue's destination rather than leaving it at the cache (else the cache-access reaction fires again).
                    dest_pos = int(cache_retr_dest[cn, out_class])
                else:
                    # hit or miss/completion produces the class at the SAME cache node; counted per produced class for hit/miss probability reconstruction.
                    dest_pos = int(smap.ph_off[cn, out_class])
                    cache_prod[cn, out_class] += 1.0
                nvec[dest_pos] += 1.0
            # out_class < 0 is a delayed hit: absorbed (produces nothing).
            cache_changed = True
        else:
            src_row = fromIdx[kfire]
            dest_pos = _fire_reaction(kfire, nvec, S, routing, src_row, jsq_flags, R, fcr, fromIR, sq_k, balk, sig, buffers, mi, rr, smap, fcr_buf, sn=sn)

            # a self-looping class re-enters the same node/class; at a buffered station point dest_pos at the source slot so _update_buffers rejoins the ordered buffer.
            if (dest_pos < 0 and not is_renege[kfire] and not is_retry[kfire]
                    and not is_poll_sw[kfire] and _isslc[fromIR[kfire][1]]):
                dest_pos = src_row

        svc_changed = False
        if is_retry[kfire]:
            # a successful retry moves an orbiting job to the free server; population is unchanged, only the orbit shrinks.
            _ind, _r = fromIR[kfire]
            for _b in range(len(buffers[_ind])):
                if buffers[_ind][_b] == _r + 1:
                    del buffers[_ind][_b]
                    break
        elif is_renege[kfire]:
            # reneging removes a waiting job (frees no server, promotes nothing); memoryless patience makes all waiting jobs of the class exchangeable.
            _ind, _r = fromIR[kfire]
            for _b in range(len(buffers[_ind])):
                if buffers[_ind][_b] == _r + 1:
                    del buffers[_ind][_b]
                    break
        elif is_phase[kfire] and bufph_node[fromIR[kfire][0]]:
            # a buffered-PH phase transition moves one in-service job between its own phases; frees no server, changes only svcph.
            _ind, _r = fromIR[kfire]
            svcph[_ind][_r, int(phase_from[kfire])] -= 1.0
            svcph[_ind][_r, int(phase_to[kfire])] += 1.0
            svc_changed = True
        elif is_cache_rx[kfire]:
            # The cache access already updated the cache contents and moved the job
            # to the hit/miss class above; there is no job buffer at a cache node.
            pass
        else:
            svc_changed = _update_buffers(kfire, nvec, buffers, fromIR, dest_pos, mi, R,
                                          sn, smap, svcph, bufph_node, is_buf_svc, dep_phase)

        # polling controller advance (DEP/SWITCH/parked-arrival) mirrors State.afterEventStation and forces a full propensity refresh.
        poll_changed = False
        if poll_on:
            src_node = fromIR[kfire][0]
            if is_poll_sw[kfire]:
                pind = int(poll_sw_node[kfire])
                pinfo = poll_info[pind]
                ctrl = buffers[pind]
                posS = int(ctrl[1]); swkS = int(ctrl[2])
                D0S = pinfo['sw_d0'][posS]
                KswS = int(pinfo['ksw'][posS])
                w = np.zeros(KswS + 1)
                for kd in range(KswS):
                    if (kd + 1) != swkS and D0S[swkS - 1, kd] > 0:
                        w[kd] = D0S[swkS - 1, kd]
                w[KswS] = max(0.0, -float(np.sum(D0S[swkS - 1, :])))  # absorption
                pick = _draw_from_dist(w)
                if pick < KswS and (pick + 1) != swkS:
                    buffers[pind] = [2, posS, pick + 1, 0]   # internal phase advance
                else:
                    nbufS = _class_counts(nvec, smap, pind)
                    qS, mdS, bgS = polling_next(pinfo, posS, nbufS, R, arrived=True)
                    buffers[pind] = _poll_land_ctrl(pinfo, qS, mdS, bgS)
                poll_changed = True
            elif is_poll[src_node] and kfire < n_dep_only and not is_phase[kfire]:
                pinfo = poll_info[src_node]
                ctrl = buffers[src_node]
                posD = int(ctrl[1]); ctrD = int(ctrl[3])
                nbufD = _class_counts(nvec, smap, src_node)   # after the departure
                ptype = pinfo['ptype']
                if ptype == PollingType.EXHAUSTIVE:
                    ctrnextD = 0;          goonD = nbufD[posD] > 0
                elif ptype == PollingType.GATED:
                    ctrnextD = ctrD - 1;   goonD = ctrnextD > 0
                elif ptype == PollingType.KLIMITED:
                    ctrnextD = ctrD - 1;   goonD = ctrnextD > 0 and nbufD[posD] > 0
                else:  # DECREMENTING
                    ctrnextD = ctrD;       goonD = nbufD[posD] > ctrD
                if goonD:
                    buffers[src_node] = [1, posD, 0, ctrnextD]
                else:
                    qD, mdD, bgD = polling_next(pinfo, posD, nbufD, R, arrived=False)
                    buffers[src_node] = _poll_land_ctrl(pinfo, qD, mdD, bgD)
                poll_changed = True
            if dest_pos is not None and dest_pos >= 0:
                jnd = int(smap.node[dest_pos])
                if is_poll[jnd]:
                    ctrlA = buffers[jnd]
                    if len(ctrlA) >= 1 and int(ctrlA[0]) == 0:   # parked server woken
                        pinfoA = poll_info[jnd]
                        nbufA = _class_counts(nvec, smap, jnd)    # includes the arrival
                        qA, mdA, bgA = polling_next(pinfoA, int(ctrlA[1]), nbufA, R, arrived=True)
                        buffers[jnd] = _poll_land_ctrl(pinfoA, qA, mdA, bgA)
                        poll_changed = True

        # WAITQ: admit parked jobs whose regions this firing may have relieved.
        # A release changes populations at arbitrary destination nodes.
        n_released = 0
        if fcr.on and fcr.any_waitq:
            n_released, rel_changed = _fcr_release_cascade(fcr, nvec, buffers, fcr_buf,
                                                           mi, R, sn, smap, svcph, bufph_node)
            svc_changed = svc_changed or rel_changed

        # a polling move, WAITQ release, or svcph update changes rates outside the fired reaction's static dependency set, so any of them forces a full refresh.
        refresh_set = range(numReactions) if (poll_changed or n_released > 0 or svc_changed or cache_changed) else D[kfire]
        for i in refresh_set:
            a_old = Ak[i]
            Ak[i] = a[i](nvec, buffers, svcph)
            a_new = Ak[i]
            if i == kfire:
                new_tau = (sim_time - np.log(np.random.random()) / a_new) if a_new > 0.0 else float('inf')
            elif a_new <= 0.0:
                new_tau = float('inf')
            elif a_old <= 0.0:
                # Reaction was dormant, now active - no residual time to rescale
                new_tau = sim_time - np.log(np.random.random()) / a_new
            else:
                new_tau = (a_old / a_new) * (pq.get_key(i) - sim_time) + sim_time
            pq.update(i, new_tau)

        n += 1

    if total_time > 0:
        QN /= total_time
        UN /= total_time
        TN /= total_time

    # class-dependent stations report utilization as T*S/peak against the declared per-class peak rate, matching the analytic solvers and serial SSA.
    cd_cell = getattr(sn, 'cdscaling', None)
    cd_peak = getattr(sn, 'cdscalingpeak', None)
    jd_cell = getattr(sn, 'jdscaling', None)
    jd_peak = getattr(sn, 'jdscalingpeak', None)
    if (cd_cell is not None and cd_peak is not None) or (jd_cell is not None and jd_peak is not None):
        for ist in range(M):
            has_cd = cd_cell is not None and ist < len(cd_cell) and cd_cell[ist] is not None
            has_jd = jd_cell is not None and ist < len(jd_cell) and jd_cell[ist] is not None
            if not (has_cd or has_jd):
                continue
            for kk in range(K):
                rate = sn.rates[ist, kk]
                # effective peak = product of declared cd and jd peaks (missing = 1)
                peak = 1.0
                if has_cd:
                    peak *= cd_peak[ist, kk]
                if has_jd:
                    peak *= jd_peak[ist, kk]
                if np.isfinite(rate) and rate > 0 and np.isfinite(peak) and peak > 0:
                    UN[ist, kk] = TN[ist, kk] / rate / peak
                else:
                    UN[ist, kk] = 0.0

    for kk in range(K):
        XN[kk] = TN[int(sn.refstat[kk]), kk]
        for ist in range(M):
            RN[ist, kk] = QN[ist, kk] / TN[ist, kk] if TN[ist, kk] > 0 else 0.0
        if XN[kk] > 0 and kk < len(NK) and np.isfinite(NK[kk]):
            CN[kk] = NK[kk] / XN[kk]

    for arr in (QN, UN, RN, XN, TN, CN):
        np.nan_to_num(arr, copy=False, nan=0.0)

    # measured hit/miss probabilities are written onto each Cache node's nodeparam for getAvgNode reconstruction; arrays are sized to nclasses to also cover internal retrieval classes.
    for ind in range(sn.nnodes):
        if is_cache_node[ind]:
            npm = sn.nodeparam[ind]
            ahp = np.full(R, np.nan)
            amp = np.full(R, np.nan)
            hitcl = np.atleast_1d(_np_get(npm, 'hitclass', [])).astype(int)
            misscl = np.atleast_1d(_np_get(npm, 'missclass', [])).astype(int)
            for r in range(R):
                if is_cache_read[ind, r] and r < len(hitcl) and hitcl[r] >= 0:
                    hc = int(hitcl[r]); mc = int(misscl[r])
                    tot = cache_prod[ind, hc] + cache_prod[ind, mc]
                    if tot > 0:
                        ahp[r] = cache_prod[ind, hc] / tot
                        amp[r] = cache_prod[ind, mc] / tot
            if isinstance(npm, dict):
                npm['actualhitprob'] = ahp
                npm['actualmissprob'] = amp
            else:
                npm.actualhitprob = ahp
                npm.actualmissprob = amp

    ret = SolverSSAReturn()
    ret.Q = QN
    ret.U = UN
    ret.R = RN
    ret.T = TN
    ret.C = CN.reshape(1, K)
    ret.X = XN.reshape(1, K)
    ret.A = _arvr_from_tput(sn, TN)
    ret.total_time = total_time
    ret.samples = n
    ret.method = 'nrm'
    return ret


# =============================================================================
# Stochastic Petri net (Place / Transition) via the Next-Reaction Method
#
# A Place holds a per-class token count (a population slot of the state vector);
# a timed Transition mode is a reaction whose stoichiometry column is the arc
# incidence -- input (enabling) arcs consume, output (firing) arcs produce.
# Enabling is a propensity gate (all input places at or above their arc weight,
# every inhibitor place strictly below its threshold); a single-server mode
# fires at its exponential rate, an infinite/k-server mode at that rate times
# its enabling degree. Each firing applies the mode's stoichiometry once, the
# atomic GSPN firing shared by the exact CTMC (single server), JMT and GreatSPN.
# IMMEDIATE modes fire in zero time and are resolved by vanishing-marking
# elimination (see _spn_collapse). A non-exponential firing distribution needs
# per-mode in-flight phase state the reaction network does not carry and is
# rejected here (the featset routes such nets to CTMC/JMT).
# =============================================================================
def _spn_is_immediate(t) -> bool:
    """True if a per-mode timing strategy denotes an immediate transition. The
    strategy may arrive as a TimingStrategy enum (value 1), or the legacy string
    'IMMEDIATE'."""
    if t is None:
        return False
    name = getattr(t, 'name', None)
    if name is not None:
        return name == 'IMMEDIATE'
    try:
        return int(t) == 1
    except (TypeError, ValueError):
        return str(t).upper() == 'IMMEDIATE'


class _SpnRx:
    """One transition-mode reaction of a stochastic Petri net."""
    __slots__ = ('node', 'mode', 'Svec', 'en_slot', 'en_w', 'inh_slot',
                 'inh_thr', 'base_rate', 'nservers', 'weight', 'prio')


def _spn_build_mode(sn, smap, ind, m, tp, NS):
    """Assemble the reaction record of transition ``ind`` mode ``m``.

    Enabling/firing/inhibiting are (nnodes, nclasses) matrices; a nonzero entry
    at (place, class) decodes to the slot ``ph_off[place, class]``.
    """
    R = sn.nclasses
    rec = _SpnRx()
    rec.node = ind
    rec.mode = m
    Svec = np.zeros(NS)
    en_slot = []
    en_w = []
    en = np.atleast_2d(np.asarray(tp.enabling[m], dtype=float))
    fir = np.atleast_2d(np.asarray(tp.firing[m], dtype=float))
    inh = np.atleast_2d(np.asarray(tp.inhibiting[m], dtype=float))
    for p in range(sn.nnodes):
        for c in range(R):
            w = en[p, c]
            if w > 0.0:
                slot = int(smap.ph_off[p, c])
                en_slot.append(slot)
                en_w.append(w)
                Svec[slot] -= w
            fw = fir[p, c]
            if fw > 0.0:
                Svec[int(smap.ph_off[p, c])] += fw
    inh_slot = []
    inh_thr = []
    for p in range(sn.nnodes):
        for c in range(R):
            thr = inh[p, c]
            if not np.isinf(thr):
                inh_slot.append(int(smap.ph_off[p, c]))
                inh_thr.append(thr)
    rec.Svec = Svec
    rec.en_slot = np.asarray(en_slot, dtype=int)
    rec.en_w = np.asarray(en_w, dtype=float)
    rec.inh_slot = np.asarray(inh_slot, dtype=int)
    rec.inh_thr = np.asarray(inh_thr, dtype=float)
    rec.base_rate = 0.0
    if not _spn_is_immediate(tp.timingstrategies[m]):
        fK = tp.firingphases[m] if m < len(np.atleast_1d(tp.firingphases)) else np.nan
        proc = tp.firingproc[m] if m < len(tp.firingproc) else None
        if np.isnan(fK) or fK != 1.0 or proc is None or len(proc) < 2:
            raise RuntimeError(
                "Transition node %d mode %d has non-exponential firing, which the "
                "NRM SPN path does not support." % (ind, m))
        D1 = np.atleast_2d(np.asarray(proc[1], dtype=float))
        rec.base_rate = float(np.sum(D1))
    ns = float(np.atleast_1d(tp.nmodeservers)[m])
    if np.isinf(ns):
        ns = float(GlobalConstants.MaxInt)
    rec.nservers = ns
    rec.weight = float(np.atleast_1d(tp.fireweight)[m])
    rec.prio = float(np.atleast_1d(tp.firingprio)[m])
    return rec


def _spn_en_degree(nvec, rx) -> float:
    """Enabling degree of a mode: min over input arcs of floor(tokens/weight),
    zeroed by any active inhibitor arc. A mode with no input arc is single."""
    for i in range(rx.inh_slot.shape[0]):
        if nvec[rx.inh_slot[i]] >= rx.inh_thr[i]:
            return 0.0
    if rx.en_slot.shape[0] == 0:
        return 1.0
    d = np.inf
    for i in range(rx.en_slot.shape[0]):
        d = min(d, np.floor(nvec[rx.en_slot[i]] / rx.en_w[i]))
    return d


def _spn_prop(nvec, rx) -> float:
    """Propensity of a timed mode: exponential rate times the effective server
    count, min(enabling degree, mode servers)."""
    d = _spn_en_degree(nvec, rx)
    eff = min(d, rx.nservers)
    if eff <= 0.0:
        return 0.0
    return rx.base_rate * eff


def _spn_weighted_draw(w) -> int:
    """Index drawn in proportion to the nonnegative weight vector ``w``."""
    tot = float(np.sum(w))
    if tot <= 0.0:
        return 0
    c = np.cumsum(w) / tot
    idx = np.where(c > np.random.random())[0]
    return int(idx[0]) if idx.size > 0 else (len(w) - 1)


def _spn_collapse(nvec, imm, maxsteps) -> None:
    """Vanishing-marking elimination. Fire enabled immediate transitions until
    the marking is tangible: highest firing priority first, ties resolved in
    proportion to firing weight. Immediate firings take zero time and advance no
    clock, so the timed race only ever samples from tangible markings."""
    if not imm:
        return
    steps = 0
    while True:
        enabled = [m for m in range(len(imm)) if _spn_en_degree(nvec, imm[m]) >= 1.0]
        if not enabled:
            return
        prios = np.array([imm[m].prio for m in enabled])
        maxp = prios.max()
        top = [enabled[i] for i in range(len(enabled)) if prios[i] == maxp]
        if len(top) == 1:
            pick = top[0]
        else:
            w = np.array([imm[m].weight for m in top])
            pick = top[_spn_weighted_draw(w)]
        nvec += imm[pick].Svec
        steps += 1
        if steps > maxsteps:
            raise RuntimeError(
                "Immediate-transition livelock: the vanishing-marking collapse "
                "did not reach a tangible marking.")


def _solver_ssa_nrm_spn(sn, options, samples) -> SolverSSAReturn:
    """Next-Reaction Method run for a stochastic Petri net."""
    R = sn.nclasses
    I = sn.nnodes
    M = sn.nstations
    K = sn.nclasses
    smap = _Smap(sn, I, R)
    NS = smap.NS

    rx: List[_SpnRx] = []      # timed modes
    imm: List[_SpnRx] = []     # immediate modes
    # consumers[(place,class)] indexes timed modes consuming from that place/class; Place throughput is their aggregate firing rate.
    consumers = {}
    for ind in range(I):
        if sn.nodetype[ind] != NodeType.TRANSITION:
            continue
        tp = sn.nodeparam[ind]
        nmodes = int(getattr(tp, 'nmodes', 0) if not isinstance(tp, dict)
                     else tp.get('nmodes', 0))
        _fmod = getattr(tp, 'firingdep', None) if not isinstance(tp, dict) else tp.get('firingdep', None)
        for m in range(nmodes):
            # Marking-dependent firing rates change the propensity with the marking;
            # the NRM SSA engine does not apply the g(marking) multiplier (unlike
            # CTMC and LDES), so reject rather than silently simulate the nominal rate.
            if _fmod is not None and m < len(_fmod) and _fmod[m] is not None:
                raise RuntimeError(
                    "SolverSSA does not support marking-dependent firing rates "
                    "(set_firing_rate_dependence); use SolverCTMC or SolverLDES.")
            rec = _spn_build_mode(sn, smap, ind, m, tp, NS)
            if _spn_is_immediate(tp.timingstrategies[m]):
                imm.append(rec)
            else:
                rx.append(rec)
                ridx = len(rx) - 1
                for a in range(rec.en_slot.shape[0]):
                    p = int(smap.node[rec.en_slot[a]])
                    c = int(smap.cls[rec.en_slot[a]])
                    consumers.setdefault((p, c), []).append(ridx)
    # Source arrivals get their own reaction (Poisson thinning per routed Place-class edge) since a Source is not a Transition; producers[(node,class)] indexes these for the Source's throughput report.
    from ....constants import ProcessType
    exp_id = int(ProcessType.EXP.value)
    rtnodes = np.asarray(sn.rtnodes, dtype=float)
    procid = (np.atleast_2d(np.asarray(sn.procid, dtype=object))
              if getattr(sn, 'procid', None) is not None else None)
    producers = {}
    for ind in range(I):
        if sn.nodetype[ind] != NodeType.SOURCE:
            continue
        ist = int(sn.nodeToStation[ind])
        for r in range(R):
            lam = sn.rates[ist, r]
            if np.isnan(lam) or lam <= 0.0:
                continue
            if procid is not None and ist < procid.shape[0] and r < procid.shape[1]:
                v = procid[ist, r]
                if v is not None and not (isinstance(v, float) and np.isnan(v)):
                    vid = int(v.value) if hasattr(v, 'value') else int(v)
                    if vid != exp_id:
                        raise RuntimeError(
                            "Source node %d class %d has a non-exponential arrival, "
                            "which the NRM SPN path does not support." % (ind, r))
            found_place = False
            for jnd in range(I):
                if sn.nodetype[jnd] != NodeType.PLACE:
                    continue
                for s in range(R):
                    p = rtnodes[ind * R + r, jnd * R + s]
                    if p <= 0.0:
                        continue
                    found_place = True
                    rec = _SpnRx()
                    rec.node = ind
                    rec.mode = -1   # arrival, not a transition mode
                    Svec = np.zeros(NS)
                    Svec[int(smap.ph_off[jnd, s])] += 1.0
                    rec.Svec = Svec
                    rec.en_slot = np.asarray([], dtype=int)
                    rec.en_w = np.asarray([], dtype=float)
                    rec.inh_slot = np.asarray([], dtype=int)
                    rec.inh_thr = np.asarray([], dtype=float)
                    rec.base_rate = float(lam) * float(p)   # Poisson thinning
                    rec.nservers = 1.0                      # constant propensity
                    rec.weight = 1.0
                    rec.prio = 1.0
                    rx.append(rec)
                    producers.setdefault((ind, r), []).append(len(rx) - 1)
            if not found_place:
                raise RuntimeError(
                    "Source node %d class %d does not route to any Place; the NRM "
                    "SPN path needs a Source->Place arc." % (ind, r))

    nR = len(rx)
    if nR == 0:
        raise RuntimeError(
            "Stochastic Petri net has no timed reaction; nothing to simulate.")

    # finite-capacity Place DROP: precompute per-slot/per-place caps and each reaction's deposited slots so the run-loop clamp only touches what just grew (JMT/CTMC loss semantics).
    pcap_slot = np.full(NS, np.inf)          # per-(place,class) slot cap
    place_total_caps = []                    # (total_cap, [slots]) per capped place
    for ind in range(I):
        if sn.nodetype[ind] != NodeType.PLACE or not sn.isstateful[ind]:
            continue
        ist = int(sn.nodeToStation[ind])
        slots_here = []
        for c in range(R):
            slot = int(smap.ph_off[ind, c])
            slots_here.append(slot)
            cc = sn.classcap[ist, c] if ist < np.asarray(sn.classcap).shape[0] else np.inf
            if np.isfinite(cc):
                pcap_slot[slot] = float(cc)
        tcap = sn.cap[ist] if ist < np.asarray(sn.cap).flatten().shape[0] else np.inf
        if np.isfinite(tcap):
            place_total_caps.append((float(tcap), slots_here))
    _has_place_caps = np.isfinite(pcap_slot).any() or len(place_total_caps) > 0
    dep_slots = [np.where(rec.Svec > 0.0)[0] for rec in rx]

    def _apply_place_caps(nv, deposited):
        # Drop tokens that a firing pushed above a place's per-class or total cap.
        for j in deposited:
            if nv[j] > pcap_slot[j]:
                nv[j] = pcap_slot[j]
        for tcap, slots in place_total_caps:
            excess = float(np.sum(nv[slots])) - tcap
            if excess > 0.0:
                for j in deposited:
                    if excess <= 0.0:
                        break
                    if j in slots and nv[j] > 0.0:
                        d = min(excess, nv[j])
                        nv[j] -= d
                        excess -= d

    # Initial marking: token counts per (place, class) from the initial state.
    nvec = np.zeros(NS)
    for ind in range(I):
        if sn.nodetype[ind] != NodeType.PLACE or not sn.isstateful[ind]:
            continue
        state_i = np.atleast_2d(np.asarray(sn.state[int(sn.nodeToStateful[ind])], dtype=float))
        _, nir = toMarginalAggr(sn, ind, state_i)
        nir = np.atleast_1d(np.asarray(nir).flatten())
        for c in range(R):
            v = nir[c] if c < nir.shape[0] else 0.0
            if np.isinf(v):
                raise RuntimeError("Infinite marking at a Place is not supported.")
            nvec[int(smap.ph_off[ind, c])] = v

    max_imm_steps = 100000  # livelock guard for the vanishing-marking collapse

    # ------------------------------------------------------------------
    # Next-Reaction Method run loop
    # ------------------------------------------------------------------
    _spn_collapse(nvec, imm, max_imm_steps)
    Ak = np.array([_spn_prop(nvec, rx[k]) for k in range(nR)])
    Pk = -np.log(np.random.random(nR))
    Tk = np.zeros(nR)
    with np.errstate(divide='ignore', invalid='ignore'):
        tau = (Pk - Tk) / Ak
    tau[Ak == 0.0] = np.inf

    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))
    CN = np.zeros(K)
    XN = np.zeros(K)
    NK = np.asarray(sn.njobs).flatten()
    total_time = 0.0

    n = 0
    while n < samples:
        kfire = int(np.argmin(tau))
        dt = tau[kfire]
        if np.isinf(dt):
            raise RuntimeError(
                "Deadlock: no transition is enabled. Quitting nrm method.")
        total_time += dt

        # a Place is an INF station: utilization is its mean token count, throughput is the summed firing rate of modes consuming from it.
        for ist in range(M):
            ind = int(sn.stationToNode[ist])
            for c in range(K):
                tokens = _class_pop(nvec, smap, ind, c)
                QN[ist, c] += tokens * dt
                UN[ist, c] += tokens * dt
                depr = 0.0
                for ridx in consumers.get((ind, c), ()):
                    depr += Ak[ridx]
                # A Source station has no consuming transition; its throughput is
                # the aggregate arrival rate it injects (producers).
                for ridx in producers.get((ind, c), ()):
                    depr += Ak[ridx]
                TN[ist, c] += depr * dt

        # fire the selected timed mode, then collapse any immediate transitions the new marking enabled, applying the DROP clamp before the immediate cascade sees it.
        nvec += rx[kfire].Svec
        if _has_place_caps:
            _apply_place_caps(nvec, dep_slots[kfire])
        _spn_collapse(nvec, imm, max_imm_steps)

        # advance Gibson & Bruck clocks with pre-firing propensities, then refresh every propensity (a firing plus its immediate cascade can touch any place).
        Tk = Tk + Ak * dt
        Ak = np.array([_spn_prop(nvec, rx[k]) for k in range(nR)])
        Pk[kfire] = Pk[kfire] - np.log(np.random.random())
        with np.errstate(divide='ignore', invalid='ignore'):
            tau = (Pk - Tk) / Ak
        tau[Ak == 0.0] = np.inf

        n += 1

    if total_time > 0:
        QN /= total_time
        UN /= total_time
        TN /= total_time
    for c in range(K):
        XN[c] = TN[int(sn.refstat[c]), c]
        for ist in range(M):
            RN[ist, c] = QN[ist, c] / TN[ist, c] if TN[ist, c] > 0 else 0.0
        if XN[c] > 0 and c < len(NK) and np.isfinite(NK[c]):
            CN[c] = NK[c] / XN[c]
    for arr in (QN, UN, RN, XN, TN, CN):
        np.nan_to_num(arr, copy=False, nan=0.0)

    ret = SolverSSAReturn()
    ret.Q = QN
    ret.U = UN
    ret.R = RN
    ret.T = TN
    ret.C = CN.reshape(1, K)
    ret.X = XN.reshape(1, K)
    ret.A = _arvr_from_tput(sn, TN)
    ret.total_time = total_time
    ret.samples = n
    ret.method = 'nrm'
    return ret


def _arvr_from_tput(sn, TN):
    """Best-effort arrival-rate matrix from throughput (mirrors basic SSA)."""
    try:
        from ...sn import sn_get_arvr_from_tput
        return sn_get_arvr_from_tput(sn, TN)
    except Exception:
        return None


# =============================================================================
# State-space method (mirror of Solver_ssa_nrm_space)
# =============================================================================
def _buffer_hash_all(bufs) -> str:
    return '|'.join('%d:[%s]' % (ind, ','.join(str(x) for x in bufs[ind]))
                    for ind in range(len(bufs)))


def _svcph_hash_all(svcph) -> str:
    if svcph is None:
        return ''
    return '|'.join('%d:%s' % (ind, np.array2string(np.asarray(arr).flatten(), precision=0))
                    for ind, arr in enumerate(svcph) if arr is not None)


def _hash_state(nvec, bufs, svcph=None) -> str:
    # svcph must enter the state cache key: two states with identical populations/buffers but different in-service phases have different propensities.
    return (np.array2string(np.asarray(nvec).flatten(), precision=12) + '|'
            + _buffer_hash_all(bufs) + '|' + _svcph_hash_all(svcph))


def solver_ssa_nrm_space(sn, options: Optional[SolverSSAOptions] = None):
    """
    State-space variant: returns (pi, outspace, depRates) keyed on the aggregate
    state and the FCFS/LCFS buffer contents. Mirrors Solver_ssa_nrm_space.
    """
    if options is None:
        options = SolverSSAOptions()
    if options.seed and options.seed > 0:
        np.random.seed(options.seed)

    samples = options.samples
    R = sn.nclasses
    I = sn.nnodes
    _isslc_sp = np.asarray(sn.isslc).flatten() if sn.isslc is not None else np.zeros(R)
    # See solver_ssa_nrm: _build_nrm_problem populates the phase fields from
    # sn.proc, so they are restored once the run has finished.
    import copy as _copy
    _orig = {f: _copy.deepcopy(getattr(sn, f))
             for f in ('mu', 'phi', 'pie', 'proc', 'phasessz', 'phaseshift',
                       'phases')}
    try:
        return _solver_ssa_nrm_space_run(sn, options, samples, R, I, _isslc_sp)
    finally:
        for f, v in _orig.items():
            setattr(sn, f, v)


def _solver_ssa_nrm_space_run(sn, options, samples, R, I, _isslc_sp):
    """Body of :func:`solver_ssa_nrm_space`, with sn's phase fields populated."""
    (S, D, a, fromIdx, fromIR, nvec0, buffers0, mi, fcr, balk, is_renege, sig,
     is_retry, rr, smap, n_dep_only, is_phase, _is_poll, _poll_info,
     _is_poll_sw, _poll_sw_node, svcph0, bufph_node, is_buf_svc, dep_phase,
     phase_from, phase_to) = _build_nrm_problem(sn)

    numReactions = S.shape[1]
    buffers = [deque(b) for b in buffers0]
    svcph = [(arr.copy() if arr is not None else None) for arr in svcph0]
    react_cache = {}

    Ak = np.array([a[i](nvec0, buffers, svcph) for i in range(numReactions)])
    nvec = nvec0.copy()
    react_cache[_hash_state(nvec, buffers, svcph)] = Ak.copy()
    Pk = -np.log(np.random.random(numReactions))
    Tk = np.zeros(numReactions)
    tau = (Pk - Tk) / np.where(Ak == 0, np.nan, Ak)
    tau[Ak == 0] = np.inf

    times = [0.0]
    states = []
    buffer_states = [[deque(b) for b in buffers]]

    routing = _build_routing(S)
    jsq_flags = _build_jsq_flags(sn, fromIR, routing)
    sq_k = _build_sq_k(sn, fromIR, routing)

    n = 0
    while n < samples:
        kfire = int(np.argmin(tau))
        dt = tau[kfire]
        if np.isinf(dt):
            raise RuntimeError("Deadlock. Quitting nrm method.")
        times.append(times[-1] + dt)

        src_row = fromIdx[kfire]
        dest_pos = _fire_reaction(kfire, nvec, S, routing, src_row, jsq_flags, R, fcr, fromIR, sq_k, balk, sig, buffers, mi, rr, smap, sn=sn)
        # Self-looping class re-enters the same node/class (see solver_ssa_nrm).
        if (dest_pos < 0 and not is_renege[kfire] and not is_retry[kfire]
                and _isslc_sp[fromIR[kfire][1]]):
            dest_pos = src_row
        svc_changed = False
        if is_retry[kfire]:
            # A successful retry moves one orbiting job into the free server;
            # the population is unchanged, so only the orbit shrinks.
            _ind, _r = fromIR[kfire]
            for _b in range(len(buffers[_ind])):
                if buffers[_ind][_b] == _r + 1:
                    del buffers[_ind][_b]
                    break
        elif is_renege[kfire]:
            # reneging removes a waiting job (frees no server, promotes nothing); memoryless patience makes all waiting jobs of the class exchangeable.
            _ind, _r = fromIR[kfire]
            for _b in range(len(buffers[_ind])):
                if buffers[_ind][_b] == _r + 1:
                    del buffers[_ind][_b]
                    break
        elif is_phase[kfire] and bufph_node[fromIR[kfire][0]]:
            # A buffered-PH phase transition moves one in-service job between its
            # phases; the buffer is untouched and only svcph changes.
            _ind, _r = fromIR[kfire]
            svcph[_ind][_r, int(phase_from[kfire])] -= 1.0
            svcph[_ind][_r, int(phase_to[kfire])] += 1.0
            svc_changed = True
        elif is_cache_rx[kfire]:
            # The cache access already updated the cache contents and moved the job
            # to the hit/miss class above; there is no job buffer at a cache node.
            pass
        else:
            svc_changed = _update_buffers(kfire, nvec, buffers, fromIR, dest_pos, mi, R,
                                          sn, smap, svcph, bufph_node, is_buf_svc, dep_phase)

        Tk += Ak * dt
        # an svcph move sits outside the static stoichiometry (like a polling move), so it forces a full propensity refresh.
        refresh_set = range(numReactions) if svc_changed else D[kfire]
        for i in refresh_set:
            Ak[i] = a[i](nvec, buffers, svcph)
        Pk[kfire] -= np.log(np.random.random())

        react_cache[_hash_state(nvec, buffers, svcph)] = Ak.copy()
        states.append(nvec.copy())
        buffer_states.append([deque(b) for b in buffers])

        with np.errstate(divide='ignore', invalid='ignore'):
            tau = (Pk - Tk) / Ak
        tau[Ak == 0] = np.inf
        n += 1

    # Empirical state probabilities over unique (nvec, buffers) states
    dt_arr = np.diff(np.asarray(times))
    num_intervals = len(states)  # == samples
    all_states = [nvec0] + states  # aligned with times / buffer_states

    keys = []
    for i in range(num_intervals):
        keys.append(_hash_state(all_states[i], buffer_states[i]))

    uniq_index = {}
    outspace_rows = []
    outspace_buffers = []
    ic = np.zeros(num_intervals, dtype=int)
    for i, key in enumerate(keys):
        if key not in uniq_index:
            uniq_index[key] = len(outspace_rows)
            outspace_rows.append(all_states[i])
            outspace_buffers.append(buffer_states[i])
        ic[i] = uniq_index[key]

    outspace = np.array(outspace_rows)
    time_accum = np.zeros(len(outspace_rows))
    for i in range(num_intervals):
        time_accum[ic[i]] += dt_arr[i]
    pi = time_accum / np.sum(time_accum)

    # depRates stays keyed by flat (node,class): phase-expanded departure reactions of one class are summed back into that column; phase transitions are excluded.
    num_states = outspace.shape[0]
    depRates = np.zeros((num_states, I * R))
    for st in range(num_states):
        a_state = react_cache.get(_hash_state(outspace[st], outspace_buffers[st]))
        if a_state is None:
            continue
        for j in range(n_dep_only):
            if is_phase[j]:
                continue
            ind, r = fromIR[j]
            depRates[st, ind * R + r] += a_state[j]

    return pi, outspace, depRates


__all__ = [
    'IndexedMinHeap',
    'solver_ssa_nrm',
    'solver_ssa_nrm_space',
]
