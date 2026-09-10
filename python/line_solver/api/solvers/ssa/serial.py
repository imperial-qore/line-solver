"""
afterEvent/sync-action SSA engine (Python-native).

Faithful port of the MATLAB stochastic-simulation core
``matlab/src/solvers/SSA/solver_ssa.m`` + ``solver_ssa_findenabled.m`` +
``solver_ssa_analyzer_serial.m`` (and the JAR ``Solver_ssa``/``Solver_ssa_findenabled``).

Unlike the legacy aggregate-count simulator, the engine drives the *same*
``State.afterEvent`` machinery used by the CTMC solver. It performs a lazy
random walk over the reachable state space: at the current per-stateful state it
enumerates the enabled sync actions (active ``afterEvent`` -> passive
``afterEvent`` -> routing probability), Gillespie-samples one micro-transition,
advances, and accumulates a time-weighted occupancy distribution plus per-state
arrival/departure rates. Each ``afterEvent`` output row is treated as a separate
candidate ``(next_state, rate)`` -- enumerate-all-rows + rate-proportional
selection is distributionally identical to MATLAB's ``isSimulation`` single-row
collapse, so no ``is_simulation`` flag is needed.

This module reuses (does not reimplement):
- ``refresh_sync`` / ``refresh_global_sync`` (``lang/sync.py``),
- ``after_event`` / ``after_global_event`` (``api/state``),
- ``toMarginal`` (``api/state/marginal.py``),
- ``_refresh_phase_fields`` (``api/solvers/ctmc/handler.py``) to populate the
  phase fields the afterEvent handlers require (getStruct does not).
"""

import math
import time
from typing import Optional

import numpy as np

from ....constants import EventType, SchedStrategy
from ....lang.base import NodeType
from ....lang.sync import refresh_sync, refresh_global_sync
from ...state.after_event import after_event, after_event_init
from ...state.after_fj_event import after_fj_event
from ...state.after_global_event import after_global_event
from ...state.marginal import toMarginal, fromMarginal, _ext_class_disabled
from ...state.ctmc_ssg import _is_bas_station

# Entries the per-run after_event memo may hold before it stops growing. The key
# is one NODE state, not a global one, so a few thousand covers most models;
# the cap only bounds the pathological case.
_AE_MEMO_MAX = 50000
from ...sn import sn_region_members


def _gd_factor_now(sn, cur_state):
    """Evaluate the global (Whittle) rate scaling phi(n) on the CURRENT state.

    Returns the (nstations, nclasses) matrix of scalings. SolverCTMC tabulates
    phi once per state of the enumerated space; a simulator has one state at a
    time, so the same factorization applies with the table collapsed to a single
    row. Within a state phi is a constant multiplying every station service rate.
    """
    M = int(sn.nstations)
    K = int(sn.nclasses)
    npop = np.zeros((M, K))
    for ind in range(int(sn.nnodes)):
        if not sn.isstateful[ind] or not sn.isstation[ind]:
            continue
        isf = int(sn.nodeToStateful[ind])
        ist = int(sn.nodeToStation[ind])
        nir = np.asarray(toMarginal(sn, ind, np.atleast_2d(cur_state[isf]))[1],
                         dtype=float).ravel()[:K]
        npop[ist, :len(nir)] = nir
    v = np.asarray(sn.gdscaling(npop), dtype=float)
    if v.size == 1:
        out = np.full((M, K), float(v.ravel()[0]))
    elif v.shape == (M,) or v.shape == (M, 1):
        out = np.tile(v.reshape(M, 1), (1, K))
    else:
        out = v.reshape(M, K)
    if not np.all(np.isfinite(out)) or np.any(out < 0):
        raise ValueError("The global dependence handle returned a non-finite or "
                         "negative scaling.")
    return out


def _sched_name(sn, ist):
    sched = sn.sched.get(ist) if isinstance(sn.sched, dict) else sn.sched[ist]
    return sched.name if hasattr(sched, 'name') else str(sched)


def _recompute_capacity(sn, cutoff):
    """Bound open-class populations by the cutoff (mirrors solver_ssa.m:59-93).

    Sets ``sn.classcap`` (M x R), ``sn.cap`` (M,) and fills infinite
    ``sn.nservers`` so the afterEvent ARV handlers block arrivals consistently
    with the CTMC truncation. Pass-and-swap stations keep their total
    order-independent buffer capacity (not the per-class sum).

    Returns the (M x R) boolean predicate ``is_physical_cap``. A refusal at
    (ist,r) is PHYSICAL when the buffer that binds it was declared by the user,
    and an ARTIFACT when the binding bound is the state-space cutoff that this
    routine folds into sn.classcap/sn.cap. The distinction is a throughput
    invariant: a job refused by a real buffer still departed the node it left and
    is offered load, whereas a job refused by the cutoff never existed.
    """
    M, R = sn.nstations, sn.nclasses
    nj = np.asarray(sn.njobs, dtype=float).reshape(-1)
    classcap = np.array(sn.classcap, dtype=float).reshape(M, R)
    cap = np.array(sn.cap, dtype=float).reshape(-1).astype(float).copy()
    nservers = np.array(sn.nservers, dtype=float).reshape(-1).astype(float)
    chains = np.asarray(sn.chains) if getattr(sn, 'chains', None) is not None else None
    fjclassmap = np.asarray(sn.fjclassmap).ravel() if getattr(sn, 'fjclassmap', None) is not None \
        and np.size(sn.fjclassmap) else None
    capacityc = np.zeros((M, R))
    is_physical_cap = np.ones((M, R), dtype=bool)
    for ist in range(M):
        for r in range(R):
            if fjclassmap is not None and r < len(fjclassmap) and fjclassmap[r] >= 0:
                # FJ auxiliary sibling class: adapter-set classcap is the
                # capacity bound (tasksPerLink can exceed the chain population).
                capacityc[ist, r] = classcap[ist, r]
            elif np.isinf(nj[r]):
                capacityc[ist, r] = min(cutoff, classcap[ist, r])
                # only open classes can be bound by the cutoff; closed classes take
                # min(chainpop, classcap) below with no cutoff involved
                if cutoff < classcap[ist, r]:
                    is_physical_cap[ist, r] = False
            elif chains is not None:
                c = None
                for ci in range(chains.shape[0]):
                    if r < chains.shape[1] and chains[ci, r] > 0:
                        c = ci
                        break
                pop = 0.0
                if c is not None:
                    for kk in range(R):
                        if kk < chains.shape[1] and chains[c, kk] > 0 and np.isfinite(nj[kk]):
                            pop += nj[kk]
                # closed classes: enumerate up to the chain population, but never
                # beyond the class capacity at this station (finite-buffer stations)
                capacityc[ist, r] = min(pop, classcap[ist, r])
            else:
                capacityc[ist, r] = min(nj[r], classcap[ist, r])
        # never raise the station capacity above its configured total capacity
        cap_sum = float(min(np.sum(capacityc[ist, :]), cap[ist]))
        # Both conjuncts needed: without the truncation check, a closed class
        # whose declared classcap exceeds its chain population also leaves
        # sum(capacityc) < cap with no cutoff involved, misreporting capacity.
        if not np.all(is_physical_cap[ist, :]) and np.sum(capacityc[ist, :]) < cap[ist]:
            is_physical_cap[ist, :] = False
        if _sched_name(sn, ist) == 'PAS' and np.isfinite(cap[ist]):
            cap_sum = float(cap[ist])   # PAS: single OI buffer of total size cap
        if _sched_name(sn, ist) == 'EXT':
            # Source: keep a finite server count so int(nservers) stays well-defined.
            nservers[ist] = float(max(1, R))
            cap[ist] = cap_sum
            continue
        if np.isinf(nservers[ist]):
            nservers[ist] = cap_sum
        cap[ist] = cap_sum
    sn.classcap = capacityc
    sn.cap = cap
    sn.nservers = nservers
    return is_physical_cap


def _state_nir_vector(sn, cur_state, M, R):
    """Per-station per-class marginal counts as a length-(M*R) vector."""
    vec = np.zeros(M * R)
    for ist in range(M):
        ind = int(sn.stationToNode[ist])
        isf = int(sn.nodeToStateful[ind])
        _ni, nir, _sir, _kir = toMarginal(sn, ind, cur_state[isf])
        vec[ist * R:(ist + 1) * R] = np.ravel(nir)[:R]
    return vec


def _state_key(cur_state):
    """Hashable identity of the full per-stateful state (variable widths ok)."""
    return tuple(tuple(np.ravel(s).tolist()) for s in cur_state)


def _eval_state_dep_prob(pprob, cur_state, isf_a, row_a, isf_p, row_p, nstateful):
    """Evaluate a state-dependent routing probability (mirrors solver_ssa.m:208)."""
    before = [np.atleast_2d(cur_state[i]) for i in range(nstateful)]
    after = [b.copy() for b in before]
    after[isf_a] = np.atleast_2d(row_a)
    if isf_p >= 0:
        after[isf_p] = np.atleast_2d(row_p)
    val = pprob(before, after)
    if isinstance(val, np.ndarray):
        val = float(val.ravel()[0]) if val.size else 0.0
    return float(val)


def solver_ssa_serial(sn, options, seed=None):
    """Run the afterEvent/sync-action Gillespie simulation.

    Returns ``(pi, SSq, arvRates, depRates, total_time, event_log, state_log)``
    where ``pi`` is the normalized time-occupancy over unique visited states,
    ``SSq`` the (nuniq x M*R) marginal-count matrix, and ``arvRates``/``depRates``
    the per-state (nuniq x nstateful x R) rate tensors.

    When ``options.record_events`` is set, ``event_log`` is the chronological
    per-transition trace ``[(t, src_st, src_k, dst_st, dst_k), ...]`` consumed by
    ``SolverSSA.sample`` (station indices; ``-1`` where there is no departing /
    arriving station), and ``state_log`` the matching pre-fire (M x R) marginal
    counts. Both are ``None`` otherwise.
    """
    if seed is None:
        seed = int(getattr(options, 'seed', 0) or 0)
    if seed and seed > 0:
        np.random.seed(int(seed))

    cutoff = getattr(options, 'cutoff', np.inf)
    if isinstance(cutoff, np.ndarray):
        cutoff = float(np.max(cutoff)) if cutoff.size else np.inf
    cutoff = float(cutoff) if cutoff is not None else np.inf
    # cutoff == inf: buffers grow dynamically (is_simulation=True), no truncation.
    # A finite cutoff caps the per-class population (matching CTMC).
    # Populate phase fields (mu/phi/pie/phasessz/phaseshift) the handlers need.
    from ..ctmc.handler import _refresh_phase_fields
    _refresh_phase_fields(sn, immediate_as_rate=True)
    is_physical_cap = _recompute_capacity(sn, cutoff)
    # Loop-invariant after_event context; built after the refreshes above so it
    # reflects the same sn the Gillespie loop passes to after_event.
    aectx = after_event_init(sn)

    # The enabled-transition enumeration is a PURE function of (node, node
    # state, event, class): after_event reads sn and its arguments, never
    # sn.state, and this loop never writes sn.state back. A sample path revisits
    # the same node state constantly -- every step re-enumerates every sync
    # action -- so the answer is memoised for the run, the same reaction cache
    # the NRM engine keeps (nrm.py react_cache, JAR Solver_ssa_nrm_space).
    # On the delayed-hit retrieval fixture of test_ssa_nrm_cache that is a 98%
    # hit rate and an 8x speedup with a bit-identical sample path.
    # Callers only READ the stored arrays -- they copy through
    # atleast_2d/ravel/astype before use -- so every hit shares one tuple.
    # The cap bounds memory on models whose per-node state space is large; past
    # it the engine simply stops memoising and stays correct.
    _ae_memo = {}

    def _after_event_cached(ind, inspace, event, job_class, no_promote=False):
        row = np.ascontiguousarray(inspace, dtype=float)
        key = (ind, event, int(job_class), bool(no_promote), row.shape, row.tobytes())
        out = _ae_memo.get(key)
        if out is None:
            out = after_event(sn, ind, inspace, event, job_class, True,
                              ctx=aectx, no_promote=no_promote)
            if len(_ae_memo) < _AE_MEMO_MAX:
                _ae_memo[key] = out
        return out
    # Global (Whittle) rate scaling declared through set_global_dependence.
    _has_gd = getattr(sn, 'gdscaling', None) is not None

    nstateful = sn.nstateful
    M, R = sn.nstations, sn.nclasses
    nnodes = sn.nnodes

    # see _kb/06-solver-catalog.md (SSA: "Python serial engine: state-
    # construction traps") -- nvars must be populated before refresh_sync
    from ...state.ctmc_ssg import _space_local_vars
    from ...state.polling import (polling_width as _polling_width,
                                  polling_init as _polling_init,
                                  _srvclass_of)
    if sn.nvars is None:
        sn.nvars = np.zeros((nnodes, 1), dtype=int)
        # BUG-83: sn.isbasblocking set for the blocking station under both
        # upstream and destination declaration forms; rejects polling+BAS collision.
        from ...state.ctmc_ssg import _populate_isbasblocking
        _populate_isbasblocking(sn)
    _local_var_seed = {}
    _init_nir = {}
    _bas_marker = set()   # stateful indices whose row carries a trailing BAS marker
    for _ind in range(nnodes):
        if not sn.isstateful[_ind]:
            continue
        _isf = int(sn.nodeToStateful[_ind])
        # see _kb/06-solver-catalog.md (SSA: "Python serial engine: state-
        # construction traps") -- use the compact marginal as-is, not via toMarginal
        if sn.isstation[_ind] and sn.state is not None:
            _st = np.ravel(np.atleast_2d(sn.state[_isf])[0]).astype(float)
            if _st.size == R:
                _init_nir[_isf] = _st[:R]
            else:
                _init_nir[_isf] = np.ravel(
                    toMarginal(sn, _ind, np.atleast_2d(sn.state[_isf])[0])[1])[:R]
        _lv = _space_local_vars(sn, _ind, first_only=True)
        if _lv is not None and np.atleast_2d(_lv).size > 0:
            _lv = np.atleast_2d(_lv)
            sn.nvars[_ind, 0] = int(_lv.shape[1])
            _local_var_seed[_isf] = _lv[0].astype(float)
        # see _kb/06-solver-catalog.md (SSA: "Python serial engine: state-
        # construction traps") -- True BAS marker column, seeded to 0 (unblocked)
        if sn.isstation[_ind] and _is_bas_station(sn, int(sn.nodeToStation[_ind])):
            sn.nvars[_ind, 0] = int(sn.nvars[_ind, 0]) + 1
            _bas_marker.add(_isf)
        # Polling controller [pos, swk, ctr]: width only reserved here; value seeded below.
        _pw = _polling_width(sn, _ind)
        if _pw > 0:
            sn.nvars[_ind, 0] = int(sn.nvars[_ind, 0]) + _pw

    sync = refresh_sync(sn)
    if getattr(sn, 'gsync', None) is None or len(sn.gsync) == 0:
        sn.gsync = refresh_global_sync(sn)
    gsync = sn.gsync if sn.gsync is not None else []
    # Iterate refresh_global_sync's VALUES, not the dict itself (keys lack .active).
    gsync_events = list(gsync.values()) if isinstance(gsync, dict) else list(gsync)
    # fork firing synchronizations (native fork-join)
    fjsync = sn.fjsync if getattr(sn, 'fjsync', None) else []
    FJ = len(fjsync)

    acts = [(a.active.node, a.active.event, a.active.job_class,
             a.passive.node, a.passive.event, a.passive.job_class, a.passive.prob)
            for a in sync]

    # Station states use full capacity width (buffer sized to cap, zero-padded),
    # matching the CTMC state encoding, so afterEvent handlers need no buffer growth.
    phasessz = np.asarray(sn.phasessz, dtype=int).reshape(sn.nstations, R)
    # phaseshift: native is (nstations,R); MATLAB is (nstations,R+1) ([0,cumsum]);
    # first R columns are identical in both.
    phaseshift = np.asarray(sn.phaseshift, dtype=int).reshape(sn.nstations, -1)[:, :R]
    # sn.nvars and the per-node local-variable seed rows were computed above,
    # before refresh_sync; reuse them here.
    nvars = sn.nvars
    cur_state = [None] * nstateful
    for ind in range(nnodes):
        if not sn.isstateful[ind]:
            continue
        isf = int(sn.nodeToStateful[ind])
        if not sn.isstation[ind] and sn.nodetype[ind] == NodeType.TRANSITION:
            # Transition node initial state [idle,phases,fired,vars], mirroring
            # ctmc_ssg.py so State.afterGlobalEvent can split the row.
            nparam = sn.nodeparam[ind]
            # Marking-dependent firing rates change the propensity with the marking;
            # the serial SSA engine does not apply the g(marking) multiplier (unlike
            # CTMC and LDES), so reject rather than silently simulate the nominal rate.
            _fmod = getattr(nparam, 'firingdep', None) if not isinstance(nparam, dict) else nparam.get('firingdep', None)
            if _fmod is not None and any(g is not None for g in _fmod):
                raise RuntimeError(
                    "SolverSSA does not support marking-dependent firing rates "
                    "(set_firing_rate_dependence); use SolverCTMC or SolverLDES.")
            nmodes_t = int(getattr(nparam, 'nmodes', 0) if not isinstance(nparam, dict)
                           else nparam.get('nmodes', 0))
            fphases = np.asarray(getattr(nparam, 'firingphases', np.array([])), dtype=float)
            fp_t = getattr(nparam, 'firingproc', None)
            fK_t = np.ones(nmodes_t, dtype=int)
            for _m in range(nmodes_t):
                if _m < fphases.size and not np.isnan(fphases[_m]):
                    fK_t[_m] = int(fphases[_m])
                elif fp_t is not None and _m < len(fp_t) and fp_t[_m] is not None:
                    fK_t[_m] = int(np.atleast_2d(np.asarray(fp_t[_m][0])).shape[0])
            fK_t = np.maximum(fK_t, 1)
            from ....constants import GlobalConstants as _GC
            _maxint = int(_GC.MaxInt) if hasattr(_GC, 'MaxInt') else (2 ** 31 - 1)
            nms = np.asarray(getattr(nparam, 'nmodeservers', np.ones(nmodes_t)), dtype=float)
            if nms.size != nmodes_t:
                nms = np.broadcast_to(nms, (nmodes_t,)).astype(float).copy()
            idle0 = np.array([_maxint if not np.isfinite(nms[_m]) else int(nms[_m])
                              for _m in range(nmodes_t)], dtype=float)
            V_t = int(np.sum(nvars[ind])) if nvars is not None else 0
            cur_state[isf] = np.concatenate([
                idle0,
                np.zeros(int(np.sum(fK_t))),
                np.zeros(nmodes_t),
                np.zeros(V_t),
            ])
            continue
        if not sn.isstation[ind]:
            # Non-station node (Cache/Router): a cache is never empty, so seed
            # with the first enumerated var configuration (mirrors CTMC init-state
            # selection). Round-robin routers get their first outlink pointer.
            srv0 = np.atleast_2d(sn.state[isf])[0].astype(float).copy()
            lv0 = _local_var_seed.get(isf)
            if lv0 is not None and lv0.size > 0:
                cur_state[isf] = np.concatenate([srv0, lv0])
            else:
                cur_state[isf] = srv0
            continue
        ist = int(sn.nodeToStation[ind])
        sumK = int(np.sum(phasessz[ist, :]))
        V = int(np.sum(nvars[ind])) if nvars is not None else 0
        if _sched_name(sn, ist) == 'EXT':
            # Started Source: each enabled class's renewal active in entry phase.
            # phasessz is the column WIDTH (always >=1) and never the enabled
            # test; a class that does not arrive here keeps a zero column, as in
            # MATLAB State.fromMarginal's EXT branch.
            row = np.zeros(sumK + V)
            for r in range(R):
                if (phasessz[ist, r] > 0 and int(phaseshift[ist, r]) < row.size
                        and not _ext_class_disabled(sn, ist, r)):
                    row[int(phaseshift[ist, r])] = 1.0
            # Seed RR outlink pointer in trailing V columns; a zero pointer
            # deadlocks the simulation (every RR routing probability would be 0).
            lv0 = _local_var_seed.get(isf)
            if lv0 is not None and lv0.size > 0:
                row[sumK:sumK + lv0.size] = lv0
            cur_state[isf] = row
            continue
        if _sched_name(sn, ist) == 'PAS':
            # PAS: ordered-list state of fixed width = total OI buffer capacity.
            captot = sn.cap[ist]
            cur_state[isf] = np.zeros(int(captot) if np.isfinite(captot) else 0)
            if isf in _bas_marker:
                cur_state[isf] = np.concatenate([cur_state[isf], np.zeros(1)])
            continue
        # Minimal initial state for the initial marginal; the afterEvent ARV
        # handler grows the buffer on demand, so no capacity-width pre-sizing.
        nir0 = _init_nir.get(isf)
        if nir0 is None:
            nir0 = np.ravel(toMarginal(sn, ind, np.atleast_2d(sn.state[isf])[0])[1])[:R]
        srv0 = np.atleast_2d(fromMarginal(sn, ind, nir0))[0].astype(float)
        # RR/WRR outlink pointer appended after the server portion, matching
        # the CTMC state encoding, so afterEvent DEP can read/advance it.
        lv0 = _local_var_seed.get(isf)
        if lv0 is not None and lv0.size > 0:
            cur_state[isf] = np.concatenate([srv0, lv0])
        else:
            cur_state[isf] = srv0
        # BAS blocked marker: the station starts unblocked. Trails the cache/round-robin
        # vars and precedes the polling controller, matching the CTMC column order.
        if isf in _bas_marker:
            cur_state[isf] = np.concatenate([cur_state[isf], np.zeros(1)])
        # Polling controller seeded from the chosen row (api/state/polling.polling_init).
        if _polling_width(sn, ind) > 0:
            _Kp = np.array(phasessz[ist, :], dtype=int)
            _Ksp = np.array(phaseshift[ist, :], dtype=int)
            _nbuf0 = srv0[:R].copy()
            _srvcls0 = _srvclass_of(srv0[R:R + sumK], _Kp, _Ksp, R)
            cur_state[isf] = np.concatenate([
                cur_state[isf], _polling_init(sn, ind, _nbuf0, _srvcls0)])

    samples = int(getattr(options, 'samples', 10000))
    timespan = getattr(options, 'timespan', None)
    t_end = float(timespan[1]) if (timespan is not None and len(timespan) > 1) else np.inf
    # Wall-clock budget (options.timeout); loop breaks and yields interim
    # statistics once exceeded.
    _tmo = float(getattr(options, 'timeout', float('inf')))
    _tmo_start = time.time()

    # see _kb/06-solver-catalog.md (SSA: "Python serial engine: state-
    # construction traps") -- ssa_max_pop safeguard against unbounded buffers
    _max_pop = float(getattr(options, 'ssa_max_pop', 1e6))
    if not math.isfinite(_max_pop) or _max_pop <= 0:
        _max_pop = float('inf')
    _POP_CHECK_EVERY = 256   # amortize the toMarginal cost
    # Source stations carry an (intentionally) infinite open-class population;
    # the guard must bound only the queueable in-system population, so mask them.
    _pop_st_mask = np.ones(M, dtype=bool)
    for _ist in range(M):
        _nd = int(sn.stationToNode[_ist])
        if _nd < len(sn.nodetype) and sn.nodetype[_nd] == NodeType.SOURCE:
            _pop_st_mask[_ist] = False

    tally = {}            # key -> [dwell, nir_vec, depRates, arvRates, dlyRates, visits, startRates, preemptRates]
    total_time = 0.0
    cur_time = 0.0

    # see _kb/09-ldes-and-cache.md (SSA: "delayed-hit rate from the merge transition").
    # Cache row layout is [srv(R) | var(V)], so the srv block ends V columns
    # from the right; keep the slice per node rather than assuming V == 0.
    cache_srv_slice = {}
    for _ind in range(nnodes):
        if _ind < len(sn.nodetype) and sn.nodetype[_ind] == NodeType.CACHE:
            _csf = int(sn.nodeToStateful[_ind])
            if _csf >= 0:
                _cv = int(np.sum(sn.nvars[_ind])) if sn.nvars is not None else 0
                cache_srv_slice[_csf] = _cv

    # Per-transition event log for SolverSSA.sample; only built when requested
    # (tally-based analysis needs no extra bookkeeping).
    record_events = bool(getattr(options, 'record_events', False))
    event_log = [] if record_events else None
    state_log = [] if record_events else None
    # (time, station, class, kind) rows with kind 0 = START and 1 = PREEMPT
    tag_log = [] if record_events else None
    sf2st = np.asarray(sn.statefulToStation, dtype=int).reshape(-1)

    def _sf_to_station(isf):
        return int(sf2st[isf]) if 0 <= isf < sf2st.size else -1

    # see _kb/06-solver-catalog.md (SSA: "Serial engine: finite capacity
    # regions (FCR), fork-join, immediate feedback")
    _fcr_on = bool(getattr(sn, 'nregions', 0) and int(sn.nregions) > 0 and getattr(sn, 'region', None))
    _fcr = []
    _fcr_buf = []
    if _fcr_on:
        from ....lang.base import DropStrategy as _BaseDrop
        _rr = (np.atleast_2d(np.asarray(sn.regionrule, dtype=float))
               if getattr(sn, 'regionrule', None) is not None and np.size(sn.regionrule) > 0
               else None)
        for f in range(int(sn.nregions)):
            Rmat = np.asarray(sn.region[f], dtype=float)
            Mrows = Rmat.shape[0]
            memvec = -np.ones(Mrows)
            if (getattr(sn, 'regionmaxmem', None) and len(sn.regionmaxmem) > f
                    and sn.regionmaxmem[f] is not None and np.size(sn.regionmaxmem[f]) > 0):
                mv_ = np.asarray(sn.regionmaxmem[f], dtype=float).ravel()[:Mrows]
                memvec[:mv_.size] = mv_
            member_mask = sn_region_members(sn, f, Rmat, memvec)
            members = [i for i in range(Mrows) if member_mask[i]]
            classcap = np.full(R, np.inf)
            for r in range(R):
                cv = [Rmat[i, r] for i in members if Rmat[i, r] != -1]
                if cv:
                    classcap[r] = min(cv)
            gv = [Rmat[i, R] for i in members if Rmat[i, R] != -1]
            globalcap = min(gv) if gv else np.inf
            memcap = np.inf
            mvals = [memvec[i] for i in members if memvec[i] != -1]
            if mvals:
                memcap = min(mvals)
            szrow = (np.asarray(sn.regionsz, dtype=float)[f].ravel()
                     if getattr(sn, 'regionsz', None) is not None and np.size(sn.regionsz) > 0
                     else np.ones(R))
            linA = None
            linb = None
            if (getattr(sn, 'regionlincon', None) and len(sn.regionlincon) > f
                    and sn.regionlincon[f] is not None):
                linA = np.atleast_2d(np.asarray(sn.regionlincon[f][0], dtype=float))
                linb = np.asarray(sn.regionlincon[f][1], dtype=float).ravel()
            iswaitq = np.zeros(R, dtype=bool)
            for r in range(R):
                if _rr is not None and f < _rr.shape[0] and r < _rr.shape[1]:
                    iswaitq[r] = (_rr[f, r] != int(_BaseDrop.DROP))
            _fcr.append((members, member_mask, classcap, globalcap, memcap, szrow, linA, linb, iswaitq))
            _fcr_buf.append([])

    def _fcr_violates(fi, xn):
        (members, member_mask, classcap, globalcap, memcap, szrow, linA, linb, iswaitq) = _fcr[fi]
        if (np.any(xn > classcap) or xn.sum() > globalcap
                or (np.isfinite(memcap) and float(xn.dot(szrow[:R])) > memcap)):
            return True
        if linA is not None:
            return bool(np.any(linA @ xn.reshape(-1, 1) > linb.reshape(-1, 1)))
        return False

    def _fcr_regionpop(fi, state_cells):
        members = _fcr[fi][0]
        xf = np.zeros(R)
        for i in members:
            ind_i = int(sn.stationToNode[i])
            isf_i = int(sn.stationToStateful[i])
            xf += np.ravel(toMarginal(sn, ind_i, state_cells[isf_i])[1])[:R]
        return xf

    def _fcr_release(state_cells):
        # strict-FIFO head-of-line release of parked tokens
        progress = True
        while progress:
            progress = False
            for fi in range(len(_fcr_buf)):
                if not _fcr_buf[fi]:
                    continue
                x = _fcr_regionpop(fi, state_cells)
                tok = _fcr_buf[fi][0]
                dest = tok // R
                r = tok % R
                xn = x.copy()
                xn[r] += 1
                if _fcr_violates(fi, xn):
                    continue
                isf_d = int(sn.nodeToStateful[dest])
                od, _rd, _pd, _sd, _ppd = _after_event_cached(dest, state_cells[isf_d], EventType.ARV, r)
                od = np.atleast_2d(od)
                if od.size == 0:
                    continue
                state_cells[isf_d] = od[0].astype(float)
                _fcr_buf[fi].pop(0)
                progress = True
        return state_cells

    # see _kb/06-solver-catalog.md (SSA: "Warmup discard, and cache hit/miss accounting")
    _wf = float(getattr(options, 'warmupfrac', 0.0) or 0.0)
    n_drop = int(math.floor(min(max(_wf, 0.0), 0.99) * samples)) if _wf > 0 else 0

    from line_solver.api.io import console as _console
    _console.loop('drawing the sample path: %d samples requested', samples)
    _console_every = max(1, samples // 20)
    for _it in range(samples):
        if (_it + 1) % _console_every == 0:
            _console.iter_line((_it + 1) // _console_every,
                               'simulated %d of %d samples (%.0f%%)',
                               _it + 1, samples, 100.0 * (_it + 1) / samples)
        if math.isfinite(_tmo) and _tmo > 0 and (time.time() - _tmo_start) > _tmo:
            break
        if math.isfinite(_max_pop) and (_it % _POP_CHECK_EVERY) == 0:
            _popvec = _state_nir_vector(sn, cur_state, M, R)
            _perst = np.array([np.sum(_popvec[i * R:(i + 1) * R]) for i in range(M)])
            _perst[~_pop_st_mask] = 0.0   # exclude infinite-population sources
            _tot = float(np.sum(_perst))
            if _tot > _max_pop:
                _ist_max = int(np.argmax(_perst))
                _nd = int(sn.stationToNode[_ist_max])
                _name = sn.nodenames[_nd] if getattr(sn, 'nodenames', None) is not None \
                    and _nd < len(sn.nodenames) else 'station %d' % _ist_max
                raise RuntimeError(
                    "SSA aborted: in-system population reached %s (cap %d) after "
                    "%d events, dominated by '%s'. The model is unstable or a "
                    "server is stalled (check rho<1 and that every service/firing "
                    "distribution is valid). Raise options.ssa_max_pop to override."
                    % (('%.3g' % _tot), int(_max_pop), _it, _name))
        cand_rates = []
        cand_states = []
        cand_meta = [] if record_events else None
        # derived tags of each candidate transition, as [(statefulIndex, class, kind)]
        # rows with kind 0 = START and 1 = PREEMPT; sn.sync carries none, so the
        # trace can only report them from here
        cand_tags = [] if record_events else None
        dep_acc = np.zeros((nstateful, R))
        arv_acc = np.zeros((nstateful, R))
        dly_acc = np.zeros((nstateful, R))
        # Derived START/PREEMPT rates, sampled exactly like the three above: the
        # rate at which the transitions enabled in the current state start a
        # class-r service, or push a class-r job in service back into the buffer.
        start_acc = np.zeros((nstateful, R))
        preempt_acc = np.zeros((nstateful, R))

        # FCR: current aggregate per-class population of each region.
        _xcur = None
        cand_fcr = {} if _fcr_on else None  # cand index -> (region, class, dest, isSwitch)
        if _fcr_on:
            _xcur = [_fcr_regionpop(fi, cur_state) for fi in range(len(_fcr))]

        # Global (Whittle) rate scaling: constant within a state, so one
        # evaluation serves every transition out of it.
        _gd_now = _gd_factor_now(sn, cur_state) if _has_gd else None

        # ---- regular sync actions (mirror solver_ssa_findenabled.m) ----
        for (na, ea, ca, npn, ep, cp, pprob) in acts:
            if not sn.isstateful[na]:
                continue
            isf_a = int(sn.nodeToStateful[na])
            # see _kb/06-solver-catalog.md (SSA: "Serial engine: finite
            # capacity regions (FCR), fork-join, immediate feedback")
            no_promote = False
            if (ea == EventType.DEP and npn == na and sn.isstation[na]
                    and getattr(sn, 'immfeed', None) is not None):
                _ist_if = int(sn.nodeToStation[na])
                if (0 <= _ist_if < sn.immfeed.shape[0]
                        and 0 <= cp < sn.immfeed.shape[1]
                        and sn.immfeed[_ist_if, cp]):
                    no_promote = True
            oa, ra, _pa, sa_tag, pa_tag = _after_event_cached(na, cur_state[isf_a], ea, ca,
                                                              no_promote=no_promote)
            oa = np.atleast_2d(oa)
            ra = np.ravel(ra)
            if oa.size == 0 or ra.size == 0:
                continue
            # PHASE matters as much as DEP: refresh_sync emits phase moves as
            # active station events, so phase-type service would otherwise
            # advance unscaled.
            if _has_gd and sn.isstation[na] and ea in (EventType.DEP, EventType.PHASE):
                ra = ra * _gd_now[int(sn.nodeToStation[na]), ca]
            is_local = (npn >= nnodes)
            isf_p = int(sn.nodeToStateful[npn]) if (not is_local and sn.isstateful[npn]) else -1
            for ia in range(oa.shape[0]):
                rate_ia = ra[ia] if ia < ra.size else 0.0
                if not np.isfinite(rate_ia) or rate_ia == 0:
                    continue
                row_a = oa[ia].astype(float)
                # A delayed hit is the ONLY cache transition that empties the
                # node: the request merges onto the in-flight fetch and is held
                # in block B, so it departs later in the hit class and is
                # otherwise indistinguishable there from a true hit.
                is_merge_a = False
                if isf_a in cache_srv_slice and ea == EventType.READ:
                    _cv = cache_srv_slice[isf_a]
                    _pre = np.ravel(cur_state[isf_a])
                    _e0, _e1 = row_a.size - _cv, _pre.size - _cv
                    is_merge_a = (float(np.sum(row_a[_e0 - R:_e0]))
                                  - float(np.sum(_pre[_e1 - R:_e1]))) == -1.0
                if is_local:
                    # Passive LOCAL probability is the routing mass leaving the
                    # system; without it, feedback networks drain too fast.
                    if callable(pprob):
                        ps_loc = 1.0
                    else:
                        ps_loc = float(pprob) if pprob is not None else 1.0
                    eff_loc = rate_ia * ps_loc
                    if eff_loc == 0:
                        continue
                    nstate = [s.copy() for s in cur_state]
                    nstate[isf_a] = row_a.copy()
                    cand_rates.append(eff_loc)
                    cand_states.append(nstate)
                    if cand_meta is not None:
                        cand_meta.append((isf_a, ca, -1, -1) if ea == EventType.DEP else None)
                        cand_tags.append(_tag_rows(isf_a, sa_tag, pa_tag, ia, R))
                    if ea == EventType.DEP:
                        dep_acc[isf_a, ca] += eff_loc
                    if is_merge_a:
                        dly_acc[isf_a, ca] += eff_loc
                    continue
                if isf_p < 0:
                    continue
                jp = int(sn.nodeToStation[npn]) if npn < nnodes else -1
                # see _kb/06-solver-catalog.md (SSA: "Serial engine: finite
                # capacity regions (FCR)")
                if _fcr_on:
                    if jp >= 0:
                        ja = int(sn.nodeToStation[na]) if na < nnodes else -1
                        _blocked = False
                        _mark = None
                        for _fi in range(len(_fcr)):
                            member_mask = _fcr[_fi][1]
                            iswaitq = _fcr[_fi][8]
                            if (jp < member_mask.size and member_mask[jp]
                                    and (ja < 0 or ja >= member_mask.size or not member_mask[ja])):
                                xn = _xcur[_fi].copy()
                                xn[cp] += 1
                                if _fcr_violates(_fi, xn):
                                    if iswaitq[cp]:
                                        _mark = (_fi, cp, npn, 0)  # park in FIFO
                                    else:
                                        _mark = (_fi, cp, npn, 2)  # DROP: destroyed
                                    break
                            elif (jp < member_mask.size and member_mask[jp]
                                  and 0 <= ja < member_mask.size and member_mask[ja]
                                  and cp != ca):
                                # exit + gated re-entry; DROP destroys on refusal
                                _mark = (_fi, cp, npn, 1 if iswaitq[cp] else 3)
                                break
                        if _blocked:
                            continue
                        if _mark is not None:
                            if callable(pprob):
                                raise RuntimeError('WAITQ finite capacity regions are not supported together with state-dependent routing in SolverSSA.')
                            ps_b = float(pprob) if pprob is not None else 1.0
                            eff_b = rate_ia * ps_b
                            if eff_b <= 0:
                                continue
                            nstate = [s.copy() for s in cur_state]
                            nstate[isf_a] = row_a.copy()
                            cand_rates.append(eff_b)
                            cand_states.append(nstate)
                            cand_fcr[len(cand_states) - 1] = _mark
                            if cand_meta is not None:
                                cand_meta.append((isf_a, ca, isf_p, cp) if ea == EventType.DEP else None)
                                cand_tags.append(_tag_rows(isf_a, sa_tag, pa_tag, ia, R))
                            if ea == EventType.DEP:
                                dep_acc[isf_a, ca] += eff_b
                                arv_acc[isf_p, cp] += eff_b
                            continue
                pin = row_a if npn == na else cur_state[isf_p]
                op, _rp, pp, sp_tag, pp_tag = _after_event_cached(npn, pin, ep, cp)
                op = np.atleast_2d(op)
                pp = np.ravel(pp)
                if op.size == 0:
                    # see _kb/06-solver-catalog.md (SSA: "True BAS in the
                    # serial engine mirrors the CTMC generator arc")
                    if ea == EventType.DEP and _is_bas_station(sn, int(sn.nodeToStation[na])):
                        cv = np.ravel(cur_state[isf_a]).astype(float)
                        if cv.size > 0 and cv[-1] == 0:
                            bv = cv.copy(); bv[-1] = 1
                            nstate = [s.copy() for s in cur_state]
                            nstate[isf_a] = bv
                            cand_rates.append(rate_ia)
                            cand_states.append(nstate)
                            if cand_meta is not None:
                                cand_meta.append(None)
                                cand_tags.append([])
                    continue
                for ip in range(op.shape[0]):
                    row_p = op[ip].astype(float)
                    if callable(pprob):
                        ps = _eval_state_dep_prob(pprob, cur_state, isf_a, row_a, isf_p, row_p, nstateful)
                    else:
                        ps = float(pprob) if pprob is not None else 1.0
                    ps *= float(pp[ip]) if ip < pp.size else 1.0
                    if ps == 0:
                        continue
                    eff = rate_ia * ps
                    nstate = [s.copy() for s in cur_state]
                    nstate[isf_a] = row_a.copy()
                    nstate[isf_p] = row_p.copy()
                    cand_rates.append(eff)
                    cand_states.append(nstate)
                    # START/PREEMPT tags of this arc, weighted like the rate it
                    # carries: they annotate the transition itself. Written for
                    # EVERY action, not only for departures -- a retrial or a
                    # polling switchover starts service without being a DEP, and
                    # the arrival half of a departure is where most starts happen.
                    _accumulate_tags(start_acc, preempt_acc, isf_a, sa_tag, pa_tag, ia, eff, R)
                    _accumulate_tags(start_acc, preempt_acc, isf_p, sp_tag, pp_tag, ip, eff, R)
                    if cand_meta is not None:
                        cand_meta.append((isf_a, ca, isf_p, cp) if ea == EventType.DEP else None)
                        cand_tags.append(_tag_rows(isf_a, sa_tag, pa_tag, ia, R)
                                         + _tag_rows(isf_p, sp_tag, pp_tag, ip, R))
                    if ea == EventType.DEP:
                        # see _kb/06-solver-catalog.md (SSA: "Python serial
                        # engine: state-construction traps") -- ARV convention
                        refused = (jp >= 0 and _sched_name(sn, jp) != 'EXT'
                                   and np.array_equal(row_p, np.ravel(pin).astype(float)))
                        # see _kb/06-solver-catalog.md (SSA: "Python serial
                        # engine: state-construction traps") -- offered-load convention
                        if not refused or is_physical_cap[jp, cp]:
                            dep_acc[isf_a, ca] += eff
                            arv_acc[isf_p, cp] += eff
                    if is_merge_a:
                        dly_acc[isf_a, ca] += eff

        # ---- SPN global sync events (mirror solver_ssa.m:244-282) ----
        for glevent in gsync_events:
            if not getattr(glevent, 'active', None):
                continue
            gind = int(glevent.active[0].node)
            glspace = []
            for isf in range(nstateful):
                glspace.append(np.ravel(cur_state[isf]).astype(float) if cur_state[isf] is not None
                               else np.zeros(0))
            res = after_global_event(sn, gind, glspace, glevent, False)
            outrate = np.ravel(res.outrate)
            outprob = np.ravel(res.outprob)
            for io in range(outrate.size):
                eff = float(outrate[io]) * (float(outprob[io]) if io < outprob.size else 1.0)
                if eff <= 0:
                    continue
                gl_io = res.outglspace[io]
                nstate = [s.copy() for s in cur_state]
                for isf in range(nstateful):
                    new_row = np.asarray(gl_io[isf], dtype=float).ravel()
                    if new_row.size and not np.array_equal(new_row, np.ravel(cur_state[isf])):
                        nstate[isf] = new_row
                cand_rates.append(eff)
                cand_states.append(nstate)
                if cand_meta is not None:
                    # SPN token moves have no single (src,dst) station pair; the
                    # sample-path trace covers queueing transitions only.
                    cand_meta.append(None)
                    cand_tags.append([])
                # ModeEvent carries the place node index and class directly;
                # use them rather than decoding a linear index (old lin%nnodes decode
                # only coincided for R=1).
                for pev in getattr(glevent, 'passive', []) or []:
                    pev_node = int(pev.node)
                    pev_class = int(getattr(pev, 'job_class', 0))
                    if pev_node >= len(sn.nodeToStateful):
                        continue
                    p_isf = int(sn.nodeToStateful[pev_node])
                    if p_isf < 0:
                        continue
                    if pev.event == EventType.PRE:
                        dep_acc[p_isf, pev_class] += eff
                    elif pev.event == EventType.POST:
                        arv_acc[p_isf, pev_class] += eff

        # Fork firing (native fork-join, solver_ssa.m:311-343): consumes the
        # parent held at the Fork, emits one sibling per branch atomically.
        for k in range(FJ):
            glspace = [np.ravel(cur_state[isf]).astype(float) if cur_state[isf] is not None
                       else np.zeros(0) for isf in range(nstateful)]
            fj_states, fj_rate, fj_prob = after_fj_event(sn, fjsync[k], glspace, True)
            if not fj_states:
                continue
            eff = float(fj_rate[0]) * (float(fj_prob[0]) if len(fj_prob) else 1.0)
            if eff <= 0:
                continue
            gl_io = fj_states[0]
            nstate = [s.copy() for s in cur_state]
            for isf in range(nstateful):
                new_row = np.asarray(gl_io[isf], dtype=float).ravel()
                if not np.array_equal(new_row, np.ravel(cur_state[isf])):
                    nstate[isf] = new_row
            cand_rates.append(eff)
            cand_states.append(nstate)
            if cand_meta is not None:
                cand_meta.append(None)
            fjentry = fjsync[k]
            isf_fork = int(sn.nodeToStateful[int(fjentry['fork'])])
            dep_acc[isf_fork, int(fjentry['class'])] += eff
            branchheads = np.asarray(fjentry['branchheads']).ravel()
            auxclasses = np.asarray(fjentry['auxclasses']).ravel()
            for b in range(len(branchheads)):
                isf_bh = int(sn.nodeToStateful[int(branchheads[b])])
                arv_acc[isf_bh, int(auxclasses[b])] += eff

        if not cand_rates:
            break   # deadlock / absorbing
        rates = np.asarray(cand_rates, dtype=float)
        tot = float(rates.sum())
        if tot <= 0:
            break

        # GILLESPIE DIRECT, IN THE REFERENCE'S DRAW ORDER: the ACTION first, the
        # holding time second. `solver_ssa.m:578` selects the transition against
        # `cumsum(enabled_rates)/tot_rate` and only then, at line 615, draws
        # `dt = -log(rand)/tot_rate`. Drawing them the other way round consumes
        # the same MT19937 stream in the opposite order, and from the second step
        # onward the two engines are on DIFFERENT sample paths -- distributionally
        # identical, so nothing looks wrong until a seeded golden is compared. The
        # streams themselves already agree: `rng(s,'twister')` and
        # `np.random.seed(s)` both init_genrand(s) and both take 53 bits per
        # double, so seed 1 gives 0.417022004702574 in either (s = 0 is MATLAB's
        # one special case, remapped to 5489).
        #
        # The comparison mirrors the reference exactly, normalized cumulative
        # against a raw uniform, rather than `rand*tot` against the unnormalized
        # sum: the two differ in the last bits and that is enough to pick a
        # different transition at a boundary.
        cum = np.cumsum(rates) / tot
        sel = int(np.searchsorted(cum, np.random.random()))
        if sel >= len(cand_states):
            sel = len(cand_states) - 1

        dt = -np.log(np.random.random()) / tot

        if n_drop and _it == n_drop:
            # see _kb/06-solver-catalog.md (SSA: "Warmup discard, and cache
            # hit/miss accounting")
            tally.clear()
            total_time = 0.0
            if record_events:
                event_log.clear()
                state_log.clear()

        key = _state_key(cur_state)
        ent = tally.get(key)
        if ent is None:
            tally[key] = [dt, _state_nir_vector(sn, cur_state, M, R),
                          dep_acc.ravel().copy(), arv_acc.ravel().copy(),
                          dly_acc.ravel().copy(), 1,
                          start_acc.ravel().copy(), preempt_acc.ravel().copy()]
        else:
            ent[0] += dt
            if cache_srv_slice:
                # see _kb/09-ldes-and-cache.md (SSA: the merge rate needs EVERY
                # visit). Unlike a DEP rate, the merge rate is random given the
                # state, so one visit is a single Bernoulli draw whose variance
                # does not shrink with the sample count.
                ent[4] += dly_acc.ravel()
                ent[5] += 1
        total_time += dt
        cur_time += dt

        if record_events:
            meta = cand_meta[sel]
            if meta is not None:
                src_isf, src_k, dst_isf, dst_k = meta
                event_log.append((total_time,
                                  _sf_to_station(src_isf), int(src_k),
                                  _sf_to_station(dst_isf), int(dst_k)))
                state_log.append(_state_nir_vector(sn, cur_state, M, R).reshape(M, R))
            # PREEMPT before START at the same instant: the victim leaves the
            # server before the job that displaced it takes it. Emitted even when
            # meta is None, since a START can ride on a non-DEP action.
            if sel < len(cand_tags):
                for (isf_t, cls_t, kind_t) in sorted(cand_tags[sel], key=lambda x: -x[2]):
                    tag_log.append((total_time, _sf_to_station(isf_t), int(cls_t), int(kind_t)))

        cur_state = cand_states[sel]

        # FCR WAITQ bookkeeping: park blocked entries, release FIFO heads,
        # resolve pending class-switch re-entries
        if _fcr_on:
            _mark = cand_fcr.get(sel)
            if _mark is not None and _mark[3] == 0:
                _fi, _cls, _dest, _ = _mark
                _fcr_buf[_fi].append(_dest * R + _cls)
            # _mark[3] == 2: DROP, the refused job was destroyed (active only)
            cur_state = _fcr_release(cur_state)
            if _mark is not None and _mark[3] in (1, 3):
                _fi, _cls, _dest, _ = _mark
                x = _fcr_regionpop(_fi, cur_state)
                xn = x.copy()
                xn[_cls] += 1
                _admitted = False
                if not _fcr_violates(_fi, xn):
                    isf_d = int(sn.nodeToStateful[_dest])
                    od, _rd, _pd, _sd2, _pd2 = _after_event_cached(_dest, cur_state[isf_d], EventType.ARV, _cls)
                    od = np.atleast_2d(od)
                    if od.size > 0:
                        cur_state[isf_d] = od[0].astype(float)
                        _admitted = True
                if not _admitted and _mark[3] == 1:
                    _fcr_buf[_fi].append(_dest * R + _cls)
                # _mark[3] == 3 refused: DROP, the switching job is destroyed

        if cur_time > t_end:
            break

    # ---- assemble per-unique-state arrays ----
    nuniq = len(tally)
    pi = np.zeros(nuniq)
    SSq = np.zeros((nuniq, M * R))
    depRates = np.zeros((nuniq, nstateful, R))
    arvRates = np.zeros((nuniq, nstateful, R))
    dlyRates = np.zeros((nuniq, nstateful, R))
    startRates = np.zeros((nuniq, nstateful, R))
    preemptRates = np.zeros((nuniq, nstateful, R))
    for s, (_key, ent) in enumerate(tally.items()):
        pi[s] = ent[0]
        SSq[s, :] = ent[1]
        depRates[s] = ent[2].reshape(nstateful, R)
        arvRates[s] = ent[3].reshape(nstateful, R)
        dlyRates[s] = ent[4].reshape(nstateful, R) / max(ent[5], 1)
        # the tag rates are a deterministic function of the state too, so one
        # visit gives them exactly, as for the departure and arrival rates
        if len(ent) > 7:
            startRates[s] = ent[6].reshape(nstateful, R)
            preemptRates[s] = ent[7].reshape(nstateful, R)
    if pi.sum() > 0:
        pi = pi / pi.sum()
    return (pi, SSq, arvRates, depRates, dlyRates, total_time, event_log, state_log,
            startRates, preemptRates, tag_log)


def _tag_rows(isf, start_tag, preempt_tag, row, R):
    """One (statefulIndex, class, kind) row per tagged job on this successor,
    kind 0 = START and 1 = PREEMPT."""
    out = []
    if isf is None or isf < 0:
        return out
    for kind, tag in ((0, start_tag), (1, preempt_tag)):
        if tag is None:
            continue
        arr = np.atleast_2d(np.asarray(tag, dtype=float))
        if row >= arr.shape[0]:
            continue
        for r in range(min(R, arr.shape[1])):
            for _ in range(int(arr[row, r])):
                out.append((isf, r, kind))
    return out


def _accumulate_tags(start_acc, preempt_acc, isf, start_tag, preempt_tag, row, w, R):
    """Add the START/PREEMPT counts of one successor row to the per-state rate
    accumulators, weighted by the rate of the arc that carries them."""
    if isf is None or isf < 0 or w == 0:
        return
    for tag, acc in ((start_tag, start_acc), (preempt_tag, preempt_acc)):
        if tag is None:
            continue
        arr = np.atleast_2d(np.asarray(tag, dtype=float))
        if row >= arr.shape[0]:
            continue
        for r in range(min(R, arr.shape[1])):
            if arr[row, r] != 0:
                acc[isf, r] += w * arr[row, r]


def _analyze(sn, pi, SSq, arvRates, depRates, user_cap=None, user_classcap=None):
    """Compute XN/UN/QN/RN/TN/CN (port of solver_ssa_analyzer_serial.m:42-112).

    ``user_cap``/``user_classcap`` are the USER capacities snapshotted before the
    solver_ssa preamble folds the state-space cutoff into sn.cap/classcap (which
    makes every open class look capacity-constrained). They drive the BUG-12
    canDropClass guard so utilization uses carried (not offered) load at
    finite-capacity stations. Fall back to the (folded) sn values if not given.
    """
    M, R = sn.nstations, sn.nclasses
    S = np.asarray(sn.nservers, dtype=float).reshape(-1)
    rates = np.asarray(sn.rates, dtype=float).reshape(M, R)
    NK = np.asarray(sn.njobs, dtype=float).reshape(-1)
    refstat = np.asarray(sn.refstat, dtype=int).reshape(-1)
    st2sf = np.asarray(sn.stationToStateful, dtype=int).reshape(-1)
    ucap = np.asarray(user_cap if user_cap is not None else sn.cap,
                      dtype=float).reshape(-1)
    uclasscap = np.asarray(user_classcap if user_classcap is not None else sn.classcap,
                           dtype=float).reshape(M, R)

    XN = np.zeros(R)
    QN = np.zeros((M, R))
    UN = np.zeros((M, R))
    TN = np.zeros((M, R))
    RN = np.zeros((M, R))
    CN = np.zeros(R)

    for r in range(R):
        ref_sf = int(st2sf[int(refstat[r])])
        XN[r] = float(pi @ depRates[:, ref_sf, r])

    for ist in range(M):
        isf = int(st2sf[ist])
        sname = _sched_name(sn, ist)
        for r in range(R):
            TN[ist, r] = float(pi @ depRates[:, isf, r])
            QN[ist, r] = float(pi @ SSq[:, ist * R + r])
        if sname == 'INF':
            UN[ist, :] = QN[ist, :]
        else:
            # see _kb/06-solver-catalog.md (SSA: "Utilization: pre-preamble
            # capacity snapshot, carried vs offered load")
            eff_c = S[ist]
            if getattr(sn, 'lldscaling', None) is not None:
                _lld = np.asarray(sn.lldscaling)
                if ist < _lld.shape[0]:
                    _peak = float(np.max(_lld[ist]))
                    if _peak > eff_c:
                        eff_c = _peak
            # Class-dependent stations divide by the declared per-class peak
            # rate (sn.cdscalingpeak), so Util = T*S/peak (matching CTMC/NC).
            _cd = getattr(sn, 'cdscaling', None)
            _is_cd = _cd is not None and (
                (_cd.get(ist) if isinstance(_cd, dict)
                 else (_cd[ist] if ist < len(_cd) else None)) is not None)
            # Joint-dependent stations use the same Util=T*S/peak convention; the
            # effective peak is the product of the declared cd and jd peaks.
            _jd = getattr(sn, 'jdscaling', None)
            _is_jd = _jd is not None and (
                (_jd.get(ist) if isinstance(_jd, dict)
                 else (_jd[ist] if ist < len(_jd) else None)) is not None)
            # A global (Whittle) dependence rescales the service rate the same
            # way, so the peak it declares normalizes the utilization too.
            _is_gd = getattr(sn, 'gdscaling', None) is not None
            for r in range(R):
                mu = rates[ist, r]
                if np.isfinite(mu) and mu > 0:
                    if _is_cd or _is_jd or _is_gd:
                        cdiv = 1.0
                        if _is_cd:
                            cdiv *= sn.cdscalingpeak[ist, r]
                        if _is_jd:
                            cdiv *= sn.jdscalingpeak[ist, r]
                        if _is_gd:
                            cdiv *= sn.gdscalingpeak[ist, r]
                    else:
                        cdiv = eff_c
                    if not (cdiv > 0):
                        UN[ist, r] = 0.0
                        continue
                    # see _kb/06-solver-catalog.md (SSA: "Utilization:
                    # pre-preamble capacity snapshot, carried vs offered load")
                    can_drop = (not np.isfinite(NK[r])) and (
                        np.isfinite(ucap[ist])
                        or (r < uclasscap.shape[1] and np.isfinite(uclasscap[ist, r])))
                    if can_drop:
                        UN[ist, r] = TN[ist, r] / mu / cdiv
                    else:
                        arv = float(pi @ arvRates[:, isf, r])
                        UN[ist, r] = arv / mu / cdiv

    for r in range(R):
        for ist in range(M):
            RN[ist, r] = QN[ist, r] / TN[ist, r] if TN[ist, r] > 0 else 0.0
        CN[r] = NK[r] / XN[r] if XN[r] > 0 else 0.0

    for arr in (QN, UN, TN, RN):
        arr[~np.isfinite(arr)] = 0.0
    XN[~np.isfinite(XN)] = 0.0
    CN[~np.isfinite(CN)] = 0.0
    return XN, UN, QN, RN, TN, CN


def _compute_cache_hitprob(sn, pi, depRates, dlyRates=None):
    """Populate sn.nodeparam[cache].actualhitprob/actualmissprob/
    actualdelayedhitprob from the simulated cache departure rates.

    Mirrors the CTMC analyzer (_compute_cache_hit_miss_probs): a Cache is not a
    station, so its hit/miss split is read off the stateful-indexed departure
    rates, actualhitprob = pi @ depRates[:, cache_sf, hit] / (hit + miss). The
    SSA cache READ already emits the correct per-state hit/miss rates, but
    without this the downstream node-throughput expansion
    (sn_get_node_tput_from_tput) has no actualhitprob to use and falls back to
    the static routing matrix, yielding an even 0.5/0.5 hit/miss split.

    With a retrieval system the hit-class departure rate is (true hits +
    delayed hits), because a request merged onto an in-flight fetch is released
    in the hit class when that fetch completes. `dlyRates` carries the rate of
    the merge transitions, which is what separates the two.
    """
    if sn.nodeparam is None or sn.nodetype is None or depRates is None:
        return
    K = sn.nclasses
    n2sf = sn.nodeToStateful
    for ind in range(sn.nnodes):
        if ind >= len(sn.nodetype) or sn.nodetype[ind] != NodeType.CACHE:
            continue
        if ind not in sn.nodeparam:
            continue
        npar = sn.nodeparam[ind]
        hitclass = getattr(npar, 'hitclass', None)
        missclass = getattr(npar, 'missclass', None)
        if hitclass is None or missclass is None:
            continue
        hitclass = np.atleast_1d(hitclass).flatten()
        missclass = np.atleast_1d(missclass).flatten()
        cache_sf = int(n2sf[ind]) if (n2sf is not None and ind < len(n2sf)) else -1
        if not (0 <= cache_sf < depRates.shape[1]):
            continue
        ahp = np.zeros(K)
        amp = np.zeros(K)
        adhp = np.zeros(K)
        has_delayed = False
        for oc in range(len(hitclass)):
            h = int(hitclass[oc])
            m = int(missclass[oc])
            if h < 0 or m < 0 or h >= K or m >= K:
                continue
            t_hit = float(pi @ depRates[:, cache_sf, h])
            t_miss = float(pi @ depRates[:, cache_sf, m])
            t_dly = 0.0
            if dlyRates is not None and oc < dlyRates.shape[2]:
                t_dly = float(pi @ dlyRates[:, cache_sf, oc])
            if t_hit + t_miss > 0:
                # t_hit already contains the released delayed hits, so carve
                # them out rather than adding a fourth share.
                ahp[oc] = max(t_hit - t_dly, 0.0) / (t_hit + t_miss)
                amp[oc] = t_miss / (t_hit + t_miss)
                adhp[oc] = t_dly / (t_hit + t_miss)
                if t_dly > 0:
                    has_delayed = True
        npar.actualhitprob = ahp
        npar.actualmissprob = amp
        if has_delayed:
            npar.actualdelayedhitprob = adhp


def solver_ssa_run(sn, options, method='serial', seed=None):
    """Run the serial engine and return a populated SolverSSAReturn."""
    from .handler import SolverSSAReturn
    import copy as _copy
    t0 = time.time()
    # see _kb/06-solver-catalog.md (SSA: "Python serial engine: state-
    # construction traps and cross-solver sn mutation")
    _orig = {f: _copy.deepcopy(getattr(sn, f))
             for f in ('nservers', 'classcap', 'cap', 'mu', 'phi', 'pie', 'proc')}
    try:
        (pi, SSq, arvRates, depRates, dlyRates, total_time, event_log, state_log,
         startRates, preemptRates, tag_log) = solver_ssa_serial(sn, options, seed=seed)
        _compute_cache_hitprob(sn, pi, depRates, dlyRates)
        # Pre-preamble USER capacities so the BUG-12 canDropClass guard sees
        # real finite-capacity stations, not the cutoff-derived truncation.
        XN, UN, QN, RN, TN, CN = _analyze(sn, pi, SSq, arvRates, depRates,
                                          user_cap=_orig['cap'],
                                          user_classcap=_orig['classcap'])
        M, R = sn.nstations, sn.nclasses
        AN = np.zeros((M, R))
        st2sf = np.asarray(sn.stationToStateful, dtype=int).reshape(-1)
        StartN = np.zeros((M, R))
        PreemptN = np.zeros((M, R))
        for ist in range(M):
            isf = int(st2sf[ist])
            for r in range(R):
                AN[ist, r] = float(pi @ arvRates[:, isf, r])
                # same time average as TN, over the derived tag rates
                StartN[ist, r] = float(pi @ startRates[:, isf, r])
                PreemptN[ist, r] = float(pi @ preemptRates[:, isf, r])
    finally:
        for f, v in _orig.items():
            setattr(sn, f, v)
    res = SolverSSAReturn()
    res.Q, res.U, res.R, res.T, res.C, res.X, res.A = QN, UN, RN, TN, CN, XN, AN
    res.total_time = total_time
    res.runtime = time.time() - t0
    res.method = method
    res.samples = int(getattr(options, 'samples', 0))
    res.event_log = event_log
    res.state_log = state_log
    # Derived START/PREEMPT rates: annotations on the transitions the engine
    # already fires, so they add no getAvgTable column.
    res.startRate = StartN
    res.preemptRate = PreemptN
    res.tag_log = tag_log
    _tmo = float(getattr(options, 'timeout', float('inf')))
    res.timedOut = math.isfinite(_tmo) and _tmo > 0 and res.runtime > _tmo
    if res.timedOut:
        import warnings as _warnings
        _warnings.warn("SolverSSA stopped after the wall-clock time budget "
                       "(timeout=%gs) was exceeded; returning the interim solution." % _tmo)
    return res


def solver_ssa_parallel(sn, options):
    """Worker-count-invariant parallel SSA (mirrors MATLAB
    solver_ssa_analyzer_parallel). The simulation budget is split into a FIXED
    number of independent replications R (default 8, MATLAB's config.nreplicas
    default): replication r simulates ceil(samples/R) events, seeded
    deterministically with base_seed+r, and the per-replication point estimates
    are averaged. Because R and the per-replication seed/budget are fixed, the
    returned averages depend only on (seed, samples, R) and are independent of
    any worker/thread count, matching MATLAB."""
    import copy as _copy
    nrep = int(getattr(options, 'nreplicas', 0)
               or getattr(options, 'numthreads', 0) or 8)
    nrep = max(1, nrep)
    base_seed = int(getattr(options, 'seed', 0) or 0)
    total_samples = int(getattr(options, 'samples', 0) or 0)
    per_rep_samples = (int(np.ceil(total_samples / nrep))
                       if total_samples > 0 else total_samples)
    accums = None
    t0 = time.time()
    for rep in range(nrep):
        rep_options = _copy.copy(options)
        rep_options.samples = per_rep_samples
        r = solver_ssa_run(sn, rep_options, method='parallel', seed=base_seed + rep)
        fields = (r.Q, r.U, r.R, r.T, r.C, r.X, r.A)
        if accums is None:
            accums = [np.asarray(f, dtype=float).copy() for f in fields]
        else:
            for i, f in enumerate(fields):
                accums[i] += np.asarray(f, dtype=float)
    from .handler import SolverSSAReturn
    res = SolverSSAReturn()
    res.Q, res.U, res.R, res.T, res.C, res.X, res.A = [a / nrep for a in accums]
    res.runtime = time.time() - t0
    res.method = 'parallel'
    res.samples = per_rep_samples * nrep
    return res
