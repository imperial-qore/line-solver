"""
Main afterEvent dispatch for CTMC state space generation.

This module provides the core event dispatch mechanism used by the
sync-action-based CTMC solver. Given a stateful node, an event type,
and a current state, it computes all possible successor states with
their rates and probabilities.

Port from JAR State.java:afterEvent/afterEventHashed.
"""

import numpy as np
from typing import Dict, Tuple, Optional
from ...constants import EventType, ProcessType
from ...lang.base import NodeType
from ...lang.base import SchedStrategy


def build_space_hash(sn) -> Dict[int, Dict[tuple, int]]:
    """
    Build hash maps for O(1) state lookup for each stateful node.

    Args:
        sn: NetworkStruct with sn.space populated

    Returns:
        Dict mapping stateful node index -> {state_tuple: row_index}
    """
    hash_maps = {}
    if not hasattr(sn, 'space') or sn.space is None:
        return hash_maps

    for isf in range(sn.nstateful):
        if isf in sn.space and sn.space[isf] is not None:
            space = np.atleast_2d(sn.space[isf])
            row_map = {}
            for row_idx in range(space.shape[0]):
                key = tuple(space[row_idx].astype(float))
                row_map[key] = row_idx
            hash_maps[isf] = row_map
    return hash_maps


def get_hash(sn, ind, outspace, hash_maps=None):
    """
    Find row indices of output states in the precomputed state space.

    Args:
        sn: NetworkStruct
        ind: Node index (0-based)
        outspace: Output state matrix, shape (n_rows, n_cols)
        hash_maps: Optional precomputed hash maps from build_space_hash()

    Returns:
        np.ndarray of shape (n_rows,) with hash indices, or -1 for not found
    """
    isf = int(sn.nodeToStateful[ind])
    outspace = np.atleast_2d(outspace)
    n_rows = outspace.shape[0]
    outhash = -1 * np.ones(n_rows, dtype=int)

    if hash_maps is not None and isf in hash_maps:
        row_map = hash_maps[isf]
        for i in range(n_rows):
            key = tuple(outspace[i].astype(float))
            if key in row_map:
                outhash[i] = row_map[key]
    else:
        # Fallback: linear search
        if isf in sn.space and sn.space[isf] is not None:
            space = np.atleast_2d(sn.space[isf])
            for i in range(n_rows):
                for j in range(space.shape[0]):
                    if np.array_equal(outspace[i], space[j]):
                        outhash[i] = j
                        break

    return outhash


def after_event_init(sn):
    """
    Precompute the loop-invariant setup of after_event.

    Hot callers (the SSA serial Gillespie loop re-evaluates every sync at
    every step) pass the returned ctx to after_event to skip the per-call
    derivation of ismkvmodclass and lldscaling; semantics are identical.

    Build ctx AFTER any caller-side refresh of sn (phase fields, capacity
    recompute) so it reflects the same sn passed to after_event.

    Port of MATLAB State.afterEventInit.
    """
    M = sn.nstations
    R = sn.nclasses
    ctx = {}

    lldscaling = sn.lldscaling if hasattr(sn, 'lldscaling') and sn.lldscaling is not None else None
    if lldscaling is None or (isinstance(lldscaling, np.ndarray) and lldscaling.size == 0):
        lldlimit = max(int(sn.nclosedjobs), 1) if hasattr(sn, 'nclosedjobs') else 1
        lldscaling = np.ones((M, lldlimit))
    else:
        lldlimit = lldscaling.shape[1] if isinstance(lldscaling, np.ndarray) and lldscaling.ndim >= 2 else 1
    ctx['lldscaling'] = lldscaling
    ctx['lldlimit'] = lldlimit

    ismkvmodclass_by_ist = {}
    if hasattr(sn, 'procid') and sn.procid is not None:
        map_val = int(ProcessType.MAP.value) if hasattr(ProcessType.MAP, 'value') else int(ProcessType.MAP)
        mmpp_val = int(ProcessType.MMPP2.value) if hasattr(ProcessType.MMPP2, 'value') else int(ProcessType.MMPP2)
        for ist in range(M):
            ismkvmodclass = np.zeros(R)
            for r in range(R):
                if isinstance(sn.procid, dict):
                    pt = sn.procid.get(ist, {}).get(r, None)
                else:
                    try:
                        pt = sn.procid[ist, r]
                    except (TypeError, KeyError, IndexError):
                        pt = None
                if pt is not None:
                    pt_val = int(pt.value) if hasattr(pt, 'value') else int(pt)
                    if pt_val in (map_val, mmpp_val):
                        ismkvmodclass[r] = 1
            ismkvmodclass_by_ist[ist] = ismkvmodclass
    else:
        for ist in range(M):
            ismkvmodclass_by_ist[ist] = np.zeros(R)
    ctx['ismkvmodclass'] = ismkvmodclass_by_ist
    return ctx


def after_event(sn, ind, inspace, event, job_class, is_simulation=False, ctx=None, no_promote=False):
    """
    Compute successor states after an event at a stateful node.

    Decomposes the state vector, dispatches to the appropriate handler
    based on node type, and returns possible output states with rates
    and probabilities.

    Port from JAR State.java:afterEvent().

    Args:
        sn: NetworkStruct
        ind: Node index (0-based)
        inspace: Input state row(s), shape (n_cols,) or (n_rows, n_cols)
        event: EventType
        job_class: Job class index (0-based)
        ctx: Optional loop-invariant context from after_event_init(sn)

    Returns:
        Tuple of (outspace, outrate, outprob):
        - outspace: np.ndarray of output state rows
        - outrate: np.ndarray of transition rates
        - outprob: np.ndarray of transition probabilities
    """
    inspace = np.atleast_2d(inspace).astype(float)
    M = sn.nstations
    R = sn.nclasses

    # Join synchronization stations on FJ-augmented structs carry a plain
    # per-class count state handled by a dedicated branch (they are stations,
    # but must bypass the buffer/phase slicing and after_event_station below)
    if getattr(sn, 'isfjaugmented', False) and sn.nodetype[ind] == NodeType.JOIN:
        from .after_event_join import after_event_join
        outspace, outrate, outprob = after_event_join(
            sn, ind, inspace, event, job_class, is_simulation)
        if outspace is None:
            return _no_successor(R)
        # a Join holds no server: nothing starts or is preempted there
        n = np.atleast_2d(outspace).shape[0]
        return outspace, outrate, outprob, np.zeros((n, R)), np.zeros((n, R))

    # Determine node type and decompose state
    if sn.isstation[ind]:
        ist = int(sn.nodeToStation[ind])
        # Pass-and-swap / order-independent stations use a dedicated ordered-list
        # representation; handle them before the server/buffer split.
        if sn.sched[ist] == SchedStrategy.PAS:
            Vp = int(np.sum(sn.nvars[ind])) if sn.nvars is not None else 0
            from .after_event_station import after_event_station_pas
            outspace, outrate, outprob = after_event_station_pas(
                sn, ind, ist, inspace, event, job_class, R, Vp)
            # after_event returns FIVE values since the start/preempt tags were
            # added; the PAS handler still returns three, and a PAS station
            # tags nothing (its service is exponential and order-independent,
            # so no job starts or is preempted at a phase boundary). Padded
            # here, exactly as the Join branch above pads its own three.
            n = np.atleast_2d(outspace).shape[0]
            return outspace, outrate, outprob, np.zeros((n, R)), np.zeros((n, R))
        K = np.array(sn.phasessz[ist], dtype=int)
        Ks = np.array(sn.phaseshift[ist], dtype=int)

        if K[job_class] == 0:
            return _no_successor(R)

        V = int(np.sum(sn.nvars[ind])) if sn.nvars is not None else 0
        sumK = int(np.sum(K))
        n_rows = inspace.shape[0]
        n_cols = inspace.shape[1]

        nt = sn.nodetype[ind]
        nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
        place_val = int(NodeType.PLACE.value) if hasattr(NodeType.PLACE, 'value') else int(NodeType.PLACE)

        if nt_val == place_val:
            # Place node: [buffer(R) | server(sumK)]
            space_var = np.zeros((n_rows, 0))
            space_buf = inspace[:, :R]
            space_srv = inspace[:, R:R + sumK] if n_cols >= R + sumK else np.zeros((n_rows, sumK))
        else:
            # Regular station: [buf | srv | var]
            if V > 0:
                space_var = inspace[:, -V:]
            else:
                space_var = np.zeros((n_rows, 0))
            srv_start = n_cols - sumK - V
            srv_end = n_cols - V
            space_srv = inspace[:, srv_start:srv_end]
            buf_end = srv_start
            space_buf = inspace[:, :buf_end] if buf_end > 0 else np.zeros((n_rows, 0))

        hasOnlyExp = int(np.max(K)) == 1

        # Build parameter dicts (precomputed by after_event_init when ctx given)
        if ctx is not None:
            ismkvmodclass = ctx['ismkvmodclass'][ist]
            lldscaling = ctx['lldscaling']
            lldlimit = ctx['lldlimit']
        else:
            ismkvmodclass = np.zeros(R)
            if hasattr(sn, 'procid') and sn.procid is not None:
                map_val = int(ProcessType.MAP.value) if hasattr(ProcessType.MAP, 'value') else int(ProcessType.MAP)
                mmpp_val = int(ProcessType.MMPP2.value) if hasattr(ProcessType.MMPP2, 'value') else int(ProcessType.MMPP2)
                for r in range(R):
                    if isinstance(sn.procid, dict):
                        pt = sn.procid.get(ist, {}).get(r, None)
                    else:
                        try:
                            pt = sn.procid[ist, r]
                        except (TypeError, KeyError, IndexError):
                            pt = None
                    if pt is not None:
                        pt_val = int(pt.value) if hasattr(pt, 'value') else int(pt)
                        if pt_val in (map_val, mmpp_val):
                            ismkvmodclass[r] = 1

            lldscaling = sn.lldscaling if hasattr(sn, 'lldscaling') and sn.lldscaling is not None else None
            if lldscaling is None or (isinstance(lldscaling, np.ndarray) and lldscaling.size == 0):
                lldlimit = max(int(sn.nclosedjobs), 1) if hasattr(sn, 'nclosedjobs') else 1
                lldscaling = np.ones((M, lldlimit))
            else:
                lldlimit = lldscaling.shape[1] if isinstance(lldscaling, np.ndarray) and lldscaling.ndim >= 2 else 1

        cdscaling = sn.cdscaling if hasattr(sn, 'cdscaling') else None

        # Fold joint-dependence handles (sn.jdscaling, non-product-form eta_i)
        # into an effective per-station handle eta_i(n)*beta_i(n), so the state
        # machinery (which reads a single cdscaling) applies both. cd and jd are
        # evaluated identically; when only one is present the product reproduces
        # the single-mechanism case.
        jdscaling = getattr(sn, 'jdscaling', None)
        if jdscaling is not None and len(jdscaling) > 0:
            def _cd_entry(container, i):
                if container is None:
                    return None
                if isinstance(container, dict):
                    return container.get(i, None)
                return container[i] if i < len(container) else None
            eff = list(cdscaling) if (cdscaling is not None and not isinstance(cdscaling, dict)) \
                else ([None] * M)
            if isinstance(cdscaling, dict):
                eff = [cdscaling.get(i, None) for i in range(M)]
            for i in range(M):
                jd_h = _cd_entry(jdscaling, i)
                if jd_h is None:
                    continue
                cd_h = eff[i] if i < len(eff) else None
                if cd_h is None:
                    eff[i] = jd_h
                else:
                    eff[i] = (lambda cf, jf: (lambda n: np.asarray(cf(n), dtype=float) * np.asarray(jf(n), dtype=float)))(cd_h, jd_h)
            cdscaling = eff

        from .after_event_station import after_event_station
        return after_event_station(
            sn, ind, inspace, event, job_class,
            M, R, ist, K, Ks, hasOnlyExp,
            sn.mu, sn.phi, sn.pie, sn.proc, ismkvmodclass,
            lldscaling, lldlimit, cdscaling,
            sn.cap if hasattr(sn, 'cap') else None,
            sn.classcap if hasattr(sn, 'classcap') else None,
            V, space_buf, space_srv, space_var, is_simulation, no_promote)

    elif sn.isstateful[ind]:
        V = int(np.sum(sn.nvars[ind])) if sn.nvars is not None else 0
        n_rows = inspace.shape[0]
        n_cols = inspace.shape[1]

        nt = sn.nodetype[ind]
        nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
        router_val = int(NodeType.ROUTER.value) if hasattr(NodeType.ROUTER, 'value') else int(NodeType.ROUTER)
        cache_val = int(NodeType.CACHE.value) if hasattr(NodeType.CACHE, 'value') else int(NodeType.CACHE)
        trans_val = int(NodeType.TRANSITION.value) if hasattr(NodeType.TRANSITION, 'value') else int(NodeType.TRANSITION)
        fork_val = int(NodeType.FORK.value) if hasattr(NodeType.FORK, 'value') else int(NodeType.FORK)

        # State variable extraction (always at the end)
        if V > 0:
            space_var = inspace[:, -V:]
        else:
            space_var = np.zeros((n_rows, 0))

        if nt_val == fork_val:
            # Stateful Fork (FJ tag-augmented copies only): per-class count of
            # parent jobs held before the fork firing (sn.fjsync)
            srv_start = max(0, n_cols - R - V)
            srv_end = n_cols - V
            space_srv = inspace[:, srv_start:srv_end] if srv_end > srv_start else np.zeros((n_rows, R))
            space_buf = np.zeros((n_rows, 0))
            from .after_event_fork import after_event_fork
            return _untagged(after_event_fork(
                sn, ind, event, job_class, space_buf, space_srv, space_var, is_simulation), R)

        if nt_val == trans_val:
            # Transition: [idle(nmodes) | phases(sumK) | fired(nmodes) | var(V)]
            nparam = sn.nodeparam[ind] if sn.nodeparam is not None and ind in sn.nodeparam else None
            nmodes = nparam.get('nmodes', 0) if isinstance(nparam, dict) else getattr(nparam, 'nmodes', 0)

            # Get K (firing phases per mode) from nodeparam. firingphases[m] is
            # NaN for a non-Markovian firing distribution (e.g. Pareto, which has
            # no finite phase-type representation). Resolve NaN entries in the
            # FLOAT domain from firingproc shape, defaulting to 1 -- casting NaN
            # straight to int yields a garbage sentinel (INT_MIN) that corrupts
            # sumK and stalls the transition (never fires -> unbounded token
            # accumulation -> OOM). Mirrors ctmc_ssg.py and ssa/serial.py.
            fp = nparam.get('firingphases', None) if isinstance(nparam, dict) else getattr(nparam, 'firingphases', None)
            firing_proc = nparam.get('firingproc', None) if isinstance(nparam, dict) else getattr(nparam, 'firingproc', None)
            fphases = np.atleast_1d(np.asarray(fp, dtype=float)) if fp is not None else np.full(nmodes, np.nan)
            K = np.ones(nmodes, dtype=int)
            for _m in range(nmodes):
                if _m < fphases.size and not np.isnan(fphases[_m]):
                    K[_m] = int(fphases[_m])
                elif firing_proc is not None and _m < len(firing_proc) and firing_proc[_m] is not None:
                    K[_m] = int(np.atleast_2d(np.asarray(firing_proc[_m][0])).shape[0])
            K = np.maximum(K, 1).astype(int)
            Ks = np.concatenate([[0], np.cumsum(K[:-1])]).astype(int) if len(K) > 1 else np.array([0], dtype=int)
            sumK = int(np.sum(K))

            space_buf = inspace[:, :nmodes]
            space_srv = inspace[:, nmodes:nmodes + sumK]
            # ``fired`` counters live between ``phases`` and ``var``. The
            # PHASE handler must propagate them unchanged so the rebuilt
            # output state matches sn.space[isf] for the hash lookup.
            fired_start = nmodes + sumK
            fired_end = fired_start + nmodes
            if fired_end <= n_cols - V:
                space_fired = inspace[:, fired_start:fired_end]
            else:
                space_fired = np.zeros((n_rows, nmodes))

            from .after_event_transition import after_event_transition
            return _untagged(after_event_transition(
                sn, ind, event, job_class, inspace,
                K, Ks, space_buf, space_srv, space_var,
                space_fired=space_fired), R)

        else:
            # Router, Cache, other: [srv(R) | var]
            space_buf = np.zeros((n_rows, 0))
            srv_start = max(0, n_cols - R - V)
            srv_end = n_cols - V
            if srv_start >= 0 and srv_end > srv_start:
                space_srv = inspace[:, srv_start:srv_end]
            else:
                space_srv = np.zeros((n_rows, R))

            if nt_val == router_val:
                from .after_event_router import after_event_router
                return _untagged(after_event_router(
                    sn, ind, event, job_class,
                    space_buf, space_srv, space_var), R)

            elif nt_val == cache_val:
                from .after_event_cache import after_event_cache
                return _untagged(after_event_cache(
                    sn, ind, event, job_class, R,
                    space_buf, space_srv, space_var, is_simulation), R)

    # Stateless node: no state change
    return _no_successor(R)


def _no_successor(R):
    """The no-successor result, with tag matrices of the right width."""
    return (np.zeros((0, 0)), np.zeros((0, 0)), np.ones((1, 1)),
            np.zeros((0, R)), np.zeros((0, R)))


def _untagged(res, R):
    """Widen a three-field handler result with all-zero tags: the node types
    that reach this path (Router, Fork, Cache, Transition) hold no server, so
    they can neither start nor preempt a service."""
    outspace, outrate, outprob = res[0], res[1], res[2]
    n = np.atleast_2d(outspace).shape[0] if np.asarray(outspace).size > 0 else 0
    return outspace, outrate, outprob, np.zeros((n, R)), np.zeros((n, R))


def after_event_hashed(sn, ind, inhash, event, job_class, hash_maps=None):
    """
    Hash-based wrapper for afterEvent.

    Looks up the state from the hash index, calls afterEvent, then
    converts output states back to hash indices.

    Port from JAR State.java:afterEventHashed().

    Args:
        sn: NetworkStruct
        ind: Node index (0-based)
        inhash: Input state hash (row index in sn.space[isf])
        event: EventType
        job_class: Job class index (0-based)
        hash_maps: Optional precomputed hash maps

    Returns:
        Tuple of (outhash, outrate, outprob):
        - outhash: np.ndarray of output hash indices (-1 if not found)
        - outrate: np.ndarray of rates
        - outprob: np.ndarray of probabilities
    """
    R = int(getattr(sn, 'nclasses', 0) or 0)
    empty = (np.array([-1]), np.array([0.0]), np.array([0.0]),
             np.zeros((0, R)), np.zeros((0, R)))
    if inhash < 0:
        return empty

    isf = int(sn.nodeToStateful[ind])
    if isf not in sn.space or sn.space[isf] is None:
        return empty

    space = np.atleast_2d(sn.space[isf])
    if inhash >= space.shape[0]:
        return empty

    inspace = space[inhash:inhash + 1]

    outspace, outrate, outprob, outstart, outpreempt = after_event(
        sn, ind, inspace, event, job_class)

    if outspace.size == 0:
        return empty

    outhash = get_hash(sn, ind, outspace, hash_maps)
    outrate = outrate.ravel() if outrate.size > 0 else np.array([0.0])
    outprob = outprob.ravel() if outprob.size > 0 else np.array([0.0])

    # the tags travel alongside the hashed successor: both CTMC and SSA reach
    # the state machine through here
    return outhash, outrate, outprob, outstart, outpreempt
