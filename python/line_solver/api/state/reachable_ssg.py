"""
Reachability-based CTMC state-space generator for native fork-join models.

Port of matlab/src/lang/+State/reachableSpaceGenerator.m (restricted to the
subset needed by the FJ path). Fork firings break per-chain population
conservation, so the population-lattice enumeration in ctmc_ssg cannot
generate FJ state spaces; this BFS explores only states reachable from the
initial state using the sync actions and the fork firing synchronizations
(sn.fjsync). It returns the same tuple shape as ctmc_ssg so the downstream
generator/metrics path is unchanged.
"""

import numpy as np
from ...constants import EventType


def _pad_left(row, width):
    row = np.asarray(row, dtype=float).ravel()
    if row.shape[0] < width:
        return np.concatenate([np.zeros(width - row.shape[0]), row])
    return row


class _NodeSpace:
    """Growable per-node raw state space with row lookup (left zero-pad)."""

    def __init__(self, init_row):
        init_row = np.asarray(init_row, dtype=float).ravel()
        self.width = init_row.shape[0]
        self.rows = [init_row]
        self._index = {tuple(init_row): 0}

    def find_or_add(self, row):
        row = np.asarray(row, dtype=float).ravel()
        if row.shape[0] > self.width:
            # widen existing rows (left zero-pad) to the new width
            self.width = row.shape[0]
            self.rows = [_pad_left(r, self.width) for r in self.rows]
            self._index = {tuple(r): i for i, r in enumerate(self.rows)}
        elif row.shape[0] < self.width:
            row = _pad_left(row, self.width)
        key = tuple(row)
        idx = self._index.get(key, -1)
        if idx >= 0:
            return idx
        self.rows.append(row)
        idx = len(self.rows) - 1
        self._index[key] = idx
        return idx

    def get(self, idx):
        return self.rows[idx]

    def matrix(self):
        return np.array([_pad_left(r, self.width) for r in self.rows], dtype=float)


def _widen_ordered_buffer(sn, isf, row):
    """Left-pad an ordered-buffer station's seed row to its capacity width.

    `sn.state[isf]` is sized from the INITIAL MARKING: `fromMarginal` gives an
    ordered buffer `max(1, n - S)` waiting slots, so an empty FCFS station comes
    back with exactly ONE. The lattice generator never uses that width -- it
    enumerates every marginal up to the capacity and the widest row sets the
    space -- but the reachable walk has no such pass, and `after_event_station`
    REFUSES an arrival that finds no free slot instead of appending one. That
    refusal is deliberate: under the lattice a full row means the state-space
    cutoff, and the state above the cutoff is absent rather than blocked.

    Seeded at the marginal width, an ordered-buffer station therefore cannot hold
    more than `S + 1` jobs for the whole walk. On a fork-join model the sibling
    arrival then comes back empty and `after_fj_event` drops the entire firing,
    leaving a Fork-occupied vanishing state with no immediate outgoing arc. The
    reference reaches the full space by growing the buffer inside the arrival
    handler (`State.afterEventStation`, "this section dynamically grows the
    number of elements in the buffer"); widening the seed instead keeps this
    port's stricter width rule and costs one zero-padding per station.

    The buffer is right-aligned and leads the row, so padding with zeros on the
    left adds empty slots and changes nothing else. A station whose capacity is
    not finite is left alone: there the walk is bounded by the cutoff, not by the
    capacity, and the cutoff is applied to whole states rather than to one row.
    """
    from ...lang.base import SCHED_BUFFER_CLASS_TAG, SCHED_PREEMPT

    stf2st = getattr(sn, 'statefulToStation', None)
    if stf2st is None:
        return row
    ist = int(np.asarray(stf2st).ravel()[isf])
    if ist < 0:
        return row  # a Fork, Join, Router or Cache carries no buffer/server split
    sched = sn.sched[ist] if getattr(sn, 'sched', None) is not None else None
    if sched in SCHED_BUFFER_CLASS_TAG:
        percol = 1
    elif sched in SCHED_PREEMPT:
        percol = 2  # the preempt families store [class, phase] PAIRS
    else:
        return row  # per-class-count buffers are R wide whatever the occupancy

    S = float(np.asarray(sn.nservers).ravel()[ist])
    cap = float(np.asarray(sn.cap).ravel()[ist])
    if not np.isfinite(S) or not np.isfinite(cap):
        return row
    srvw = int(np.sum(np.asarray(sn.phasessz[ist], dtype=int)))
    nvar = int(np.sum(sn.nvars[int(np.asarray(sn.statefulToNode).ravel()[isf])])) \
        if getattr(sn, 'nvars', None) is not None else 0
    bufw = int(row.shape[0]) - srvw - nvar
    if bufw <= 0:
        return row
    want = max(1, int(cap) - int(S)) * percol
    if want <= bufw:
        return row
    return np.concatenate([np.zeros(want - bufw), row])


def _initial_raw_states(sn):
    """Build the initial raw local state of every stateful node.

    `sn.state[isf]` ALREADY IS that raw row: `Network.init_default` builds it
    through `initFromMarginal`, so the encoding of each node -- per-class counts
    at a Delay, counts plus service positions at an FCFS station, the per-class
    parent and buffered vectors at a Fork and a Join -- is the one the successor
    functions expect. It is seeded verbatim, as MATLAB's
    `State.reachableSpaceGenerator` does (`space{i} = sn.state{i}`).

    Passing that row back through `fromMarginal` instead, as this did, treats it
    as a per-class marginal. That is only the same vector at a node whose raw
    encoding IS the marginal. At an FCFS station it is longer than the class
    count, and `fromMarginal` then compares a state vector of one width against a
    `classcap` row of another: `fj_tiny_closed` died there with "operands could
    not be broadcast together with shapes (4,) (3,)".
    """
    R = int(sn.nclasses)
    nstateful = int(sn.nstateful)

    init = [None] * nstateful
    for isf in range(nstateful):
        if sn.state is not None and isf < len(sn.state) and sn.state[isf] is not None:
            row = np.asarray(sn.state[isf], dtype=float).ravel()
        else:
            row = np.zeros(R)
        init[isf] = _widen_ordered_buffer(sn, isf, row.copy())
    return init


def reachable_ssg(sn, options=None):
    """
    Generate the reachable CTMC state space for an FJ-augmented struct.

    Returns (state_space, state_space_aggr, state_space_hashed, sn) with
    sn.space populated (dict isf -> ndarray of per-node raw rows), matching
    the ctmc_ssg contract used by solver_ctmc_basic.
    """
    from .after_event import after_event
    from .after_fj_event import after_fj_event
    from ...lang.sync import refresh_sync
    from .ctmc_ssg import _build_state_space_aggr

    nstateful = int(sn.nstateful)
    nnodes = int(sn.nnodes)
    local = nnodes  # LOCAL passive sentinel

    sync = refresh_sync(sn)
    A = len(sync)
    fjsync = sn.fjsync if getattr(sn, 'fjsync', None) else []
    FJ = len(fjsync)

    # per-node growable spaces, seeded with the initial raw state
    init_rows = _initial_raw_states(sn)
    spaces = [_NodeSpace(init_rows[isf]) for isf in range(nstateful)]

    init_hashed = tuple(0 for _ in range(nstateful))
    SSh = [init_hashed]
    seen = {init_hashed}
    stack = [init_hashed]

    def _glspace(hashed):
        return [spaces[isf].get(int(hashed[isf])) for isf in range(nstateful)]

    def _hash_successor(cur_hashed, changes):
        """changes: dict isf -> new raw row. Returns new hashed tuple."""
        new_h = list(cur_hashed)
        for isf, row in changes.items():
            new_h[isf] = spaces[isf].find_or_add(row)
        return tuple(new_h)

    def _push(new_h):
        if new_h not in seen:
            seen.add(new_h)
            SSh.append(new_h)
            stack.append(new_h)

    while stack:
        h = stack.pop()
        glspace = _glspace(h)

        for a in range(A):
            act = sync[a]
            node_a = int(act.active.node)
            class_a = int(act.active.job_class)
            event_a = act.active.event
            node_p = int(act.passive.node)
            class_p = int(act.passive.job_class)
            event_p = act.passive.event

            if not sn.isstateful[node_a]:
                continue
            isf_a = int(sn.nodeToStateful[node_a])
            is_local = (node_p >= nnodes)

            out_a, rate_a, _, _, _ = after_event(sn, node_a, np.atleast_2d(glspace[isf_a]),
                                           event_a, class_a)
            out_a = np.atleast_2d(out_a)
            if out_a.size == 0:
                continue
            rate_a = np.asarray(rate_a).ravel()

            for ia in range(out_a.shape[0]):
                if ia < len(rate_a) and rate_a[ia] == 0:
                    continue
                new_a_row = out_a[ia, :]
                if is_local:
                    _push(_hash_successor(h, {isf_a: new_a_row}))
                    continue
                if not sn.isstateful[node_p]:
                    continue
                isf_p = int(sn.nodeToStateful[node_p])
                state_p_row = new_a_row if node_p == node_a else glspace[isf_p]
                out_p, rate_p, _, _, _ = after_event(sn, node_p, np.atleast_2d(state_p_row),
                                               event_p, class_p)
                out_p = np.atleast_2d(out_p)
                if out_p.size == 0:
                    continue
                for ip in range(out_p.shape[0]):
                    changes = {isf_a: new_a_row, isf_p: out_p[ip, :]}
                    _push(_hash_successor(h, changes))

        # fork firing synchronizations (multi-node atomic successors)
        for k in range(FJ):
            fj_states, fj_rate, fj_prob = after_fj_event(sn, fjsync[k], glspace, False)
            for io in range(len(fj_states)):
                if io < len(fj_prob) and fj_prob[io] <= 0:
                    continue
                gl_io = fj_states[io]
                changes = {}
                for isf in range(nstateful):
                    new_row = np.asarray(gl_io[isf], dtype=float).ravel()
                    if not np.array_equal(_pad_left(new_row, spaces[isf].width),
                                          _pad_left(glspace[isf], spaces[isf].width)):
                        changes[isf] = new_row
                _push(_hash_successor(h, changes))

    # finalize sn.space
    sn.space = {isf: spaces[isf].matrix() for isf in range(nstateful)}

    # build hashed matrix and concatenated raw state space
    n_states = len(SSh)
    state_space_hashed = np.array([list(hh) for hh in SSh], dtype=int)
    widths = [spaces[isf].width for isf in range(nstateful)]
    state_space = np.zeros((n_states, int(np.sum(widths))))
    for s in range(n_states):
        col = 0
        for isf in range(nstateful):
            w = widths[isf]
            state_space[s, col:col + w] = spaces[isf].get(int(state_space_hashed[s, isf]))
            col += w

    state_space_aggr = _build_state_space_aggr(sn, state_space_hashed)
    return state_space, state_space_aggr, state_space_hashed, sn
