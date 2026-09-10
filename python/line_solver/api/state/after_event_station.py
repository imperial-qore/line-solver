"""
Station node event handler for afterEvent dispatch.

Handles ARV, DEP, and PHASE events at station nodes (queues, delays, sources).
This is the largest handler, covering all scheduling strategies.

Port from JAR AfterEventStation.java.
"""

import numpy as np
from ...constants import EventType, GlobalConstants, ProcessType
from ...lang.base import (SchedStrategy, NodeType,
                          SCHED_PREEMPT as _SCHED_PREEMPT,
                          SCHED_PREEMPT_LCFS as _SCHED_PREEMPT_LCFS,
                          SCHED_PREEMPT_PRIO as _SCHED_PREEMPT_PRIO,
                          SCHED_PREEMPT_RESUME as _SCHED_PREEMPT_RESUME)
from .reply_block import (reply_width as _reply_width, reply_blocked as _reply_blocked,
                          reply_block_info as _reply_block_info, is_reply_class as _is_reply_class)


def after_event_station(sn, ind, inspace, event, job_class,
                        M, R, ist, K, Ks, hasOnlyExp,
                        mu, phi, pie, proc, ismkvmodclass,
                        lldscaling, lldlimit, cdscaling,
                        capacity, classcap, V,
                        space_buf, space_srv, space_var,
                        is_simulation=False, no_promote=False):
    """
    Handle events at a station node.

    Args:
        sn: NetworkStruct
        ind: Node index (0-based)
        inspace: Full input state, shape (n_rows, n_cols)
        event: EventType (ARV, DEP, or PHASE)
        job_class: Job class index (0-based)
        M: Number of stations
        R: Number of classes
        ist: Station index
        K: Phase counts per class, shape (R,)
        Ks: Phase shifts per class, shape (R,)
        hasOnlyExp: True if all processes are exponential (single-phase)
        mu: Service rates, mu[ist][r][k]
        phi: Completion probabilities, phi[ist][r][k]
        pie: Entry probabilities, pie[ist][r][k]
        proc: Process matrices, proc[ist][r] = [D0, D1]
        ismkvmodclass: MAP/MMPP2 flags per class, shape (R,)
        lldscaling: Load-level-dependent scaling, shape (M, lldlimit)
        lldlimit: Max load level for scaling
        cdscaling: Class-dependent scaling functions
        capacity: Station capacities, shape (M,)
        classcap: Per-class capacities, shape (M, R)
        V: Number of state variables
        space_buf: Buffer state
        space_srv: Server state
        space_var: Variable state

    Returns:
        Tuple of (outspace, outrate, outprob)
    """
    # INF/Delay nservers kept as float infinity (never int()); every downstream use (min, ni<S, ni<=S) is correct with inf.
    S = sn.nservers[ist]
    S = int(S) if np.isfinite(S) else float('inf')
    sched = sn.sched[ist]
    K = np.atleast_1d(K).astype(int)
    Ks = np.atleast_1d(Ks).astype(int)

    # breakdown status = trailing local-var col (0=down,1=up), exclusive with BAS marker/polling controller; down: no DEP/PHASE unless degraded rate set.
    is_breakdown_station = False
    if getattr(sn, 'hasbreakdown', None) is not None:
        from .ctmc_ssg import _is_breakdown_station
        is_breakdown_station = _is_breakdown_station(sn, ind)
    down_rate_scale = 1.0
    _in = np.atleast_2d(inspace)
    if is_breakdown_station and _in.size > 0 and \
            event in (EventType.DEP, EventType.PHASE):
        if _in[0, -1] == 0:  # server down
            down_rate = 0.0
            dsr = getattr(sn, 'downServiceRates', None)
            if dsr is not None:
                dsr = np.atleast_2d(np.asarray(dsr, dtype=float))
                if ist < dsr.shape[0] and job_class < dsr.shape[1]:
                    down_rate = float(dsr[ist, job_class])
            if down_rate <= 0:
                return _empty_result(sn)
            up_rate = float(np.atleast_2d(sn.rates)[ist, job_class])
            if not np.isfinite(up_rate) or up_rate <= 0:
                raise ValueError(
                    "Station '%s' declares a down-server service rate for class '%s' but "
                    "has no finite up-server service rate to rescale."
                    % (sn.nodenames[ind], sn.classnames[job_class]))
            down_rate_scale = down_rate / up_rate

    if event == EventType.FAILURE:
        # up server fails at rate breakdownMu regardless of serving state; only status column changes (memoryless service resumes on repair, no job moves).
        if is_breakdown_station and _in.size > 0 and _in[0, -1] == 1:
            outspace = _in.astype(float).copy()
            outspace[:, -1] = 0
            # only the server status changes: no job moves, so no tag
            return (outspace,
                    np.full((outspace.shape[0], 1), float(sn.breakdownMu[ist])),
                    np.ones((outspace.shape[0], 1)),
                    np.zeros((outspace.shape[0], R)), np.zeros((outspace.shape[0], R)))
        return _empty_result(sn)
    elif event == EventType.REPAIR:
        # A down server is restored at the memoryless rate repairMu.
        if is_breakdown_station and _in.size > 0 and _in[0, -1] == 0:
            outspace = _in.astype(float).copy()
            outspace[:, -1] = 1
            # REPAIR emits no START: the supported breakdown model RESUMES the
            # job the server was holding rather than restarting it
            return (outspace,
                    np.full((outspace.shape[0], 1), float(sn.repairMu[ist])),
                    np.ones((outspace.shape[0], 1)),
                    np.zeros((outspace.shape[0], R)), np.zeros((outspace.shape[0], R)))
        return _empty_result(sn)

    outspace, outrate, outprob, outstart, outpreempt = _dispatch_station_event(
        sn, ind, inspace, event, job_class, M, R, ist, S, K, Ks, hasOnlyExp,
        mu, phi, pie, proc, ismkvmodclass, lldscaling, lldlimit, cdscaling,
        capacity, classcap, V, sched, space_buf, space_srv, space_var,
        is_simulation, no_promote)

    # degraded down-server rate rescales the up-server-computed completion rate; down_rate_scale=1 (no-op) when server up or no degraded rate configured.
    if down_rate_scale != 1.0 and np.asarray(outrate).size > 0:
        outrate = np.asarray(outrate, dtype=float) * down_rate_scale
    # the tags travel with the successors: only the rate is rescaled
    return outspace, outrate, outprob, outstart, outpreempt


def _dispatch_station_event(sn, ind, inspace, event, job_class, M, R, ist, S,
                            K, Ks, hasOnlyExp, mu, phi, pie, proc,
                            ismkvmodclass, lldscaling, lldlimit, cdscaling,
                            capacity, classcap, V, sched, space_buf, space_srv,
                            space_var, is_simulation, no_promote):
    """Route an event to its per-event handler. Split out of
    after_event_station so that the server-breakdown gating and the degraded
    down-server rescaling wrap every scheduling branch without any of them
    having to know about the server status."""
    if event == EventType.ARV:
        return _handle_arv(sn, ind, inspace, job_class, M, R, ist, S, K, Ks,
                           hasOnlyExp, mu, phi, pie, proc, ismkvmodclass,
                           lldscaling, lldlimit, cdscaling, capacity, classcap,
                           V, sched, space_buf, space_srv, space_var, is_simulation)
    elif event == EventType.DEP:
        outspace, outrate, outprob, outstart, outpreempt = _handle_dep(
            sn, ind, inspace, job_class, M, R, ist, S, K, Ks,
            hasOnlyExp, mu, phi, pie, proc, ismkvmodclass,
            lldscaling, lldlimit, cdscaling, capacity, classcap,
            V, sched, space_buf, space_srv, space_var, no_promote)
        # true-BAS transfer of already-blocked job at rate 1e7 clears successor's blocked marker; complementary become-blocked edge from CTMC/SSA generator.
        _in = np.atleast_2d(inspace)
        if (_is_bas_station_marker(sn, ist) and _in.size > 0 and _in.shape[1] > 0
                and _in[0, -1] == 1 and np.asarray(outspace).size > 0):
            outspace = np.atleast_2d(outspace).astype(float).copy()
            outspace[:, -1] = 0
            outrate = np.full((outspace.shape[0], 1), 1.0e7)
        return outspace, outrate, outprob, outstart, outpreempt
    elif event == EventType.PHASE:
        return _handle_phase(sn, ind, inspace, job_class, M, R, ist, S, K, Ks,
                             hasOnlyExp, mu, phi, pie, proc, ismkvmodclass,
                             lldscaling, lldlimit, cdscaling,
                             V, sched, space_buf, space_srv, space_var)
    elif event == EventType.RETRY:
        return _handle_retry(sn, ind, inspace, job_class, ist, R, S, K, Ks,
                             pie, space_buf, space_srv, space_var, is_simulation)
    elif event == EventType.RENEGE:
        return _handle_renege(sn, ind, inspace, job_class, ist, R, K, Ks,
                              space_buf, space_srv, space_var)
    elif event == EventType.SWITCH:
        return _handle_switch(sn, ind, inspace, job_class, ist, R, K, Ks, pie,
                              space_buf, space_srv, space_var)

    return _empty_result(sn)


def _has_blocking_droprule(sn, ist, job_class):
    """True when class job_class at station ist uses a blocking drop rule.

    BAS, BBS and RSRD hold a refused arrival on the upstream server, so the
    arrival is not an enabled event and its rows are removed. DROP and WAITQ
    instead leave the unchanged state as a valid successor, which is the drop.
    """
    droprule = getattr(sn, 'droprule', None)
    if droprule is None:
        return False
    try:
        from ..sn.network_struct import DropStrategy as SnDropStrategy
        rule = int(droprule[ist, job_class])
    except (IndexError, TypeError, ValueError):
        return False
    return rule in (int(SnDropStrategy.BAS), int(SnDropStrategy.BBS),
                    int(SnDropStrategy.RSRD))


def _arrival_is_lost(sn, ist, job_class):
    """True when an arrival of job_class that finds no room at station ist is LOST,
    False when it must BLOCK the upstream instead. Mirrors MATLAB State.arrivalIsLost
    and JAR State.arrivalIsLost.

    The rule is the CLASS TYPE, not the drop rule:
      OPEN class   -> LOST. The memoryless external stream simply does not enter; the
                      caller leaves the state UNCHANGED (self-loop) so the arrival
                      event still fires and the arrival-rate statistic counts the
                      OFFERED job. A self-loop cancels on the generator diagonal, so
                      QLen/Util/Tput are unaffected.
      CLOSED class -> BLOCKED. A closed network's N jobs cannot be dropped (population
                      conservation is an invariant); the caller returns an EMPTY row
                      so the upstream departure is disabled until room frees, which is
                      what the true-BAS become-blocked edge tests for.
    An explicit blocking rule (BAS) also asks for blocking. Same open/closed predicate
    as the CTMC analyzer's canDropClass and the BUG-12 utilization guard; the
    conventions are complementary (arrival rate = offered, Util/QLen/Tput = carried).
    """
    if _has_blocking_droprule(sn, ist, job_class):
        return False
    # upstream-declared BAS records the destination side separately in sn.isbasdestination; see _kb/04-networkstruct.md isbasblocking field and BUG-83.
    isbd = getattr(sn, 'isbasdestination', None)
    if isbd is not None:
        arr2 = np.atleast_2d(np.asarray(isbd))
        if ist < arr2.shape[0] and job_class < arr2.shape[1] and bool(arr2[ist, job_class]):
            return False  # refusing here must block the upstream BAS station
    njobs = getattr(sn, 'njobs', None)
    if njobs is None:
        return True
    arr = np.asarray(njobs).flatten()
    if job_class >= arr.size:
        return True
    return bool(np.isinf(arr[job_class]))


def _reply_held_rows(sn, ind, space_var):
    """Servers held at node ind for pending synchronous calls, per state row.

    Returns the scalar 0.0 for every model without REPLY signals at this node,
    so callers can subtract it from the server count unconditionally and pay
    nothing when the feature is unused.
    """
    if _reply_width(sn, ind) == 0:
        return 0.0
    return _reply_blocked(sn, ind, space_var)[1]


def _reply_call_slot(sn, ind, job_class):
    """Local-variable index of the blocked-server counter that a departing
    job_class job at node ind increments, or -1 when this departure is not a
    synchronous call from this station."""
    rb = getattr(sn, 'replyblock', None)
    if rb is None:
        return -1
    rb = np.atleast_2d(np.asarray(rb))
    if rb.size == 0 or ind >= rb.shape[0] or job_class >= rb.shape[1]:
        return -1
    if rb[ind, job_class] <= 0:
        return -1
    return int(_reply_block_info(sn, ind).slot[job_class])


def _is_bas_station_marker(sn, ist):
    """True when station ist declares a blocked marker as the trailing state column.

    This must agree exactly with the predicate that DECLARES the column when the
    per-node state space is built (ctmc_ssg._is_bas_station: BAS on some class at a
    finite-capacity station). Gating on a broader "is blocked by some rule" test
    instead desynchronizes the reader from the declaration: BBS/RSRD stations never
    get a marker column, so the trailing column is their last server phase slot, and
    the transfer post-process below would clobber a job in service and fire the DEP at
    the 1e7 immediate rate. MATLAB and the JAR get this for free by gating directly on
    sn.nvars(ind,2R+1)==1 / sn.nvars(ind,2R)==1; the native Python sn.nvars is an
    accumulated per-node width rather than a column-addressed layout, so the shared
    predicate is what keeps reader and declarer in sync.
    """
    from .ctmc_ssg import _is_bas_station
    return _is_bas_station(sn, ist)


# ---------------------------------------------------------------------------
# ARV handler
# ---------------------------------------------------------------------------

def _handle_arv(sn, ind, inspace, job_class, M, R, ist, S, K, Ks,
                hasOnlyExp, mu, phi, pie, proc, ismkvmodclass,
                lldscaling, lldlimit, cdscaling, capacity, classcap,
                V, sched, space_buf, space_srv, space_var, is_simulation=False):
    """Handle arrival events."""
    from .marginal import toMarginalAggr
    # Place: no servers/phases; token only increments class marking, else row exceeds marking state space, drop it. _kb/11-conventions-and-gotchas.md SPN.
    if sn.nodetype[ind] == NodeType.PLACE:
        out = np.atleast_2d(np.array(inspace, dtype=float).copy())
        out[:, job_class] += 1
        n_out = out.shape[0]
        # passive action: the rate is set by the active node. FIVE values, as
        # every other branch returns since the start/preempt tags were added:
        # a Place holds no server, so nothing starts or is preempted there.
        return (out,
                -1.0 * np.ones(n_out).reshape(-1, 1),
                np.ones(n_out).reshape(-1, 1),
                np.zeros((n_out, R)), np.zeros((n_out, R)))
    # signal/catastrophe arrival removes job(s) already present and is annihilated itself; mirrors MATLAB State.afterEventStationSignal.
    if getattr(sn, 'issignal', None) is not None and bool(sn.issignal[job_class]):
        # A REPLY signal is not a negative customer: it completes a synchronous
        # call, releasing the server this station holds for the caller, and then
        # joins as an ordinary job. Only stations that actually hold a block for
        # it take this path; elsewhere a REPLY class is a plain job class and
        # falls through to the normal arrival handling below.
        if _is_reply_class(sn, job_class) and _reply_block_info(sn, ind).width > 0:
            return _handle_arv_reply(sn, ind, ist, job_class, K, Ks, S, pie,
                                     space_buf, space_srv, space_var,
                                     inspace, is_simulation)
        from .signal_removal import handle_signal_arrival
        return handle_signal_arrival(sn, ind, ist, inspace, job_class, sched,
                                     K, Ks, S, space_buf, space_srv, space_var,
                                     is_simulation)
    # ordinary SPN Place arrival: see _kb/11-conventions-and-gotchas.md (An ordinary Place must be special-cased before the generic scheduling switch).
    nodetype = getattr(sn, 'nodetype', None)
    if nodetype is not None and int(np.asarray(nodetype).ravel()[ind]) == int(NodeType.PLACE):
        is_queueing = False
        _iqp = getattr(sn, 'isqueueingplace', None)
        if _iqp is not None:
            _iqp = np.asarray(_iqp).ravel()
            if ist < _iqp.size and bool(_iqp[ist]):
                is_queueing = True
        if not is_queueing:
            rows = np.atleast_2d(inspace).astype(float)
            outspace = rows.copy()
            cap = np.inf
            if classcap is not None:
                cap = float(np.asarray(classcap)[ist, job_class])
            out_states, out_rates, out_probs = [], [], []
            for i in range(rows.shape[0]):
                new_row = rows[i].copy()
                if new_row[job_class] < cap:
                    new_row[job_class] += 1
                    out_states.append(new_row)
                    out_rates.append(-1.0)   # passive: rate set by active source
                    out_probs.append(1.0)
                else:
                    # place full: arrival is blocked and lost (no state change)
                    out_states.append(rows[i].copy())
                    out_rates.append(-1.0)
                    out_probs.append(0.0)
            return _finalize(out_states, out_rates, out_probs, inspace, None, None, R)
    ni, nir = toMarginalAggr(sn, ind, inspace, K, Ks, space_buf, space_srv, space_var)

    # Get phase entry probabilities
    pentry = _get_pie(pie, ist, job_class, K)

    out_states = []
    out_rates = []
    out_probs = []
    # START/PREEMPT tags of the arrival arcs, keyed by the index of the
    # successor they annotate: only the arcs that take or free a server say
    # anything, every other successor is a zero row. See _finalize.
    start_tags = {}
    preempt_tags = {}

    for kentry in range(int(K[job_class])):
        space_buf_k = space_buf.copy()
        space_srv_k = space_srv.copy()
        space_var_k = space_var.copy() if space_var.size > 0 else space_var

        n_rows = space_srv_k.shape[0]
        # per-row class that takes a server in this row (0 = none); read at the
        # common composition loop below, where the successors are appended
        start_row = np.zeros(n_rows, dtype=int)
        outprob_k = np.full(n_rows, pentry[kentry])
        valid = np.ones(n_rows, dtype=bool)

        # MAP job entering service resumes its persistent server-phase entry with probability 1; mirrors MATLAB afterEventStation.m:53-56.
        if (ismkvmodclass is not None and job_class < len(ismkvmodclass)
                and ismkvmodclass[job_class] and space_var.size > 0):
            _mvcol = int(np.sum(np.asarray(ismkvmodclass)[:job_class]))
            _sv = np.atleast_2d(space_var)
            if _mvcol < _sv.shape[1]:
                for i in range(n_rows):
                    if int(_sv[i if _sv.shape[0] > 1 else 0, _mvcol]) == kentry:
                        outprob_k[i] = 1.0
                    else:
                        valid[i] = False

        if sched == SchedStrategy.EXT:
            # Source: virtual arrival, rate=0, state unchanged
            outrate_k = np.zeros(n_rows)
            # Only one phase iteration needed for EXT
            for i in range(n_rows):
                if valid[i]:
                    state_row = _compose_state(space_buf_k[i] if space_buf_k.ndim >= 2 else space_buf_k,
                                               space_srv_k[i] if space_srv_k.ndim >= 2 else space_srv_k,
                                               space_var_k[i] if space_var_k.ndim >= 2 and space_var_k.size > 0 else space_var_k)
                    out_states.append(state_row)
                    out_rates.append(-1.0)
                    out_probs.append(1.0)
            return _finalize(out_states, out_rates, out_probs, inspace, None, None, R)

        elif sched in (SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.DPS,
                        SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO,
                        SchedStrategy.GPSPRIO, SchedStrategy.LPS):
            # Jobs enter service immediately
            for i in range(n_rows):
                if classcap is not None and nir[i, job_class] >= classcap[ist, job_class]:
                    # no room: DROP leaves state unchanged (job destroyed, upstream departure stays on); BAS/BBS/RSRD remove the row, holding back upstream departure.
                    if _has_blocking_droprule(sn, ist, job_class):
                        valid[i] = False
                    continue
                space_srv_k[i, int(Ks[job_class]) + kentry] += 1
                start_row[i] = job_class + 1  # the job enters service at once

        elif sched == SchedStrategy.POLLING:
            # polling controller (not the arrival) picks who is served; only a parked server (empty station, immediate switchovers) visits the arrival at once.
            from .polling import polling_info, polling_get, polling_set, polling_next
            pinfo_a = polling_info(sn, ind)
            for i in range(n_rows):
                srvclass_a = -1
                for ra in range(R):
                    if float(np.sum(space_srv_k[i, int(Ks[ra]):int(Ks[ra]) + int(K[ra])])) > 0:
                        srvclass_a = ra
                        break
                var_row = space_var_k[i] if space_var_k.ndim >= 2 else space_var_k
                _, swk_a, _ = polling_get(pinfo_a, var_row, srvclass_a)
                if srvclass_a < 0 and swk_a == 0:
                    # parked controller: the arriving job is the only work present, so it enters service directly in phase kentry rather than the buffer.
                    nbuf_a = space_buf_k[i, :R].copy()
                    nbuf_a[job_class] += 1
                    q_a, _, budget_a = polling_next(pinfo_a, job_class, nbuf_a, R, arrived=True)
                    space_srv_k[i, int(Ks[job_class]) + kentry] += 1
                    start_row[i] = job_class + 1  # a parked server takes it at once
                    space_var_k[i] = polling_set(pinfo_a, var_row, q_a, 0, budget_a)
                else:
                    space_buf_k[i, job_class] += 1

        elif sched in (SchedStrategy.SIRO, SchedStrategy.SEPT, SchedStrategy.LEPT):
            # SIRO/SEPT/LEPT share unordered per-class-count buffer; idle test on server occupancy (not total count), so a self-loop re-enters the vacated server.
            for i in range(n_rows):
                srv_count = float(np.sum(space_srv_k[i])) if space_srv_k.ndim >= 2 else float(np.sum(space_srv_k))
                if srv_count < S:
                    space_srv_k[i, int(Ks[job_class]) + kentry] += 1
                    start_row[i] = job_class + 1
                else:
                    if space_buf_k.ndim >= 2 and job_class < space_buf_k.shape[1]:
                        space_buf_k[i, job_class] += 1
                    elif space_buf_k.ndim == 1 and job_class < len(space_buf_k):
                        space_buf_k[job_class] += 1

        elif sched in (SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS):
            # Servers held by synchronous calls awaiting a REPLY are NOT
            # available to an arriving job: subtract them from the server count.
            # Computed once for every row (a zero scalar, hence free, for every
            # model without reply signals).
            Seff_rows = S - _reply_held_rows(sn, ind, space_var_k)
            for i in range(n_rows):
                ni_val = ni[i]
                # idle: server occupancy not total count; immediate-feedback self-loop (idle server, nonempty buffer) re-enters vacated server; afterEventStation.m.
                srv_count = float(np.sum(space_srv_k[i])) if space_srv_k.ndim >= 2 else float(np.sum(space_srv_k))
                Seff = Seff_rows if np.isscalar(Seff_rows) else Seff_rows[min(i, len(Seff_rows) - 1)]
                if srv_count < Seff:
                    # Idle server available: enter service directly
                    space_srv_k[i, int(Ks[job_class]) + kentry] += 1
                    start_row[i] = job_class + 1
                else:
                    # All servers busy: add to buffer
                    cap_ok = True
                    if capacity is not None and ni_val >= capacity[ist]:
                        cap_ok = False
                    if classcap is not None and nir[i, job_class] >= classcap[ist, job_class]:
                        cap_ok = False
                    if not cap_ok:
                        # refusal LOST (self-loop, open), BLOCKED (invalid row, closed/BAS), by class type not drop rule; _kb/11-conventions-and-gotchas.md (BUG-81/BUG-84).
                        if not _arrival_is_lost(sn, ist, job_class):
                            valid[i] = False
                        continue

                    # Insert into buffer (right-aligned: jobs packed to the right)
                    buf_row = space_buf_k[i] if space_buf_k.ndim >= 2 else space_buf_k
                    # Find leftmost non-zero (first existing job in buffer)
                    first_job = -1
                    for b in range(len(buf_row)):
                        if buf_row[b] > 0:
                            first_job = b
                            break
                    if first_job > 0:
                        # Insert just left of first existing job
                        buf_row[first_job - 1] = job_class + 1
                    elif first_job == 0:
                        # buffer physically full but capacity allows more: grow by one slot on the left in simulation mode; mirrors MATLAB afterEventStation.m:73-99.
                        if is_simulation:
                            if space_buf_k.ndim >= 2:
                                space_buf_k = np.hstack([np.zeros((space_buf_k.shape[0], 1)), space_buf_k])
                                space_buf_k[i, 0] = job_class + 1
                            else:
                                space_buf_k = np.concatenate([[job_class + 1], space_buf_k])
                        else:
                            # CTMC width-full (state-space cutoff, not physical capacity) always BLOCKS even open, else offered not carried. _kb/11-conventions-and-gotchas.md
                            valid[i] = False
                            continue
                    elif len(buf_row) > 0:
                        # Empty buffer: insert at rightmost position
                        buf_row[len(buf_row) - 1] = job_class + 1
                    else:
                        # zero-width buffer: same treatment as width-full (grow in simulation mode, block in CTMC mode).
                        if is_simulation:
                            if space_buf_k.ndim >= 2:
                                space_buf_k = np.hstack([np.zeros((space_buf_k.shape[0], 1)), space_buf_k])
                                space_buf_k[i, 0] = job_class + 1
                            else:
                                space_buf_k = np.concatenate([[job_class + 1], space_buf_k])
                        else:
                            valid[i] = False
                            continue

        elif sched in _SCHED_PREEMPT:
            # preemptive: idle server enters directly; all-busy preempts one in-service job into buffer as a (class,phase) pair; output states composed inline.
            for i in range(n_rows):
                ni_val = ni[i]
                if ni_val < S:
                    # Idle server: enter service directly
                    space_srv_k[i, int(Ks[job_class]) + kentry] += 1
                    state_row = _compose_state(
                        space_buf_k[i] if space_buf_k.ndim >= 2 else space_buf_k,
                        space_srv_k[i] if space_srv_k.ndim >= 2 else space_srv_k,
                        space_var_k[i] if space_var_k.ndim >= 2 and space_var_k.size > 0 else space_var_k)
                    out_states.append(state_row)
                    out_rates.append(-1.0)
                    out_probs.append(pentry[kentry])
                    start_tags[len(out_states) - 1] = job_class + 1
                else:
                    # All servers busy: generate states for each possible preemption
                    srv_total = np.sum(space_srv_k[i])
                    any_preempt = False
                    classprio_arr = sn.classprio if hasattr(sn, 'classprio') and sn.classprio is not None else np.zeros(R)

                    for classpreempt in range(R):
                        is_prio_sched = sched in _SCHED_PREEMPT_PRIO
                        # priority-awareness is a property of the declared policy, never inferred from the data; see _kb/11-conventions-and-gotchas.md Modeling traps.
                        is_prio_aware = is_prio_sched
                        if is_prio_aware:
                            # Across priority groups a strictly higher-priority arrival
                            # preempts. WITHIN a group the base discipline decides: LCFS-PR
                            # keeps the NEWEST job in service, so an equal-priority arrival
                            # preempts, whereas FCFS-PR never lets an arrival preempt
                            # (MATLAB afterEventStation.m:377-392).
                            if sched in (SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO):
                                cannot_preempt = classprio_arr[job_class] > classprio_arr[classpreempt]
                            else:
                                cannot_preempt = classprio_arr[job_class] >= classprio_arr[classpreempt]
                            if cannot_preempt:
                                continue
                        for phasepreempt in range(int(K[classpreempt])):
                            col_preempt = int(Ks[classpreempt]) + phasepreempt
                            count_preempt = space_srv_k[i, col_preempt]
                            if count_preempt > 0:
                                srv_k = space_srv_k[i].copy()
                                buf_k = space_buf_k[i].copy() if space_buf_k.ndim >= 2 else space_buf_k.copy()
                                var_k = space_var_k[i].copy() if space_var_k.ndim >= 2 and space_var_k.size > 0 else (space_var_k.copy() if space_var_k.size > 0 else np.array([]))
                                # Remove preempted job from service
                                srv_k[col_preempt] -= 1
                                # Add arriving job to service
                                srv_k[int(Ks[job_class]) + kentry] += 1
                                # preempted (class,phase) pair fills rightmost empty slot, keeping it reachable in right-aligned enumerated space; mirrors MATLAB afterEventStation.m.
                                # PI discards the phase on resume, so it stores 1 rather than the
                                # phase reached; storing the real phase would split each state into
                                # K lumpable copies (MATLAB afterEventStation.m:433-437).
                                stored_phase = phasepreempt + 1 if sched in _SCHED_PREEMPT_RESUME else 1
                                for b in range(len(buf_k) - 2, -1, -2):
                                    if buf_k[b] == 0:
                                        buf_k[b] = classpreempt + 1      # 1-based class
                                        buf_k[b + 1] = stored_phase
                                        break
                                else:
                                    # no empty pair slot: simulation grows buffer on demand (prepend); CTMC pre-sizes to capacity so this only fires on the dynamically-growing SSA path.
                                    if is_simulation:
                                        buf_k = np.concatenate(
                                            ([classpreempt + 1, stored_phase], buf_k))
                                state_row = _compose_state(buf_k, srv_k, var_k)
                                out_states.append(state_row)
                                out_rates.append(-1.0)
                                out_probs.append(pentry[kentry] * count_preempt / srv_total)
                                # the displaced job leaves the server and the
                                # arriving one takes it, on the same arc
                                start_tags[len(out_states) - 1] = job_class + 1
                                preempt_tags[len(out_states) - 1] = classpreempt + 1
                                any_preempt = True
                    if not any_preempt:
                        # no lower-priority job to preempt: arrival WAITS in buffer, not dropped (previously invalidated row, BUG-70, starving highest-priority class queue).
                        buf_k = space_buf_k[i].copy() if space_buf_k.ndim >= 2 else space_buf_k.copy()
                        srv_k = space_srv_k[i].copy()
                        var_k = (space_var_k[i].copy() if space_var_k.ndim >= 2 and space_var_k.size > 0
                                 else (space_var_k.copy() if space_var_k.size > 0 else np.array([])))
                        cap_ok = True
                        if capacity is not None and ni_val >= capacity[ist]:
                            cap_ok = False
                        if classcap is not None and nir[i, job_class] >= classcap[ist, job_class]:
                            cap_ok = False
                        placed = False
                        if cap_ok:
                            # Fill the rightmost empty (class, phase) pair, matching
                            # the preemption store above and the right-aligned space.
                            for b in range(len(buf_k) - 2, -1, -2):
                                if buf_k[b] == 0:
                                    buf_k[b] = job_class + 1
                                    buf_k[b + 1] = kentry + 1
                                    placed = True
                                    break
                            if not placed and is_simulation:
                                buf_k = np.concatenate(([job_class + 1, kentry + 1], buf_k))
                                placed = True
                        if placed:
                            state_row = _compose_state(buf_k, srv_k, var_k)
                            out_states.append(state_row)
                            out_rates.append(-1.0)
                            out_probs.append(pentry[kentry])
                        elif cap_ok and not is_simulation:
                            # buffer width exhausted but capacity permits the job: a beyond-cutoff CTMC state, BLOCK; mirrors the FCFS width-full branch above.
                            valid[i] = False
                        elif not _arrival_is_lost(sn, ist, job_class):
                            # no room and the arrival cannot be lost (closed/BAS): block the row; an open class falls through to the unchanged-state loss.
                            valid[i] = False
                        # open-class capacity DROP: leave state unchanged (self-loop), no output row appended.
            continue

        else:
            # Default: enter service
            for i in range(n_rows):
                space_srv_k[i, int(Ks[job_class]) + kentry] += 1

        # Compose output states
        for i in range(n_rows):
            if valid[i]:
                state_row = _compose_state(
                    space_buf_k[i] if space_buf_k.ndim >= 2 else space_buf_k,
                    space_srv_k[i] if space_srv_k.ndim >= 2 else space_srv_k,
                    space_var_k[i] if space_var_k.ndim >= 2 and space_var_k.size > 0 else space_var_k)
                out_states.append(state_row)
                out_rates.append(-1.0)  # Passive action
                out_probs.append(outprob_k[i])
                if start_row[i] > 0:
                    start_tags[len(out_states) - 1] = int(start_row[i])

    # balking (QUEUE_LENGTH strategy) decided on pre-arrival population; admitted branches scaled by (1-balk_prob); mirrors MATLAB.
    if (getattr(sn, 'balkingStrategy', None) is not None and len(out_states) > 0
            and int(sn.balkingStrategy[ist, job_class]) == 1):
        balk_prob = 0.0
        thr = None
        if getattr(sn, 'balkingThresholds', None) is not None:
            thr = sn.balkingThresholds[ist][job_class]
        qlen = int(round(float(np.ravel(ni)[0])))
        if thr:
            for tup in thr:
                lo, hi, p = tup[0], tup[1], tup[2]
                if qlen >= lo and qlen <= hi:
                    balk_prob = float(p)
                    break
        if balk_prob > 0.0:
            out_probs = [p * (1.0 - balk_prob) for p in out_probs]
            bw = max(len(s) for s in out_states)
            inrow = np.ravel(inspace).astype(float)
            if bw > len(inrow):
                inrow = np.concatenate([np.zeros(bw - len(inrow)), inrow])
            out_states.append(inrow)
            out_rates.append(-1.0)
            out_probs.append(balk_prob)
            # a balked job never joins: the arc carries no tag

    return _finalize(out_states, out_rates, out_probs, inspace,
                     start_tags, preempt_tags, R)


# ---------------------------------------------------------------------------
# RETRY handler (retrial orbit)
# ---------------------------------------------------------------------------

def _handle_arv_reply(sn, ind, ist, job_class, K, Ks, S, pie,
                      space_buf, space_srv, space_var, inspace,
                      is_simulation=False):
    """Arrival of a REPLY signal class at the station holding its server.

    The reply completes the synchronous call:

      1. one held server of the calling class is released (the reply block
         counter is decremented);
      2. the released server is taken by the REPLY ITSELF, never by a waiting
         job.

    Unlike a NEGATIVE or CATASTROPHE signal the reply is NOT annihilated: it is
    a job that carries the call result onward, so it is served here (typically
    Immediate) and routed on by the ordinary routing matrix. This mirrors LDES,
    where a REPLY arrival frees the blocked server and then continues routing.

    The pass-through in step 2 is deliberate: the reply is the released work of
    a call this station already paid for, so it does not queue behind the
    residents (queueing it gave QLen 0.125 and residence 0.0996 against 0 in
    LDES, and stole capacity, costing 5% of throughput). Its service is
    typically Immediate, so the server is handed back at once and the ordinary
    FCFS departure path then promotes the head of line -- which also keeps
    sum(srv) <= S, unlike admitting the reply on top of a promoted job.

    The event is passive: the rate is set by the active departure at the callee.
    Port of MATLAB State/afterEventStationReply.m.
    """
    rinfo = _reply_block_info(sn, ind)

    # The calling class this reply releases: the class whose expected reply is
    # job_class and which holds a block here. sn.syncreply is 0-based.
    callclass = -1
    syncreply = np.asarray(sn.syncreply).ravel()
    for r in rinfo.classes:
        if r < syncreply.size and int(syncreply[r]) == job_class:
            callclass = r
            break

    srv_rows = np.atleast_2d(space_srv)
    buf_rows = np.atleast_2d(space_buf) if np.asarray(space_buf).size > 0 else None
    var_rows = np.atleast_2d(space_var)

    out_states, out_rates, out_probs = [], [], []
    start_tags = {}
    pentry = _get_pie(pie, ist, job_class, K)
    for i in range(srv_rows.shape[0]):
        buf = (buf_rows[i].copy() if buf_rows is not None and buf_rows.shape[0] > i
               else (buf_rows[0].copy() if buf_rows is not None else np.array([])))
        srv = srv_rows[i].copy()
        var = var_rows[i].copy() if var_rows.shape[0] > i else var_rows[0].copy()

        if callclass >= 0:
            col = int(rinfo.slot[callclass])
            if 0 <= col < var.size and var[col] > 0:
                var[col] -= 1

        # Join: into a free server (enumerating its entry phase) or, if all
        # remaining servers are busy, at the tail of the right-aligned buffer.
        nb = float(_reply_blocked(sn, ind, var.reshape(1, -1))[1][0])
        Seff = S - nb
        if float(np.sum(srv)) < Seff:
            for kentry in range(int(K[job_class])):
                if pentry[kentry] <= 0:
                    continue
                srv_k = srv.copy()
                srv_k[int(Ks[job_class]) + kentry] += 1
                out_states.append(_compose_state(buf, srv_k, var))
                out_rates.append(-1.0)
                out_probs.append(pentry[kentry])
                # the reply took the server it just released: a service start
                _tag_start(start_tags, out_states, job_class)
            continue

        buf_new = buf.copy()
        empty = np.where(buf_new == 0)[0]
        if empty.size == 0:
            buf_new = np.concatenate([[0.0], buf_new])
            emptypos = 0
        else:
            emptypos = int(empty[-1])
        buf_new[emptypos] = job_class + 1
        out_states.append(_compose_state(buf_new, srv, var))
        out_rates.append(-1.0)
        out_probs.append(1.0)

    R = int(getattr(sn, 'nclasses', 0) or 0)
    outspace, outrate, outprob, outstart, outpreempt = _finalize(
        out_states, out_rates, out_probs, inspace, start_tags, None, R)
    if is_simulation and outprob.shape[0] > 1:
        cum = np.cumsum(outprob.ravel()) / float(np.sum(outprob))
        idx = int(np.searchsorted(cum, np.random.rand()))
        idx = min(idx, outspace.shape[0] - 1)
        outspace = outspace[idx:idx + 1]
        outrate = outrate[idx:idx + 1]
        outprob = np.ones((1, 1))
        outstart = outstart[idx:idx + 1]
        outpreempt = outpreempt[idx:idx + 1]
    return outspace, outrate, outprob, outstart, outpreempt


def _handle_renege(sn, ind, inspace, job_class, ist, R, K, Ks,
                   space_buf, space_srv, space_var):
    """Exponential-patience reneging: each waiting (queued, not-in-service) class-r
    job abandons at memoryless rate impatienceMu, so the aggregate rate is
    (waiting count) * mu; one waiting job is removed and leaves the system
    (passive LOCAL). Mirrors MATLAB afterEventStation."""
    from .marginal import toMarginal
    _, nir, sir, _ = toMarginal(sn, ind, inspace)
    waiting = float(nir[0, job_class]) - float(sir[0, job_class])
    buf = np.ravel(space_buf).astype(float)
    if waiting > 0 and buf.size > 0:
        slot = next((c for c in range(len(buf)) if int(round(buf[c])) == job_class + 1), -1)
        if slot >= 0:
            buf_k = np.concatenate([[0.0], np.delete(buf, slot)])
            srv = np.ravel(space_srv).astype(float)
            var = np.ravel(space_var).astype(float) if space_var.size > 0 else np.array([])
            out = np.concatenate([buf_k, srv, var])
            rate = waiting * float(sn.impatienceMu[ist, job_class])
            # a WAITING job abandons: no server is freed and none is taken
            return (out.reshape(1, -1), np.array([[rate]]), np.array([[1.0]]),
                    np.zeros((1, R)), np.zeros((1, R)))
    return _empty_result(sn)


def _handle_switch(sn, ind, inspace, job_class, ist, R, K, Ks, pie,
                   space_buf, space_srv, space_var):
    """Switchover of a polling server walking towards buffer job_class.

    Unlike PHASE, which carries only the internal transitions of a phase-type
    and leaves the absorption to DEP, this event carries both: a completed
    switchover moves no job and so has no departure to attach the absorption to
    (see lang/sync). It is therefore also emitted for a single-phase switchover,
    where it consists of the absorption alone.
    """
    from .polling import polling_info, polling_get, polling_set, polling_next, polling_land

    pinfo = polling_info(sn, ind)
    if pinfo is None or not bool(pinfo['has_sw'][job_class]):
        return _empty_result(sn)

    inspace = np.atleast_2d(inspace)
    space_buf = np.atleast_2d(space_buf)
    space_srv = np.atleast_2d(space_srv)
    space_var = np.atleast_2d(space_var)

    out_states, out_rates, out_probs = [], [], []
    start_tags = {}
    D0 = pinfo['sw_d0'][job_class]
    D1 = pinfo['sw_d1'][job_class]

    for row in range(inspace.shape[0]):
        pos, swk, _ = polling_get(pinfo, space_var[row], -1)
        if pos != job_class or swk == 0:
            continue  # the server is not inside the switchover into job_class

        # internal transitions of the switchover phase-type
        for kdest in range(int(pinfo['ksw'][job_class])):
            if kdest == swk - 1 or D0[swk - 1, kdest] <= 0:
                continue
            var_k = polling_set(pinfo, space_var[row], job_class, kdest + 1, 0)
            out_states.append(np.concatenate([space_buf[row], space_srv[row], var_k]))
            out_rates.append(float(D0[swk - 1, kdest]))
            out_probs.append(1.0)

        # absorption: the server arrives at buffer job_class and either opens a
        # visit there or walks on
        rate = float(np.sum(D1[swk - 1, :]))
        if rate <= 0:
            continue
        nbuf = space_buf[row, :R]
        q, mode, budget = polling_next(pinfo, job_class, nbuf, R, arrived=True)
        rows, probs = polling_land(pinfo, q, mode, budget, space_buf[row],
                                   space_srv[row], space_var[row], K, Ks, pie, ist, R)
        for j in range(len(rows)):
            # switchover over empty buffer re-entering same phase of same switchover lands on departure state; self-loop suppressed, not inflating exit rate.
            if np.array_equal(rows[j], np.ravel(inspace[row])):
                continue
            out_states.append(rows[j])
            out_rates.append(rate * probs[j])
            out_probs.append(1.0)
            # A completed switchover that opens a visit pulls a waiting class-q
            # job into the server, so it starts service just as an ARV or a DEP
            # promotion does. This is the one service start a polling station
            # reaches through neither, and leaving it untagged would break
            # startRate == TN + preemptRate there for no reason other than the
            # name of the carrier event.
            if mode == 1:
                _tag_start(start_tags, out_states, q)

    if not out_states:
        return _empty_result(sn)
    return (np.array(out_states, dtype=float),
            np.array(out_rates, dtype=float).reshape(-1, 1),
            np.array(out_probs, dtype=float).reshape(-1, 1),
            _tag_matrix(start_tags, len(out_states), R),
            np.zeros((len(out_states), R)))


def _is_retrial_station(sn, ist):
    """True if station ist has a retrial orbit configured (retrialProc non-empty)."""
    rp = getattr(sn, 'retrialProc', None)
    if rp is None or ist >= len(rp):
        return False
    return any(x is not None for x in rp[ist])


def _handle_retry(sn, ind, inspace, job_class, ist, R, S, K, Ks, pie,
                  space_buf, space_srv, space_var, is_simulation=False):
    """Exponential retrial: an orbiting (buffered) class-r job retries entry. It
    succeeds only when a server is free, in which case one orbiting job enters
    service (entry phase from pie). Aggregate rate = (orbit) * retrialMu. Mirrors
    MATLAB afterEventStation."""
    from .marginal import toMarginal
    _, nir, sir, _ = toMarginal(sn, ind, inspace)
    orbit = float(nir[0, job_class]) - float(sir[0, job_class])
    in_srv = float(np.sum(space_srv))
    Sist = int(np.ravel(S)[ist]) if np.ndim(S) else int(S)
    out_states, out_rates, out_probs = [], [], []
    start_tags = {}
    buf = np.ravel(space_buf).astype(float)
    if orbit > 0 and in_srv < Sist and buf.size > 0:
        slot = -1
        for c in range(len(buf)):
            if int(round(buf[c])) == job_class + 1:
                slot = c
                break
        if slot >= 0:
            buf_k = np.concatenate([[0.0], np.delete(buf, slot)])
            pentry = _get_pie(pie, ist, job_class, K)
            retrial_mu = float(sn.retrialMu[ist, job_class])
            # LINEAR retrial policy scales the aggregate rate with orbit size (per-job timers); CONSTANT policy does not (one controller for the whole orbit).
            retrial_rate = orbit * retrial_mu
            rpol = getattr(sn, 'retrialPolicy', None)
            if rpol is not None:
                from ...lang.base import RetrialPolicy as _RetrialPolicy
                try:
                    if int(rpol[ist, job_class]) == int(_RetrialPolicy.CONSTANT):
                        retrial_rate = retrial_mu
                except (IndexError, TypeError, ValueError):
                    pass
            srv = np.ravel(space_srv).astype(float)
            var = np.ravel(space_var).astype(float) if space_var.size > 0 else np.array([])
            for kentry in range(int(K[job_class])):
                pe = float(pentry[kentry])
                if pe <= 0:
                    continue
                srv_k = srv.copy()
                srv_k[int(Ks[job_class]) + kentry] += 1
                out_states.append(np.concatenate([buf_k, srv_k, var]))
                out_rates.append(retrial_rate * pe)
                out_probs.append(1.0)
                # A successful retry is the only way into the server at a
                # retrial station (its departures never promote from the orbit),
                # so it carries the START the invariant needs.
                _tag_start(start_tags, out_states, job_class)
            if is_simulation and len(out_states) > 1:
                rates = np.array(out_rates)
                tot = rates.sum()
                cr = np.cumsum(rates) / tot
                rnd = np.random.rand()
                fc = 1 + max([-1] + [i for i in range(len(cr)) if rnd > cr[i]])
                sampled = _tag_matrix({0: start_tags.get(fc)}, 1, R) if fc in start_tags else np.zeros((1, R))
                return (np.array(out_states[fc]).reshape(1, -1),
                        np.array([[tot]]), np.array([[1.0]]),
                        sampled, np.zeros((1, R)))
    if len(out_states) == 0:
        return _empty_result(sn)
    return (np.array(out_states), np.array(out_rates).reshape(-1, 1),
            np.array(out_probs).reshape(-1, 1),
            _tag_matrix(start_tags, len(out_states), R), np.zeros((len(out_states), R)))


# ---------------------------------------------------------------------------
# DEP handler
# ---------------------------------------------------------------------------

def _handle_dep(sn, ind, inspace, job_class, M, R, ist, S, K, Ks,
                hasOnlyExp, mu, phi, pie, proc, ismkvmodclass,
                lldscaling, lldlimit, cdscaling, capacity, classcap,
                V, sched, space_buf, space_srv, space_var, no_promote=False):
    """Handle departure events."""
    from .marginal import toMarginal
    ni, nir, sir, kir = toMarginal(sn, ind, inspace, K, Ks,
                                    space_buf=space_buf, space_srv=space_srv,
                                    space_var=space_var)
    ni = np.atleast_1d(ni)
    nir = np.atleast_2d(nir)
    sir = np.atleast_2d(sir)

    # marked (MMAP) source class departs from the carrier's phase block (mark index 1) using its per-mark D1k matrix.
    markofclass = -1
    phclass = job_class
    if (sched == SchedStrategy.EXT and getattr(sn, 'markidx', None) is not None
            and ist < sn.markidx.shape[0] and sn.markidx[ist, job_class] > 0):
        markofclass = int(sn.markidx[ist, job_class])
        _carrier = np.where(sn.markidx[ist, :] == 1)[0]
        if _carrier.size > 0:
            phclass = int(_carrier[0])

    # Check any job of this class in service (marked: the carrier block)
    has_jobs = False
    for i in range(space_srv.shape[0]):
        for k_phase in range(int(K[phclass])):
            col = int(Ks[phclass]) + k_phase
            if col < space_srv.shape[1] and space_srv[i, col] > 0:
                has_jobs = True
                break
        if has_jobs:
            break
    if not has_jobs:
        return _empty_result(sn)

    # Update round-robin pointer
    space_var = space_var.copy()
    _maybe_update_rrobin(sn, ind, job_class, R, space_var)

    # Synchronous call: a departing job of this class keeps its server here
    # until its REPLY signal returns. -1 when this is an ordinary departure.
    _reply_slot = _reply_call_slot(sn, ind, job_class)

    out_states = []
    out_rates = []
    out_probs = []
    # START tags of the departure arcs: a completion hands the freed server to a
    # waiting job, and that promotion is the service start. A departure preempts
    # nobody, so there is no PREEMPT counterpart here. Keyed by successor index.
    start_tags = {}

    for k_phase in range(int(K[phclass])):
        col = int(Ks[phclass]) + k_phase
        # Find enabled rows
        en = space_srv[:, col] > 0
        if not np.any(en):
            continue

        enabled_rows = np.where(en)[0]

        for row_idx in enabled_rows:
            space_buf_k = space_buf[row_idx].copy() if space_buf.ndim >= 2 else space_buf.copy()
            space_srv_k = space_srv[row_idx].copy()
            space_var_k = space_var[row_idx].copy() if space_var.ndim >= 2 and space_var.size > 0 else (space_var.copy() if space_var.size > 0 else np.array([]))

            # Record departure
            space_srv_k[col] -= 1

            # Synchronous call: this departing job keeps its server until its
            # REPLY signal returns here, so the server is NOT handed to a
            # waiting job; it is recorded as held in the reply block instead.
            # Mirrors LDES, which omits the markServerIdle call and records a
            # pendingReply.
            no_promote_k = no_promote
            if _reply_slot >= 0:
                space_var_k = np.atleast_1d(space_var_k).astype(float).copy()
                if _reply_slot < space_var_k.size:
                    space_var_k[_reply_slot] += 1
                no_promote_k = True

            # Get kir for this row (marked: the token lives in the carrier block)
            if kir.ndim == 3:
                kir_val = kir[row_idx, phclass, k_phase]
            else:
                kir_val = sir[row_idx, phclass] if hasOnlyExp else 1.0

            ni_val = ni[row_idx]
            nir_row = nir[row_idx]

            # Compute rate based on scheduling strategy
            if markofclass > 0:
                # Marked source: this class's departure rate from shared phase
                # k is the row sum of its per-mark matrix D1k
                p_src_m = _get_proc(proc, ist, job_class)
                D1m = np.atleast_2d(np.asarray(p_src_m[1 + markofclass], dtype=float))
                lld_m = _get_lld(lldscaling, ist, ni_val, lldlimit)
                cd_m = _get_cd(cdscaling, ist, nir_row, job_class)
                rate = float(np.sum(D1m[k_phase, :])) * lld_m * cd_m
            else:
                rate = _compute_dep_rate(sched, ist, S, job_class, k_phase, K, Ks,
                                          mu, phi, proc, kir_val, ni_val, nir_row, sir[row_idx],
                                          lldscaling, lldlimit, cdscaling, R, sn)

            if rate == 0 or (rate < 0 and not _allows_signed_rates(sn, ist, job_class)):
                continue

            if sched == SchedStrategy.EXT:
                # Source DEP->kentry, prob D1[k_phase,kentry]/rowsum (MAP/MMPP2 autocorr, else pie); weight->rate (out_prob=1); CTMC sync skips active-side out_prob.
                pentry = None
                p_src = _get_proc(proc, ist, job_class)
                if markofclass > 0 and p_src is not None and len(p_src) > 1 + markofclass:
                    # Marked source: entry phase drawn from row k of the
                    # per-mark matrix D1k over the shared (carrier) chain
                    D1_src = np.atleast_2d(np.asarray(p_src[1 + markofclass], dtype=float))
                    row_sum = float(np.sum(D1_src[k_phase, :]))
                    if row_sum > 0:
                        pentry = D1_src[k_phase, :] / row_sum
                elif p_src is not None and len(p_src) > 1 and p_src[1] is not None:
                    D1_src = np.atleast_2d(np.asarray(p_src[1], dtype=float))
                    if (D1_src.shape[0] == D1_src.shape[1]
                            and D1_src.shape[0] == int(K[phclass])
                            and k_phase < D1_src.shape[0]):
                        row_sum = float(np.sum(D1_src[k_phase, :]))
                        if row_sum > 0:
                            pentry = D1_src[k_phase, :] / row_sum
                if pentry is None:
                    pentry = _get_pie(pie, ist, job_class, K)
                for kentry in range(int(K[phclass])):
                    p_k = pentry[kentry] if kentry < len(pentry) else 1.0
                    if p_k == 0:
                        continue
                    space_srv_e = space_srv_k.copy()
                    col_entry = int(Ks[phclass]) + kentry
                    space_srv_e[col_entry] += 1
                    state_row = _compose_state(space_buf_k, space_srv_e, space_var_k)
                    out_states.append(state_row)
                    out_rates.append(rate * p_k)
                    out_probs.append(1.0)
            elif sched in (SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS):
                if (ismkvmodclass is not None and job_class < len(ismkvmodclass)
                        and ismkvmodclass[job_class]):
                    # MAP service: correlated restart phase drawn from row k of D1
                    _dep_map_fcfs(sn, ind, ist, job_class, k_phase, K, Ks, S, proc,
                                  ismkvmodclass, space_buf_k, space_srv_k, space_var_k,
                                  kir_val, ni_val, nir_row, lldscaling, lldlimit,
                                  cdscaling, pie, out_states, out_rates, out_probs,
                                  no_promote=no_promote_k, start_tags=start_tags)
                else:
                    # Buffer promotion: move head-of-line job to service
                    _dep_with_buffer_promotion(sched, sn, space_buf_k, space_srv_k, space_var_k,
                                                K, Ks, pie, ist, S, ni_val, sir[row_idx],
                                                out_states, out_rates, out_probs, rate, R,
                                                ismkvmodclass=ismkvmodclass, no_promote=no_promote_k,
                                                start_tags=start_tags)
            elif sched in _SCHED_PREEMPT:
                # Promote a buffered job back to service, at its saved phase (PR) or from pie (PI)
                _dep_preemptive(sched, sn, space_buf_k, space_srv_k, space_var_k,
                                K, Ks, pie, ist, S, ni_val,
                                out_states, out_rates, out_probs, rate, R,
                                start_tags=start_tags)
            elif sched == SchedStrategy.POLLING:
                _dep_polling(sn, ind, space_buf_k, space_srv_k, space_var_k,
                              K, Ks, pie, ist, job_class,
                              out_states, out_rates, out_probs, rate, R,
                              start_tags=start_tags)
            elif sched == SchedStrategy.SIRO:
                _dep_siro_promotion(space_buf_k, space_srv_k, space_var_k,
                                     K, Ks, pie, ist, ni_val, nir_row, sir[row_idx],
                                     out_states, out_rates, out_probs, rate, R,
                                     no_promote=no_promote, start_tags=start_tags)
            elif sched in (SchedStrategy.SEPT, SchedStrategy.LEPT):
                _dep_sept_lept_promotion(sched, sn, space_buf_k, space_srv_k, space_var_k,
                                          K, Ks, pie, ist, ni_val,
                                          out_states, out_rates, out_probs, rate, R,
                                          no_promote=no_promote, start_tags=start_tags)
            else:
                # No buffer promotion (PS, INF, DPS, GPS, etc.)
                state_row = _compose_state(space_buf_k, space_srv_k, space_var_k)
                out_states.append(state_row)
                out_rates.append(rate)
                out_probs.append(1.0)

    return _finalize(out_states, out_rates, out_probs, inspace,
                     start_tags, None, R)


def _compute_dep_rate(sched, ist, S, job_class, k_phase, K, Ks,
                       mu, phi, proc, kir_val, ni_val, nir_row, sir_row,
                       lldscaling, lldlimit, cdscaling, R, sn):
    """Compute departure rate for a given scheduling strategy."""
    # Get mu and phi values
    mu_val = _get_mu(mu, ist, job_class, k_phase)
    phi_val = _get_phi(phi, ist, job_class, k_phase)

    if mu_val == 0 or np.isnan(mu_val):
        return 0.0

    # Load-dependent scaling
    lld = _get_lld(lldscaling, ist, ni_val, lldlimit)
    cd = _get_cd(cdscaling, ist, nir_row, job_class)

    if sched == SchedStrategy.EXT:
        # Source: rate = mu * phi * lld * cd
        return mu_val * phi_val * lld * cd

    elif sched == SchedStrategy.INF:
        return mu_val * phi_val * kir_val * lld * cd

    elif sched in (SchedStrategy.PS, SchedStrategy.LPS):
        if ni_val > 0:
            return mu_val * phi_val * kir_val * min(ni_val, S) / ni_val * lld * cd
        return 0.0

    elif sched == SchedStrategy.PSPRIO:
        # PS with priority behaves as plain PS when all jobs fit in the servers; otherwise only the most urgent present cohort shares the servers.
        if ni_val <= S:
            if ni_val > 0:
                return mu_val * phi_val * kir_val * min(ni_val, S) / ni_val * lld * cd
            return 0.0
        eligible, nirprio = _priority_nir(sn, nir_row, job_class, R)
        if not eligible:
            return 0.0
        niprio = float(np.sum(nirprio))
        if niprio > 0:
            return mu_val * phi_val * kir_val * min(niprio, S) / niprio * lld * cd
        return 0.0

    elif sched == SchedStrategy.DPSPRIO:
        # Discriminatory PS with priority.
        w = _get_sched_weights(sn, ist, R)
        if ni_val <= S:
            w_sum = np.dot(w, nir_row)
            if w_sum > 0 and nir_row[job_class] > 0:
                return mu_val * phi_val * (kir_val / nir_row[job_class]) \
                    * w[job_class] * nir_row[job_class] / w_sum * lld * cd
            return 0.0
        eligible, nirprio = _priority_nir(sn, nir_row, job_class, R)
        if not eligible:
            return 0.0
        w_sum = np.dot(w, nirprio)
        if w_sum > 0 and nirprio[job_class] > 0:
            return mu_val * phi_val * (kir_val / nirprio[job_class]) \
                * w[job_class] * nirprio[job_class] / w_sum * lld * cd
        return 0.0

    elif sched == SchedStrategy.GPSPRIO:
        # Generalized PS with priority.
        w = _get_sched_weights(sn, ist, R)
        if ni_val <= S:
            cir = np.minimum(nir_row, 1)
            w_cir = np.dot(w, cir)
            if w_cir > 0 and nir_row[job_class] > 0:
                return mu_val * phi_val * (kir_val / nir_row[job_class]) \
                    * w[job_class] / w_cir * lld * cd
            return 0.0
        eligible, nirprio = _priority_nir(sn, nir_row, job_class, R)
        if not eligible:
            return 0.0
        cir = np.minimum(nirprio, 1)
        w_cir = np.dot(w, cir)
        if w_cir > 0 and nirprio[job_class] > 0:
            return mu_val * phi_val * (kir_val / nirprio[job_class]) \
                * w[job_class] / w_cir * lld * cd
        return 0.0

    elif sched == SchedStrategy.POLLING:
        # Single server dedicated to the buffer being visited: the completion
        # rate is that of the one job in service.
        return mu_val * phi_val * kir_val * lld * cd

    elif sched in ((SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS,
                    SchedStrategy.SIRO) + _SCHED_PREEMPT):
        # For FCFS with MAP processes, use D1 matrix
        if proc is not None:
            p = _get_proc(proc, ist, job_class)
            if p is not None and len(p) > 1:
                d1 = np.atleast_2d(p[1])
                # Sum D1 row k_phase to get total departure rate from this phase
                rate = np.sum(d1[k_phase, :]) * kir_val * lld * cd
                return rate

        return mu_val * phi_val * kir_val * lld * cd

    elif sched == SchedStrategy.DPS:
        # Discriminatory Processor Sharing
        w = _get_sched_weights(sn, ist, R)
        w_sum = np.dot(w, nir_row)
        if w_sum > 0 and nir_row[job_class] > 0:
            return mu_val * phi_val * (kir_val / nir_row[job_class]) * w[job_class] * nir_row[job_class] / w_sum * lld * cd
        return 0.0

    elif sched == SchedStrategy.GPS:
        # Generalized Processor Sharing
        w = _get_sched_weights(sn, ist, R)
        cir = np.minimum(nir_row, 1)
        w_cir = np.dot(w, cir)
        if w_cir > 0 and nir_row[job_class] > 0:
            return mu_val * phi_val * (kir_val / nir_row[job_class]) * w[job_class] / w_cir * lld * cd
        return 0.0

    # Default
    return mu_val * phi_val * kir_val * lld * cd


def _dep_with_buffer_promotion(sched, sn, space_buf_k, space_srv_k, space_var_k,
                                K, Ks, pie, ist, S, ni_val, sir_row,
                                out_states, out_rates, out_probs, rate, R,
                                ismkvmodclass=None, no_promote=False, start_tags=None):
    """Handle departure with buffer promotion for FCFS/HOL/LCFS.

    When no_promote is True (immediate-feedback self-loop), the vacated server
    is left idle and no waiting job is promoted, so the fed-back job re-entering
    via the passive arrival holds the server instead of re-queueing.
    """
    buf = space_buf_k
    has_buffer_job = False
    # Retrial orbit: a freed server is NOT filled from the orbit (which is held in
    # the buffer); orbiting jobs re-enter only through RETRY events. Mirrors MATLAB.
    if buf.size > 0 and not (no_promote or _is_retrial_station(sn, ist)):
        has_buffer_job = np.any(buf > 0)

    if not has_buffer_job:
        # No buffer jobs: just departure
        state_row = _compose_state(space_buf_k, space_srv_k, space_var_k)
        out_states.append(state_row)
        out_rates.append(rate)
        out_probs.append(1.0)
        return

    # Find the job to promote from buffer
    if sched == SchedStrategy.FCFS:
        # FCFS: rightmost non-zero = head of line
        promote_idx = -1
        for b in range(len(buf) - 1, -1, -1):
            if buf[b] > 0:
                promote_idx = b
                break
    elif sched == SchedStrategy.LCFS:
        # plain LCFS promotes newest arrival ignoring priority; LCFSPRIO is the priority-aware branch. See _kb/11-conventions-and-gotchas.md Modeling traps.
        promote_idx = -1
        for b in range(len(buf)):
            if buf[b] > 0:
                promote_idx = b
                break
    elif sched == SchedStrategy.HOL:
        # HOL: rightmost with highest priority (lowest classprio value)
        promote_idx = -1
        best_prio = float('inf')
        classprio = sn.classprio if hasattr(sn, 'classprio') and sn.classprio is not None else np.zeros(R)
        for b in range(len(buf) - 1, -1, -1):
            if buf[b] > 0:
                cls = int(buf[b]) - 1
                prio = classprio[cls] if cls < len(classprio) else float('inf')
                if prio < best_prio:
                    best_prio = prio
                    promote_idx = b
        if promote_idx == -1:
            # Fallback to FCFS
            for b in range(len(buf) - 1, -1, -1):
                if buf[b] > 0:
                    promote_idx = b
                    break
    else:
        promote_idx = -1

    if promote_idx < 0:
        state_row = _compose_state(space_buf_k, space_srv_k, space_var_k)
        out_states.append(state_row)
        out_rates.append(rate)
        out_probs.append(1.0)
        return

    start_svc_class = int(buf[promote_idx]) - 1  # Convert from 1-based to 0-based
    buf_new = buf.copy()
    # Remove from buffer: shift remaining left
    if sched == SchedStrategy.FCFS:
        # Shift left from the right end
        buf_new[promote_idx] = 0
        # Shift: move buf_new[0:promote_idx] right by prepending 0
        buf_new[1:promote_idx + 1] = buf[0:promote_idx]
        buf_new[0] = 0
    elif sched == SchedStrategy.LCFS:
        buf_new[promote_idx] = 0
    elif sched == SchedStrategy.HOL:
        buf_new[promote_idx] = 0
        buf_new[1:promote_idx + 1] = buf[0:promote_idx]
        buf_new[0] = 0

    # a promoted MAP job resumes its persistent server-phase entry; a PH job enters per pie. Mirrors MATLAB afterEventStation.m:468-474.
    if (ismkvmodclass is not None and start_svc_class < len(ismkvmodclass)
            and ismkvmodclass[start_svc_class]):
        _mv = int(np.sum(np.asarray(ismkvmodclass)[:start_svc_class]))
        _sv = np.atleast_1d(space_var_k)
        kentry = int(_sv[_mv]) if _sv.size > _mv else 0
        srv_new = space_srv_k.copy()
        srv_new[int(Ks[start_svc_class]) + kentry] += 1
        out_states.append(_compose_state(buf_new, srv_new, space_var_k))
        out_rates.append(rate)
        out_probs.append(1.0)
        _tag_start(start_tags, out_states, start_svc_class)
    else:
        pentry = _get_pie(pie, ist, start_svc_class, K)
        for kentry in range(int(K[start_svc_class])):
            srv_new = space_srv_k.copy()
            srv_new[int(Ks[start_svc_class]) + kentry] += 1
            state_row = _compose_state(buf_new, srv_new, space_var_k)
            out_states.append(state_row)
            out_rates.append(rate * pentry[kentry])
            out_probs.append(1.0)
            _tag_start(start_tags, out_states, start_svc_class)


def _dep_map_fcfs(sn, ind, ist, job_class, k_phase, K, Ks, S, proc, ismkvmodclass,
                  space_buf_k, space_srv_k, space_var_k, kir_val, ni_val, nir_row,
                  lldscaling, lldlimit, cdscaling, pie, out_states, out_rates, out_probs,
                  no_promote=False, start_tags=None):
    """MAP service completion at an FCFS station.

    Port of MATLAB State.afterEventStation case DEP / SchedStrategy.FCFS
    (afterEventStation.m:442-528). The departing job completes from modulating
    phase k_phase; the restart phase kdest is taken from row k_phase of D1
    (correlated, not from pie). The persistent server-phase variable is set to
    kdest, the completion rate to kdest is D1(k,kdest)*kir*scaling, and when a
    job is promoted from the buffer a same-class MAP job resumes in phase kdest,
    a different-class MAP job resumes from its own stored phase, and an i.i.d.
    (PH) job enters per its pie.

    `space_srv_k` already has the departing job removed at col Ks(job_class)+k.
    """
    p = _get_proc(proc, ist, job_class)
    if p is None or len(p) < 2:
        return
    D1 = np.atleast_2d(np.asarray(p[1], dtype=float))
    n_phases = D1.shape[0]
    lld = _get_lld(lldscaling, ist, ni_val, lldlimit)
    cd = _get_cd(cdscaling, ist, nir_row, job_class)
    mvcol = int(np.sum(np.asarray(ismkvmodclass)[:job_class]))

    # FCFS buffered promotion is rightmost non-zero; immediate feedback (no_promote) holds the server so no waiting job is promoted.
    buf = space_buf_k
    promote_idx = -1
    if buf.size > 0 and not no_promote:
        for b in range(len(buf) - 1, -1, -1):
            if buf[b] > 0:
                promote_idx = b
                break

    for kdest in range(n_phases):
        rate_kd = float(D1[k_phase, kdest]) * kir_val * lld * cd
        if rate_kd == 0 or (rate_kd < 0 and not _allows_signed_rates(sn, ist, job_class)):
            continue
        var_kd = space_var_k.copy() if space_var_k.size > 0 else space_var_k
        if var_kd.size > mvcol:
            var_kd[mvcol] = kdest

        if promote_idx < 0:
            # Server goes idle; the modulating phase is remembered in var_kd.
            out_states.append(_compose_state(buf, space_srv_k, var_kd))
            out_rates.append(rate_kd)
            out_probs.append(1.0)
            continue

        # Promote the head-of-line buffered job.
        scls = int(buf[promote_idx]) - 1  # 0-based class
        buf_new = buf.copy()
        buf_new[promote_idx] = 0
        buf_new[1:promote_idx + 1] = buf[0:promote_idx]
        buf_new[0] = 0

        if ismkvmodclass is not None and scls < len(ismkvmodclass) and ismkvmodclass[scls]:
            if scls == job_class:
                kentry_range = [kdest]            # resume in the phase just vacated
            else:
                scls_mv = int(np.sum(np.asarray(ismkvmodclass)[:scls]))
                kentry_range = [int(var_kd[scls_mv])] if var_kd.size > scls_mv else [0]
            pentry = None
        else:
            kentry_range = range(int(K[scls]))
            pentry = _get_pie(pie, ist, scls, K)

        for kentry in kentry_range:
            srv_new = space_srv_k.copy()
            srv_new[int(Ks[scls]) + kentry] += 1
            pr = 1.0 if pentry is None else float(pentry[kentry])
            if pr <= 0:
                continue
            out_states.append(_compose_state(buf_new, srv_new, var_kd))
            out_rates.append(rate_kd * pr)
            out_probs.append(1.0)
            _tag_start(start_tags, out_states, scls)


def _dep_sept_lept_promotion(sched, sn, space_buf_k, space_srv_k, space_var_k,
                              K, Ks, pie, ist, ni_val,
                              out_states, out_rates, out_probs, rate, R,
                              no_promote=False, start_tags=None):
    """SEPT/LEPT departure promotion (non-preemptive, size-based).

    The next job to enter service is the buffered class with the SHORTEST (SEPT)
    or LONGEST (LEPT) expected processing time, i.e. the highest (SEPT) or lowest
    (LEPT) service rate sn.rates[ist, r]. The buffer holds per-class counts (the
    SIRO format). Ties break on the lowest class index (deterministic). Mirrors
    MATLAB afterEventStation.m, where sn.schedparam encodes the class service-time
    order and the first buffered class in that order is promoted.
    """
    buf = space_buf_k
    total_buffered = 0 if no_promote else (np.sum(buf) if buf.size > 0 else 0)
    if total_buffered <= 0:
        out_states.append(_compose_state(space_buf_k, space_srv_k, space_var_k))
        out_rates.append(rate)
        out_probs.append(1.0)
        return

    rates_row = np.asarray(sn.rates[ist], dtype=float).reshape(-1)
    cands = [r for r in range(R) if (r < len(buf) and buf[r] > 0)]
    if not cands:
        out_states.append(_compose_state(space_buf_k, space_srv_k, space_var_k))
        out_rates.append(rate)
        out_probs.append(1.0)
        return

    def _rate_of(r):
        v = rates_row[r] if r < len(rates_row) else np.nan
        return v if np.isfinite(v) else (-np.inf if sched == SchedStrategy.SEPT else np.inf)

    if sched == SchedStrategy.SEPT:
        # highest rate (shortest mean); tie -> lowest class index
        pick = max(cands, key=lambda r: (_rate_of(r), -r))
    else:  # LEPT: lowest rate (longest mean); tie -> lowest class index
        pick = min(cands, key=lambda r: (_rate_of(r), r))

    buf_new = buf.copy()
    buf_new[pick] -= 1
    pentry = _get_pie(pie, ist, pick, K)
    for kentry in range(int(K[pick])):
        srv_new = space_srv_k.copy()
        srv_new[int(Ks[pick]) + kentry] += 1
        out_states.append(_compose_state(buf_new, srv_new, space_var_k))
        out_rates.append(rate * pentry[kentry])
        out_probs.append(1.0)
        _tag_start(start_tags, out_states, pick)


def _dep_polling(sn, ind, space_buf_k, space_srv_k, space_var_k,
                  K, Ks, pie, ist, job_class,
                  out_states, out_rates, out_probs, rate, R, start_tags=None):
    """Handle a departure at a polling station.

    A completion ends the visit unless the discipline still allows another job of
    the same class to be taken; when it ends, the server walks the cyclic order
    to wherever the next tangible controller state lies. space_var_k already
    carries any round-robin pointer advance applied by the caller, so it is used
    as-is rather than re-read from the input state.
    """
    from ...constants import PollingType
    from .polling import polling_info, polling_get, polling_next, polling_land

    pinfo = polling_info(sn, ind)
    if pinfo is None:
        return

    _, swk, ctr = polling_get(pinfo, space_var_k, job_class)
    if swk != 0:
        return  # no job can complete while the server is walking

    nbuf = np.ravel(space_buf_k)[:R]
    ptype = pinfo['ptype']
    if ptype == PollingType.EXHAUSTIVE:
        ctrnext = 0
        goon = nbuf[job_class] > 0
    elif ptype == PollingType.GATED:
        ctrnext = ctr - 1  # one of the gated jobs completed
        goon = ctrnext > 0
    elif ptype == PollingType.KLIMITED:
        ctrnext = ctr - 1  # one of the K permitted services used
        goon = ctrnext > 0 and nbuf[job_class] > 0
    else:  # DECREMENTING
        ctrnext = ctr  # the target level is fixed for the visit
        goon = nbuf[job_class] > ctr

    if goon:
        q, mode, budget = job_class, 1, ctrnext
    else:
        q, mode, budget = polling_next(pinfo, job_class, nbuf, R, arrived=False)

    rows, probs = polling_land(pinfo, q, mode, budget, space_buf_k, space_srv_k,
                               space_var_k, K, Ks, pie, ist, R)
    for j in range(len(rows)):
        out_states.append(rows[j])
        out_rates.append(rate * probs[j])
        out_probs.append(1.0)
        # mode 1 pulls a waiting class-q job into the server; a switchover or a
        # park starts nobody
        if mode == 1:
            _tag_start(start_tags, out_states, q)


def _dep_siro_promotion(space_buf_k, space_srv_k, space_var_k,
                          K, Ks, pie, ist, ni_val, nir_row, sir_row,
                          out_states, out_rates, out_probs, rate, R, no_promote=False,
                          start_tags=None):
    """Handle departure with SIRO buffer promotion.

    When no_promote is True (immediate-feedback self-loop) the vacated server is
    left idle and no waiting job is promoted; the fed-back job re-enters via the
    passive arrival and holds the server.
    """
    buf = space_buf_k
    total_buffered = 0 if no_promote else (np.sum(buf) if buf.size > 0 else 0)

    if total_buffered <= 0:
        state_row = _compose_state(space_buf_k, space_srv_k, space_var_k)
        out_states.append(state_row)
        out_rates.append(rate)
        out_probs.append(1.0)
        return

    # Pick class proportional to buffer occupancy
    for r in range(R):
        buf_r = buf[r] if r < len(buf) else 0
        if buf_r <= 0:
            continue
        pick_prob = buf_r / total_buffered
        buf_new = buf.copy()
        buf_new[r] -= 1

        pentry = _get_pie(pie, ist, r, K)
        for kentry in range(int(K[r])):
            srv_new = space_srv_k.copy()
            srv_new[int(Ks[r]) + kentry] += 1
            state_row = _compose_state(buf_new, srv_new, space_var_k)
            out_states.append(state_row)
            out_rates.append(rate * pick_prob * pentry[kentry])
            out_probs.append(1.0)
            _tag_start(start_tags, out_states, r)


def _dep_preemptive(sched, sn, space_buf_k, space_srv_k, space_var_k,
                    K, Ks, pie, ist, S, ni_val,
                    out_states, out_rates, out_probs, rate, R, start_tags=None):
    """Handle departure with buffer promotion for the eight preempt disciplines.

    Buffer stores [class, phase, class, phase, ...] pairs (1-based values).
    On departure a buffered job is promoted back to service. The eight
    disciplines factor as {FCFS,LCFS} x {PR,PI} x {plain,PRIO}, so three flags
    cover what MATLAB afterEventStation.m writes as eight separate cases
    (:1043-1400): LCFS scans the buffer newest-first, PRIO restricts the scan to
    the highest-priority class, and PR resumes at the saved phase whereas PI
    discards it and restarts from the entry distribution pie.
    """
    buf = space_buf_k
    # Check if any job in buffer (class values at even indices)
    has_buffer_job = False
    if buf.size > 0:
        for b in range(0, len(buf), 2):
            if buf[b] > 0:
                has_buffer_job = True
                break

    if not has_buffer_job:
        # No buffer jobs: just departure
        state_row = _compose_state(space_buf_k, space_srv_k, space_var_k)
        out_states.append(state_row)
        out_rates.append(rate)
        out_probs.append(1.0)
        return

    # Find job to promote from buffer based on scheduling strategy
    classprio = sn.classprio if hasattr(sn, 'classprio') and sn.classprio is not None else np.zeros(R)

    if sched not in _SCHED_PREEMPT:
        target_col = -1
    else:
        lcfs = sched in _SCHED_PREEMPT_LCFS
        prio_aware = sched in _SCHED_PREEMPT_PRIO
        # LCFS scans newest-first (leftmost), FCFS oldest-first (rightmost).
        cols = range(0, len(buf), 2) if lcfs else range(len(buf) - 2, -1, -2)
        best_prio = float('inf')
        target_col = -1
        for b in cols:
            if buf[b] <= 0:
                continue
            if not prio_aware:
                # priority-awareness is a property of the declared policy, never
                # inferred from the data: the plain variants promote by arrival
                # order alone, the PRIO variants are the priority-aware ones.
                # See _kb/11-conventions-and-gotchas.md Modeling traps.
                target_col = b
                break
            cls = int(buf[b]) - 1  # 0-based class
            prio = classprio[cls] if cls < len(classprio) else float('inf')
            if prio < best_prio:
                # strict <, so the first column reached in scan order wins ties
                best_prio = prio
                target_col = b

    if target_col < 0:
        state_row = _compose_state(space_buf_k, space_srv_k, space_var_k)
        out_states.append(state_row)
        out_rates.append(rate)
        out_probs.append(1.0)
        return

    start_svc_class = int(buf[target_col]) - 1       # 0-based class

    # Remove [class, phase] pair from buffer, keeping it right-aligned
    buf_new = np.concatenate([
        np.array([0.0, 0.0]),
        buf[:target_col],
        buf[target_col + 2:]
    ])

    if sched in _SCHED_PREEMPT_RESUME:
        # PR: the promoted job resumes in the phase it was preempted at
        kentry_phase = int(buf[target_col + 1]) - 1   # 0-based phase
        srv_new = space_srv_k.copy()
        srv_new[int(Ks[start_svc_class]) + kentry_phase] += 1
        state_row = _compose_state(buf_new, srv_new, space_var_k)
        out_states.append(state_row)
        out_rates.append(rate)
        out_probs.append(1.0)
        _tag_start(start_tags, out_states, start_svc_class)
        return

    # PI: the saved phase is discarded and the job restarts from pie, so the
    # departure branches into one successor per entry phase, the rate split by
    # the entry probability (MATLAB afterEventStation.m:1100-1114).
    pentry = _get_pie(pie, ist, start_svc_class, K)
    for kentry in range(int(K[start_svc_class])):
        srv_new = space_srv_k.copy()
        srv_new[int(Ks[start_svc_class]) + kentry] += 1
        state_row = _compose_state(buf_new, srv_new, space_var_k)
        out_states.append(state_row)
        out_rates.append(rate * pentry[kentry])
        out_probs.append(1.0)
        _tag_start(start_tags, out_states, start_svc_class)


# ---------------------------------------------------------------------------
# PHASE handler
# ---------------------------------------------------------------------------

def _handle_phase(sn, ind, inspace, job_class, M, R, ist, S, K, Ks,
                  hasOnlyExp, mu, phi, pie, proc, ismkvmodclass,
                  lldscaling, lldlimit, cdscaling,
                  V, sched, space_buf, space_srv, space_var):
    """Handle phase transition events (internal to service)."""
    from .marginal import toMarginal
    ni, nir, sir, kir = toMarginal(sn, ind, inspace, K, Ks,
                                    space_buf=space_buf, space_srv=space_srv,
                                    space_var=space_var)
    ni = np.atleast_1d(ni)
    nir = np.atleast_2d(nir)

    if nir[0, job_class] <= 0:
        return _empty_result(sn)

    # Get D0 matrix for phase transitions
    d0 = None
    p = _get_proc(proc, ist, job_class)
    if p is not None and len(p) > 0:
        d0 = np.atleast_2d(p[0])

    lld_val = _get_lld(lldscaling, ist, ni[0], lldlimit)
    cd_val = _get_cd(cdscaling, ist, nir[0], job_class)

    out_states = []
    out_rates = []
    out_probs = []

    for k in range(int(K[job_class])):
        col_k = int(Ks[job_class]) + k
        en = space_srv[:, col_k] > 0
        if not np.any(en):
            continue

        enabled_rows = np.where(en)[0]

        for kdest in range(int(K[job_class])):
            if kdest == k:
                continue

            col_kdest = int(Ks[job_class]) + kdest

            for row_idx in enabled_rows:
                srv_k = space_srv[row_idx].copy()
                buf_k = space_buf[row_idx].copy() if space_buf.ndim >= 2 else space_buf.copy()
                var_k = space_var[row_idx].copy() if space_var.ndim >= 2 and space_var.size > 0 else (space_var.copy() if space_var.size > 0 else np.array([]))

                if kir.ndim == 3:
                    kir_val = kir[row_idx, job_class, k]
                else:
                    kir_val = 1.0

                # Move from phase k to kdest
                srv_k[col_k] -= 1
                srv_k[col_kdest] += 1

                # MAP: keep the persistent server-phase variable consistent with
                # the modulating-phase move (MATLAB afterEventStation.m:1074-1076).
                if (ismkvmodclass is not None and job_class < len(ismkvmodclass)
                        and ismkvmodclass[job_class] and var_k.size > 0):
                    _mvcol = int(np.sum(np.asarray(ismkvmodclass)[:job_class]))
                    if _mvcol < var_k.size:
                        var_k[_mvcol] = kdest

                # Rate from D0 matrix
                if d0 is not None:
                    d0_rate = d0[k, kdest]
                else:
                    d0_rate = 0.0

                # Apply scheduling-specific scaling
                rate = _compute_phase_rate(sched, d0_rate, kir_val, ni[row_idx],
                                            nir[row_idx], S, job_class, R, sn, ist,
                                            lld_val, cd_val)

                if rate == 0:
                    continue

                state_row = _compose_state(buf_k, srv_k, var_k)
                out_states.append(state_row)
                out_rates.append(rate)
                out_probs.append(1.0)

    return _finalize(out_states, out_rates, out_probs, inspace, None, None, R)


def _compute_phase_rate(sched, d0_rate, kir_val, ni_val, nir_row, S,
                         job_class, R, sn, ist, lld, cd):
    """Compute phase transition rate based on scheduling strategy."""
    if sched == SchedStrategy.EXT:
        return d0_rate * lld * cd

    elif sched == SchedStrategy.INF:
        return d0_rate * kir_val * lld * cd

    elif sched in (SchedStrategy.PS, SchedStrategy.LPS):
        if ni_val > 0:
            return d0_rate * kir_val * min(ni_val, S) / ni_val * lld * cd
        return 0.0

    elif sched == SchedStrategy.PSPRIO:
        if ni_val <= S:
            if ni_val > 0:
                return d0_rate * kir_val * min(ni_val, S) / ni_val * lld * cd
            return 0.0
        eligible, nirprio = _priority_nir(sn, nir_row, job_class, R)
        if not eligible:
            return 0.0
        niprio = float(np.sum(nirprio))
        if niprio > 0:
            return d0_rate * kir_val * min(niprio, S) / niprio * lld * cd
        return 0.0

    elif sched == SchedStrategy.DPS:
        w = _get_sched_weights(sn, ist, R)
        w_sum = np.dot(w, nir_row)
        if w_sum > 0:
            return d0_rate * kir_val * w[job_class] / w_sum * lld * cd
        return 0.0

    elif sched == SchedStrategy.DPSPRIO:
        w = _get_sched_weights(sn, ist, R)
        if ni_val <= S:
            w_sum = np.dot(w, nir_row)
            if w_sum > 0:
                return d0_rate * kir_val * w[job_class] / w_sum * lld * cd
            return 0.0
        eligible, nirprio = _priority_nir(sn, nir_row, job_class, R)
        if not eligible:
            return 0.0
        w_sum = np.dot(w, nirprio)
        if w_sum > 0:
            return d0_rate * kir_val * w[job_class] / w_sum * lld * cd
        return 0.0

    elif sched == SchedStrategy.GPS:
        w = _get_sched_weights(sn, ist, R)
        cir = np.minimum(nir_row, 1)
        w_cir = np.dot(w, cir)
        if w_cir > 0 and nir_row[job_class] > 0:
            return d0_rate * kir_val / nir_row[job_class] * w[job_class] / w_cir * lld * cd
        return 0.0

    elif sched == SchedStrategy.GPSPRIO:
        w = _get_sched_weights(sn, ist, R)
        if ni_val <= S:
            cir = np.minimum(nir_row, 1)
            w_cir = np.dot(w, cir)
            if w_cir > 0 and nir_row[job_class] > 0:
                return d0_rate * kir_val / nir_row[job_class] * w[job_class] / w_cir * lld * cd
            return 0.0
        eligible, nirprio = _priority_nir(sn, nir_row, job_class, R)
        if not eligible:
            return 0.0
        cir = np.minimum(nirprio, 1)
        w_cir = np.dot(w, cir)
        if w_cir > 0 and nirprio[job_class] > 0:
            return d0_rate * kir_val / nirprio[job_class] * w[job_class] / w_cir * lld * cd
        return 0.0

    else:
        # FCFS, HOL, LCFS, SIRO, SEPT, LEPT: D0 * kir
        return d0_rate * kir_val * lld * cd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _allows_signed_rates(sn, ist, job_class) -> bool:
    """True when the (station, class) process is a matrix exponential.

    An ME embeds in the generator exactly as a phase-type does, except that the
    off-diagonal entries of D0 and the completion vector -A*e may be negative:
    the stationary vector is then a signed measure whose aggregates over each
    phase block are still the exact probabilities. Dropping those negative
    entries breaks the balance equations, so the usual "rate <= 0 is not a
    transition" filter must be relaxed to "rate == 0" for such a station-class.
    See sn.isph, sn_is_phasetype and _kb/04-networkstruct.md.
    """
    isph = getattr(sn, 'isph', None)
    if isph is None:
        return False
    try:
        return not bool(np.asarray(isph)[int(ist), int(job_class)])
    except (IndexError, TypeError, ValueError):
        return False


def _compose_state(buf, srv, var):
    """Compose full state row from buf, srv, var parts."""
    parts = []
    if isinstance(buf, np.ndarray) and buf.size > 0:
        parts.append(buf.ravel())
    if isinstance(srv, np.ndarray) and srv.size > 0:
        parts.append(srv.ravel())
    if isinstance(var, np.ndarray) and var.size > 0:
        parts.append(var.ravel())
    return np.concatenate(parts) if parts else np.array([])


def _finalize(out_states, out_rates, out_probs, inspace,
              start_tags=None, preempt_tags=None, R=None):
    """Finalize output arrays.

    START_TAGS and PREEMPT_TAGS map the INDEX of a successor in OUT_STATES to
    the 1-based class it tags, so only the arcs that actually start or preempt
    a job have to say anything; every other successor is a zero row. The two
    tag matrices come back with one row per successor and one column per class,
    so a caller can index them exactly like OUTSPACE.
    """
    ncls = int(R) if R is not None else 0
    if len(out_states) == 0:
        return (np.zeros((0, 0)), np.zeros((0, 0)), np.zeros((0, 0)),
                np.zeros((0, ncls)), np.zeros((0, ncls)))
    outspace = np.array(out_states)
    outrate = np.array(out_rates).reshape(-1, 1)
    outprob = np.array(out_probs).reshape(-1, 1)
    outstart = _tag_matrix(start_tags, outspace.shape[0], ncls)
    outpreempt = _tag_matrix(preempt_tags, outspace.shape[0], ncls)
    return outspace, outrate, outprob, outstart, outpreempt


def _tag_start(start_tags, out_states, cls0):
    """Tag the successor just appended to OUT_STATES as starting a class-CLS0
    (0-based) service. A no-op when the caller did not ask for tags."""
    if start_tags is None or cls0 is None or cls0 < 0:
        return
    start_tags[len(out_states) - 1] = int(cls0) + 1


def _empty_result(sn):
    """The no-successor result, with tag matrices of the right width so that a
    caller can concatenate or index them without special-casing emptiness."""
    ncls = int(getattr(sn, 'nclasses', 0) or 0)
    return (np.zeros((0, 0)), np.zeros((0, 0)), np.zeros((0, 0)),
            np.zeros((0, ncls)), np.zeros((0, ncls)))


def _tag_matrix(tags, nrows, ncls):
    """Materialize an index -> 1-based-class tag map as an (nrows x ncls) count
    matrix. An absent index is an untagged arc, i.e. a zero row."""
    out = np.zeros((nrows, ncls))
    if not tags or ncls == 0:
        return out
    for idx, cls in tags.items():
        if idx is None or cls is None:
            continue
        i = int(idx)
        c = int(cls)
        if 0 <= i < nrows and 1 <= c <= ncls:
            out[i, c - 1] += 1
    return out


def _get_pie(pie, ist, job_class, K):
    """Get phase entry probabilities."""
    if pie is not None and ist in pie and job_class in pie[ist]:
        pe = pie[ist][job_class]
        if isinstance(pe, np.ndarray):
            return pe.ravel()[:int(K[job_class])]
        elif isinstance(pe, (list, tuple)):
            return np.array(pe[:int(K[job_class])])
    # Default: all probability in phase 0
    pentry = np.zeros(int(K[job_class]))
    if len(pentry) > 0:
        pentry[0] = 1.0
    return pentry


def _get_mu(mu, ist, job_class, k_phase):
    """Get service rate for a given station, class, phase."""
    if mu is not None and ist in mu and job_class in mu[ist]:
        m = mu[ist][job_class]
        if isinstance(m, np.ndarray):
            m = m.ravel()
            if k_phase < len(m):
                return float(m[k_phase])
        elif isinstance(m, (list, tuple)):
            if k_phase < len(m):
                return float(m[k_phase])
        else:
            return float(m)
    return 0.0


def _get_phi(phi, ist, job_class, k_phase):
    """Get completion probability for a given station, class, phase."""
    if phi is not None and ist in phi and job_class in phi[ist]:
        p = phi[ist][job_class]
        if isinstance(p, np.ndarray):
            p = p.ravel()
            if k_phase < len(p):
                return float(p[k_phase])
        elif isinstance(p, (list, tuple)):
            if k_phase < len(p):
                return float(p[k_phase])
        else:
            return float(p)
    return 1.0


def _get_proc(proc, ist, job_class):
    """Get process matrices [D0, D1] for a given station and class."""
    if proc is None:
        return None
    try:
        if isinstance(proc, dict):
            if ist in proc and job_class in proc[ist]:
                return proc[ist][job_class]
        elif isinstance(proc, (list, tuple)):
            if ist < len(proc) and proc[ist] is not None:
                station_proc = proc[ist]
                if isinstance(station_proc, (list, tuple)) and job_class < len(station_proc):
                    return station_proc[job_class]
                elif isinstance(station_proc, dict) and job_class in station_proc:
                    return station_proc[job_class]
    except (IndexError, KeyError, TypeError):
        pass
    return None


def _get_lld(lldscaling, ist, ni_val, lldlimit):
    """Get load-level-dependent scaling factor."""
    if lldscaling is None or lldscaling.size == 0:
        return 1.0
    idx = int(min(max(ni_val, 1), lldlimit)) - 1
    if isinstance(lldscaling, np.ndarray) and lldscaling.ndim >= 2:
        if ist < lldscaling.shape[0] and idx < lldscaling.shape[1]:
            return float(lldscaling[ist, idx])
    return 1.0


def _get_cd(cdscaling, ist, nir_row, job_class):
    """Get class-dependent scaling factor for a class job_class completion.

    sn.cdscaling is stored as a per-station list ([None, fn, ...]) — matching
    the flat builder which indexes cdscaling[ist] — but may also appear as a
    dict keyed by station index. Handle both.

    The handle maps the per-class population vector n (one state row) to the
    1xR vector of rate scalings beta_r(n) (see fes_beta_handle); the factor is
    the component of the completing class, v[min(job_class, len(v)-1)], so a
    scalar-returning handle applies to every class. Mirrors
    State.cdclassfactor.m.
    """
    if cdscaling is None:
        return 1.0
    fn = None
    if isinstance(cdscaling, dict):
        fn = cdscaling.get(ist, None)
    elif ist < len(cdscaling):
        fn = cdscaling[ist]
    if callable(fn):
        v = np.atleast_1d(np.asarray(
            fn(np.asarray(nir_row, dtype=float).flatten()), dtype=float)).flatten()
        return float(v[min(int(job_class), v.size - 1)])
    return 1.0


def _get_sched_weights(sn, ist, R):
    """Get scheduling weights (for DPS, GPS)."""
    w = np.ones(R)
    if hasattr(sn, 'schedparam') and sn.schedparam is not None:
        for r in range(R):
            if sn.schedparam[ist, r] is not None:
                w[r] = float(sn.schedparam[ist, r])
    w_sum = np.sum(w)
    if w_sum > 0:
        w = w / w_sum
    return w


def _get_niprio(sn, ist, nir_row, job_class, R):
    """Get total jobs at the same or higher priority level."""
    if hasattr(sn, 'classprio') and sn.classprio is not None:
        my_prio = sn.classprio[job_class]
        niprio = 0
        for r in range(R):
            if sn.classprio[r] <= my_prio and nir_row[r] > 0:
                niprio += nir_row[r]
        return niprio
    return np.sum(nir_row)


def _priority_nir(sn, nir_row, job_class, R):
    """Priority masking for the *PRIO PS/DPS/GPS disciplines.

    Mirrors MATLAB State.afterEventStation (PSPRIO/DPSPRIO/GPSPRIO cases): a
    job is served only if its class is in the most urgent priority group among
    classes currently present (lower classprio value = higher priority). When
    a strictly higher-priority class is present, this class is not served.

    Returns (eligible, nirprio):
      eligible  - False when a strictly higher-priority class is present, in
                  which case the departure rate must be 0.
      nirprio   - copy of nir_row with all classes not at job_class's priority
                  level zeroed out (the cohort that shares the server).
    """
    nir = np.asarray(nir_row, dtype=float)
    classprio = getattr(sn, 'classprio', None)
    if classprio is None:
        return True, nir.copy()
    my_prio = classprio[job_class]
    present = [classprio[r] for r in range(R) if nir[r] > 0]
    if not present:
        return True, np.zeros(R)
    if my_prio != min(present):
        return False, np.zeros(R)
    nirprio = np.array([nir[r] if classprio[r] == my_prio else 0.0
                        for r in range(R)], dtype=float)
    return True, nirprio


def _maybe_update_rrobin(sn, ind, job_class, R, space_var):
    """Update round-robin pointer if applicable."""
    if space_var.size == 0:
        return
    routing_val = sn.routing[ind, job_class] if hasattr(sn, 'routing') and sn.routing is not None else None
    if routing_val is None:
        return

    rs_val = int(routing_val.value) if hasattr(routing_val, 'value') else int(routing_val)

    # Round-robin / weighted round-robin are RoutingStrategy values.
    from ...constants import RoutingStrategy
    rr_val = int(RoutingStrategy.RROBIN.value) if hasattr(RoutingStrategy.RROBIN, 'value') else int(RoutingStrategy.RROBIN)
    wrr_val = int(RoutingStrategy.WRROBIN.value) if hasattr(RoutingStrategy.WRROBIN, 'value') else int(RoutingStrategy.WRROBIN)
    if rs_val not in (rr_val, wrr_val):
        return

    def _isrr(rt):
        v = int(rt.value) if hasattr(rt, 'value') else (int(rt) if rt is not None else -1)
        return v in (rr_val, wrr_val)

    # RR pointers after MAP phase vars in local-var block; a head write corrupts a service phase, freezes the pointer, degrading RR to static routing.
    from .ctmc_ssg import _map_phases_at
    _n_map = sum(1 for _r in range(R) if _map_phases_at(sn, ind, _r) > 1)

    # WRROBIN advances a cyclic POSITION pointer through the weighted-outlink cycle, not a node index.
    if rs_val == wrr_val:
        from .routing_pointer import wrr_weighted_outlinks
        wol = wrr_weighted_outlinks(sn, ind, job_class)
        if wol is not None and len(wol) > 0:
            cyc = len(wol)
            col = _n_map + sum(1 for rr in range(job_class) if _isrr(sn.routing[ind, rr]))
            ncol = space_var.shape[1] if space_var.ndim >= 2 else space_var.shape[0]
            if col >= ncol:
                return
            for row in range(space_var.shape[0] if space_var.ndim >= 2 else 1):
                cur = space_var[row, col] if space_var.ndim >= 2 else space_var[col]
                nxt = cur + 1 if 1 <= cur < cyc else 1
                if space_var.ndim >= 2:
                    space_var[row, col] = nxt
                else:
                    space_var[col] = nxt
            return

    # Outlinks: prefer nodeparam, else derive from the connection matrix (the
    # native path does not populate nodeparam.outlinks).
    ol = None
    nparam = sn.nodeparam[ind] if sn.nodeparam is not None and ind in sn.nodeparam else None
    if nparam is not None:
        outlinks = nparam.get('outlinks', None) if isinstance(nparam, dict) else getattr(nparam, 'outlinks', None)
        if outlinks is not None and job_class < len(outlinks):
            ol = np.atleast_1d(outlinks[job_class])
    if ol is None or len(ol) == 0:
        if hasattr(sn, 'connmatrix') and sn.connmatrix is not None:
            cm = np.asarray(sn.connmatrix)
            if ind < cm.shape[0]:
                ol = np.where(cm[ind, :] > 0)[0]
    if ol is None or len(ol) == 0:
        return

    # Pointer column within the variable block: one pointer per round-robin
    # class, ordered by class index.
    def _isrr(rt):
        v = int(rt.value) if hasattr(rt, 'value') else (int(rt) if rt is not None else -1)
        return v in (rr_val, wrr_val)
    space_var_idx = _n_map + sum(1 for rr in range(job_class) if _isrr(sn.routing[ind, rr]))
    ncol = space_var.shape[1] if space_var.ndim >= 2 else space_var.shape[0]
    if space_var_idx >= ncol:
        return

    for row in range(space_var.shape[0] if space_var.ndim >= 2 else 1):
        current_val = space_var[row, space_var_idx] if space_var.ndim >= 2 else space_var[space_var_idx]
        idx = -1
        for i in range(len(ol)):
            if ol[i] == current_val:
                idx = i
                break
        next_val = ol[idx + 1] if idx >= 0 and idx < len(ol) - 1 else ol[0]
        if space_var.ndim >= 2:
            space_var[row, space_var_idx] = next_val
        else:
            space_var[space_var_idx] = next_val


# ---------------------------------------------------------------------------
# Pass-and-swap (PAS) / order-independent station handler
# ---------------------------------------------------------------------------

def passAndSwap(c, p, G):
    """Apply the pass-and-swap mechanism (Dorsman & Gardner 2024, Sect. 2.3).

    Args:
        c: 1D sequence of 1-based class indices (ordered list, c[0] oldest).
        p: 0-based position of the job whose service token completes.
        G: (nclasses x nclasses) 0-based swapping-graph adjacency.

    Returns:
        (cnew, dep_class, chain) where cnew is the ordered list after the
        transition (1-based, length len(c)-1), dep_class is the 1-based class of
        the departing job, and chain is the list of visited positions.
    """
    n = len(c)
    chain = [p]
    moving = int(c[p])          # 1-based class currently scanning
    cur = p
    while True:
        q = -1
        for j in range(cur + 1, n):
            if G[moving - 1, int(c[j]) - 1]:
                q = j
                break
        if q < 0:
            break               # no swappable successor: moving class departs
        chain.append(q)
        moving = int(c[q])
        cur = q
    dep_class = int(c[chain[-1]])
    cnew = list(c)
    for i in range(len(chain) - 1):
        cnew[chain[i + 1]] = c[chain[i]]
    del cnew[chain[0]]
    return cnew, dep_class, chain


def after_event_station_pas(sn, ind, ist, inspace, event, job_class, R, V):
    """Event handler for pass-and-swap (PAS) / order-independent stations.

    The local state is the ordered list of 1-based class indices (oldest first),
    left-aligned in the first W = cap columns and right zero-padded; the trailing
    V columns are routing variables (carried through). Service is governed by the
    total rate function mu(c) (nodeparam svcRateFun, called with 0-based class
    indices) and the swapping graph G (nodeparam swapGraph, 0-based).
    """
    inspace = np.atleast_2d(inspace).astype(float)
    npj = sn.nodeparam[ind] if (sn.nodeparam is not None and ind in sn.nodeparam) else None
    mu_fun = npj.get('svcRateFun') if isinstance(npj, dict) else None
    G = npj.get('swapGraph') if isinstance(npj, dict) else None
    if mu_fun is None:
        raise ValueError('PAS station has no service rate function mu(c); set it via set_service(lambda c: ...).')

    n_cols = inspace.shape[1]
    W = n_cols - V
    cap = int(sn.cap[ist])
    n_rows = inspace.shape[0]
    job_class_1 = job_class + 1   # 1-based departing class to match list encoding

    out_space = []
    out_rate = []
    out_prob = []

    for row in range(n_rows):
        listrow = inspace[row, :W]
        varrow = inspace[row, W:]
        c = listrow[listrow > 0].astype(int)   # 1-based ordered list (contiguous)
        n = len(c)
        if event == EventType.ARV:
            if n >= cap:
                continue                        # buffer full: arrival lost
            newc = np.concatenate([c, [job_class_1]])
            padded = np.concatenate([newc.astype(float), np.zeros(W - len(newc))])
            out_space.append(np.concatenate([padded, varrow]))
            out_rate.append(-1.0)               # passive
            out_prob.append(1.0)
        elif event == EventType.DEP:
            if n == 0:
                continue
            c0 = c - 1                          # 0-based for mu(c) and G
            mu_prev = 0.0
            for p in range(n):
                mu_cur = float(mu_fun(c0[:p + 1]))
                ratep = mu_cur - mu_prev
                mu_prev = mu_cur
                if ratep <= 0:
                    continue
                cnew, dep_class, _ = passAndSwap(c, p, G)
                if dep_class != job_class_1:
                    continue
                padded = np.concatenate([np.array(cnew, dtype=float), np.zeros(W - len(cnew))])
                out_space.append(np.concatenate([padded, varrow]))
                out_rate.append(ratep)
                out_prob.append(1.0)
        # PHASE: PAS service is exponential, no phase transitions

    if out_space:
        outspace = np.array(out_space, dtype=float)
        outrate = np.array(out_rate, dtype=float).reshape(-1, 1)
        outprob = np.array(out_prob, dtype=float).reshape(-1, 1)
    else:
        outspace = np.zeros((0, n_cols))
        outrate = np.zeros((0, 1))
        outprob = np.zeros((0, 1))
    return outspace, outrate, outprob
