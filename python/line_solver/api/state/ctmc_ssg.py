"""
CTMC State Space Generator for sync-action-based CTMC builder.

Builds per-node state spaces (sn.space) and global state space with hashing.
Port from MATLAB ctmc_ssg.m + spaceGeneratorNodes.m + spaceGenerator.m.
"""

import warnings
import numpy as np
from typing import Tuple, Optional, Dict, List
from itertools import product as iproduct, islice
from .marginal import fromMarginal, toMarginal
from .polling import polling_space as _polling_space, polling_width as _polling_width
from .routing_pointer import wrr_weighted_outlinks as _wrr_weighted_outlinks
from ...lang.base import NodeType, SchedStrategy, RoutingStrategy
from ...constants import GlobalConstants, ProcessType


def _map_phases_at(sn, ind, r):
    """Number of MAP/MMPP modulating phases for class r at station node ind,
    or 0 if the service process there is not a MAP (correlated) process.

    A MAP service process needs a persistent server-phase state variable that a
    plain phase-type process does not, because its restart phase on completion
    is correlated with the completing phase (D1 is a full matrix, not d*alpha).
    This correlation only manifests at a single shared server (FCFS); for
    per-job disciplines (PS/INF/priority) MAP behaves as the marginal PH and no
    persistent phase variable is needed (matches MATLAB afterEventStation.m,
    where the INF/PS DEP uses mu*phi and leaves space_var unchanged).
    """
    if not (hasattr(sn, 'isstation') and ind < len(sn.isstation) and sn.isstation[ind]):
        return 0
    if not hasattr(sn, 'procid') or sn.procid is None:
        return 0
    ist = int(sn.nodeToStation[ind]) if hasattr(sn, 'nodeToStation') else ind
    # Only FCFS carries a persistent server MAP phase.
    if hasattr(sn, 'sched') and sn.sched is not None:
        sched_ist = sn.sched.get(ist) if isinstance(sn.sched, dict) else (
            sn.sched[ist] if ist < len(sn.sched) else None)
        sv = int(sched_ist.value) if hasattr(sched_ist, 'value') else (
            int(sched_ist) if sched_ist is not None else -1)
        fcfs_v = int(SchedStrategy.FCFS.value) if hasattr(SchedStrategy.FCFS, 'value') else int(SchedStrategy.FCFS)
        if sv != fcfs_v:
            return 0
    try:
        pid = sn.procid[ist, r]
    except (TypeError, KeyError, IndexError):
        return 0
    if pid not in (ProcessType.MAP, ProcessType.MMPP2):
        return 0
    # phase count from the D0 block of the (D0, D1) representation
    try:
        proc_ir = sn.proc[ist][r] if not isinstance(sn.proc, dict) else sn.proc[ist][r]
        D0 = np.atleast_2d(np.asarray(proc_ir[0], dtype=float))
        return int(D0.shape[0])
    except (TypeError, KeyError, IndexError, AttributeError):
        return 0


def _assert_no_nonfcfs_map(sn):
    """Mirror MATLAB State.fromMarginal (line 39): MAP/MMPP service is only
    supported at FCFS stations, because the correlated D1 restart needs a single
    shared server with a persistent phase. Raise for any non-FCFS MAP station so
    the sync builder rejects the model exactly as MATLAB does, rather than
    silently producing an unvalidated result.
    """
    if not hasattr(sn, 'procid') or sn.procid is None:
        return
    fcfs_v = int(SchedStrategy.FCFS.value) if hasattr(SchedStrategy.FCFS, 'value') else int(SchedStrategy.FCFS)
    ext_v = int(SchedStrategy.EXT.value) if hasattr(SchedStrategy.EXT, 'value') else int(SchedStrategy.EXT)
    R = int(sn.nclasses) if hasattr(sn, 'nclasses') else 0
    for ind in range(int(sn.nnodes)):
        if not (hasattr(sn, 'isstation') and ind < len(sn.isstation) and sn.isstation[ind]):
            continue
        ist = int(sn.nodeToStation[ind]) if hasattr(sn, 'nodeToStation') else ind
        sched_ist = None
        if hasattr(sn, 'sched') and sn.sched is not None:
            sched_ist = sn.sched.get(ist) if isinstance(sn.sched, dict) else (
                sn.sched[ist] if ist < len(sn.sched) else None)
        sv = int(sched_ist.value) if hasattr(sched_ist, 'value') else (
            int(sched_ist) if sched_ist is not None else -1)
        # FCFS stations are supported; Source (EXT) MAP is an arrival process,
        # also allowed (MATLAB: sn.nodetype ~= NodeType.Source).
        if sv in (fcfs_v, ext_v):
            continue
        for r in range(R):
            try:
                pid = sn.procid[ist, r]
            except (TypeError, KeyError, IndexError):
                continue
            if pid in (ProcessType.MAP, ProcessType.MMPP2):
                raise NotImplementedError("Non-FCFS MAP stations are not supported.")


def _downstream_stations(sn, ind):
    """Node indices of the stations directly downstream of node ind, walking through
    stateless nodes (Router, ClassSwitch, ...) but stopping at the first station on
    each path, since a blocked job is held for its immediate destination. Mirrors the
    MATLAB/JAR downstreamStations."""
    conn = getattr(sn, 'connmatrix', None)
    if conn is None:
        return []
    conn = np.atleast_2d(conn)
    n = min(sn.nnodes, conn.shape[1])
    out = []
    seen = np.zeros(sn.nnodes, dtype=bool)
    queue = [ind]
    seen[ind] = True
    while queue:
        cur = queue.pop(0)
        for j in range(n):
            if conn[cur, j] != 1 or seen[j]:
                continue
            seen[j] = True
            if bool(sn.isstation[j]):
                out.append(j)   # immediate destination station
            else:
                queue.append(j)  # stateless hop, keep walking
    return out


def _declares_blocked_marker(sn, ind):
    """(tf, destmask) for the station at node ind.

    TF is True when the station is the BLOCKING (upstream) side of a true-BAS
    relation: it has a directly reachable destination station with finite capacity, and
    BAS is declared EITHER on this station (upstream form) OR on that destination
    (destination form, canonical since BUG-83). Both resolve to the same blocking
    station. DESTMASK is an (nstations x nclasses) boolean marking the (destination
    station, class) pairs whose refusals must block IND rather than be lost.
    Mirrors MATLAB/JAR declaresBlockedMarker."""
    destmask = np.zeros((int(sn.nstations), int(sn.nclasses)), dtype=bool)
    dr = getattr(sn, 'droprule', None)
    if dr is None or not bool(sn.isstation[ind]) or not _is_blockable_station_node(sn, ind):
        return False, destmask
    ist = int(sn.nodeToStation[ind])
    if ist < 0:
        return False, destmask
    tf = False
    dests = _downstream_stations(sn, ind)
    for r in range(sn.nclasses):
        here_bas = False
        try:
            here_bas = int(dr[ist, r]) == 2
        except (IndexError, TypeError, ValueError, KeyError):
            pass
        for dnode in dests:
            jst = int(sn.nodeToStation[dnode])
            if jst < 0:
                continue
            try:
                dcap = float(np.asarray(sn.cap).flatten()[jst])
            except (IndexError, TypeError, ValueError):
                continue
            if not (np.isfinite(dcap) and dcap > 0):
                continue
            there_bas = False
            try:
                there_bas = int(dr[jst, r]) == 2
            except (IndexError, TypeError, ValueError, KeyError):
                pass
            if here_bas or there_bas:
                # BAS declared upstream or on the reachable full destination. Do not
                # stop here: every such destination must be recorded, since a refusal
                # at any of them has to block IND instead of dropping.
                tf = True
                destmask[jst, r] = True
    return tf, destmask


def _is_blockable_station_node(sn, ind):
    """True when node ind is a Station that can carry a BAS blocked marker, i.e. not a
    Source and not a Cache. Mirrors the isa(node,'Station') && ~isa(node,'Source') &&
    ~isa(node,'Cache') guard in MATLAB refreshLocalVars / JAR Network.refreshLocalVars."""
    from ..sn.network_struct import NodeType as _NT
    nt = getattr(sn, 'nodetype', None)
    if nt is None:
        return True
    try:
        code = int(nt[ind])
    except (IndexError, TypeError, ValueError):
        return True
    return code != int(_NT.SOURCE) and code != int(_NT.CACHE)


def _populate_isbasblocking(sn):
    """Set sn.isbasblocking (node-indexed, 1 iff blocking side of a true-BAS relation)
    and sn.isbasdestination ((nstations,nclasses), True on the RECEIVING side).

    isbasdestination is the complement of isbasblocking: a refusal is resolved at the
    destination by _arrival_is_lost, which can only see that station's own drop rule,
    but under the upstream declaration form the BAS rule sits on the blocking station.
    Without the field an open-class arrival refused at the destination is declared
    LOST, the become-blocked edge never fires, and the upstream degenerates into an
    isolated M/M/1/K while the destination silently drops its overflow.

    Rejects a polling+BAS station (Option B, BUG-83): the polling controller and the
    BAS marker share one local-state column and a station cannot carry both."""
    from ...lang.base import SchedStrategy as _Sched
    isb = np.zeros(sn.nnodes, dtype=int)
    isbasdest = np.zeros((int(sn.nstations), int(sn.nclasses)), dtype=bool)
    for ind in range(sn.nnodes):
        if not bool(sn.isstation[ind]):
            continue
        declares, destmask = _declares_blocked_marker(sn, ind)
        isbasdest |= destmask
        if declares:
            ist = int(sn.nodeToStation[ind])
            try:
                is_polling = int(sn.sched[ist]) == int(_Sched.POLLING)
            except (IndexError, TypeError, ValueError):
                is_polling = False
            if is_polling:
                name = sn.nodenames[ind] if getattr(sn, 'nodenames', None) is not None else str(ind)
                raise RuntimeError(
                    "True BAS blocking is not supported at a polling station (%s): the polling "
                    "controller and the BAS blocked marker share one local-state column. Use a "
                    "non-polling scheduling strategy at the blocking station, or remove the BAS "
                    "drop rule." % name)
            isb[ind] = 1
    sn.isbasblocking = isb
    sn.isbasdestination = isbasdest


def _is_breakdown_station(sn, ind):
    """True when node ind carries a server subject to breakdowns (sn.hasbreakdown)."""
    hb = getattr(sn, 'hasbreakdown', None)
    if hb is None:
        return False
    hb = np.asarray(hb).ravel()
    return ind < hb.size and int(hb[ind]) == 1


def _assert_breakdown_exclusive(sn, ind, ist):
    """Reject a breakdown station that also needs the shared trailing
    local-state column.

    Like the BAS marker and the polling controller, the breakdown status
    occupies the shared trailing local-variable column, so a station cannot
    carry two of them at once. Rejecting the combination here, before the state
    space is built, is what keeps the shared column from being silently
    corrupted."""
    from ...lang.base import SchedStrategy as _Sched
    name = sn.nodenames[ind] if getattr(sn, 'nodenames', None) is not None else str(ind)
    if _is_bas_station(sn, ist):
        raise RuntimeError(
            "Station '%s' combines server breakdowns with true-BAS blocking: the breakdown "
            "status and the BAS blocked marker share one local-state column. Remove the BAS "
            "drop rule or the breakdown." % name)
    try:
        is_polling = int(sn.sched[ist]) == int(_Sched.POLLING)
    except (IndexError, TypeError, ValueError):
        is_polling = False
    if is_polling:
        raise RuntimeError(
            "Server breakdowns are not supported at a polling station (%s): the polling "
            "controller and the breakdown status share one local-state column." % name)


def _is_bas_station(sn, ist):
    """True when station ist carries a true-BAS blocked marker. Consults the dedicated
    sn.isbasblocking field (node-indexed, set by _populate_isbasblocking under BOTH the
    upstream and destination declaration forms) so declaration and enumeration agree.
    Falls back to the own-station rule only if the field is absent."""
    isb = getattr(sn, 'isbasblocking', None)
    if isb is not None:
        try:
            ind = int(sn.stationToNode[ist])
            return ind >= 0 and int(np.asarray(isb).flatten()[ind]) == 1
        except (IndexError, TypeError, ValueError):
            return False
    # Fallback (field not built): the pre-BUG-83 own-station test.
    dr = getattr(sn, 'droprule', None)
    cap = getattr(sn, 'cap', None)
    if dr is None or cap is None:
        return False
    try:
        capv = float(np.asarray(cap).flatten()[ist])
    except (IndexError, TypeError, ValueError):
        return False
    if not np.isfinite(capv) or capv <= 0:
        return False
    for r in range(sn.nclasses):
        try:
            if int(dr[ist, r]) == 2:  # DropStrategy.BAS
                return True
        except (IndexError, TypeError, ValueError, KeyError):
            pass
    return False


def ctmc_ssg(sn, cutoff, options=None):
    """
    Build per-node state spaces and global hashed state space for CTMC.

    Port from MATLAB ctmc_ssg.m.

    Args:
        sn: NetworkStruct
        cutoff: Population cutoff for open classes (scalar or matrix)
        options: Optional solver options dict

    Returns:
        Tuple of (state_space, state_space_aggr, state_space_hashed, sn):
        - state_space: Full state space (rows = global states, cols = concatenated per-node states)
        - state_space_aggr: Aggregated state space (rows = states, cols = M*K per-station per-class counts)
        - state_space_hashed: Hashed state space (rows = states, cols = nstateful, values = row indices in sn.space)
        - sn: Updated NetworkStruct with sn.space populated
    """
    if options is None:
        options = {}

    _assert_no_nonfcfs_map(sn)

    R = sn.nclasses
    M = sn.nstations

    # Expand cutoff to matrix
    if np.isscalar(cutoff):
        cutoff_matrix = np.full((M, R), cutoff, dtype=float)
    else:
        cutoff_matrix = np.atleast_2d(cutoff).astype(float)
        if cutoff_matrix.shape[0] == 1 and M > 1:
            cutoff_matrix = np.tile(cutoff_matrix, (M, 1))

    # Step 1: Build per-node state spaces
    capacityc = _space_generator_nodes(sn, cutoff_matrix, options)

    # Step 2: Build global state space from per-node spaces
    N = sn.njobs.flatten() if sn.njobs is not None else np.ones(R)
    Np = N.copy()
    is_open_class = np.isinf(Np)
    for r in range(R):
        if is_open_class[r]:
            Np[r] = np.max(capacityc[:, r])

    # State-space safeguard: read the cap once here and thread it into the
    # builders, which abort as soon as the ACTUAL number of enumerated global
    # states crosses it. A product-of-local-spaces estimate is unusable (it
    # grossly overcounts closed models, whose reachable states are a tiny
    # fraction of the local-space product), so the guard is on the real row count
    # instead -- no false positives. A model MVA/MAM cannot handle (e.g. a
    # finite-capacity region with a large token bound) can enumerate tens of
    # millions of states; materializing the generator then exhausts memory and
    # OOM-kills the process. Configurable via options.ctmc_max_states (default
    # 3e6); None/<=0 disables the guard.
    _sscap = options.get('ctmc_max_states', 3_000_000) if isinstance(options, dict) \
        else getattr(options, 'ctmc_max_states', 3_000_000)
    if _sscap is None or _sscap <= 0:
        _sscap = float('inf')

    state_space, state_space_hashed = _build_global_states(sn, Np, is_open_class, capacityc,
                                                           options=options, max_states=_sscap)

    # Step 3: Build aggregated state space (per-station per-class job counts)
    state_space_aggr = _build_state_space_aggr(sn, state_space_hashed)

    return state_space, state_space_aggr, state_space_hashed, sn


def _space_generator_nodes(sn, cutoff_matrix, options):
    """
    Build per-node state spaces sn.space[isf] for each stateful node.

    Port from MATLAB spaceGeneratorNodes.m.

    Args:
        sn: NetworkStruct
        cutoff_matrix: Cutoff matrix (M x R)
        options: Solver options

    Returns:
        capacityc: Per-node per-class capacity matrix (nnodes x R)
    """
    R = sn.nclasses
    M = sn.nstations
    N = sn.njobs.flatten() if sn.njobs is not None else np.ones(R)

    # MATLAB's spaceGeneratorNodes.m receives sn by value, so its local
    # substitution of inf nservers (for IS stations) by a finite job-count bound
    # never escapes to the caller's model. Python passes the NetworkStruct by
    # reference, so we snapshot and restore nservers to confine that mutation to
    # this function -- otherwise a later solver on the same model (e.g. MVA after
    # SSA) would see a finite IS nservers and wrongly divide IS utilization.
    import copy as _copy
    _nservers_orig = _copy.deepcopy(sn.nservers)

    sn.space = {}
    capacityc = np.zeros((sn.nnodes, R))
    # Per-node local-variable width (sum over the row gives V in after_event,
    # which slices the trailing local-var columns off the server/buffer state).
    # getStruct does not populate sn.nvars for the native path, so build it here
    # from the local-var spaces; without it the cache server/var split is wrong.
    sn.nvars = np.zeros((sn.nnodes, 1), dtype=int)
    _populate_isbasblocking(sn)

    for ind in range(sn.nnodes):
        if sn.isstation[ind]:
            ist = int(sn.nodeToStation[ind])
            isf = int(sn.nodeToStateful[ind])

            for r in range(R):
                c = _find_chain(sn, r)

                # Check visits
                if c is not None and sn.visits is not None and c in sn.visits:
                    v = sn.visits[c]
                    # Use a tolerance: visit ratios for non-visited (station,class)
                    # pairs can carry tiny floating-point noise (~1e-17) rather than
                    # an exact 0, which an `== 0` test would miss — leaving the class
                    # erroneously enabled and over-generating the state space.
                    if v is not None and isf < v.shape[0] and abs(v[isf, r]) < 1e-12:
                        capacityc[ind, r] = 0
                        continue

                # Check disabled process — Places (SPN buffers) hold tokens
                # without a service distribution, so their procid is Disabled
                # by construction. MATLAB spaceGeneratorNodes.m excludes
                # Places from this branch (line 25 of the .m). Mirror that.
                _is_place_node = False
                if hasattr(sn, 'nodetype') and sn.nodetype is not None and ind < len(sn.nodetype):
                    _nt = sn.nodetype[ind]
                    _nt_val = int(_nt.value) if hasattr(_nt, 'value') else int(_nt)
                    _place_val = int(NodeType.PLACE.value) if hasattr(NodeType.PLACE, 'value') else int(NodeType.PLACE)
                    _is_place_node = (_nt_val == _place_val)
                if _is_disabled_proc(sn, ist, r) and not _is_place_node:
                    capacityc[ind, r] = 0
                    continue

                if np.isinf(N[r]):
                    cc = cutoff_matrix[ist, r] if ist < cutoff_matrix.shape[0] else cutoff_matrix[0, r]
                    capacityc[ind, r] = min(cc, sn.classcap[ist, r])
                else:
                    # closed classes: enumerate up to the chain population, but never
                    # beyond the class capacity at this station (finite-buffer stations)
                    if c is not None:
                        chain_mask = sn.chains[c].astype(bool)
                        capacityc[ind, r] = min(np.sum(N[chain_mask]), sn.classcap[ist, r])
                    else:
                        capacityc[ind, r] = min(N[r], sn.classcap[ist, r])

                # Finite-capacity region bound: a station in a DROP/WAITQ region can
                # never hold more than the region's per-class cap (nor the region-global
                # cap) of class r, so enumerating beyond it produces only states the
                # region filter later discards. Bounding capacityc here keeps the
                # generated space small (the auto-cutoff is region-blind and can far
                # exceed any reachable population), with an identical final result.
                # -1 = unbounded. sn.region[f] is M x (K+1): cols 0..K-1 per-class,
                # col K the region-global cap at each member station. Mirrors MATLAB
                # State.spaceGeneratorNodes.
                if getattr(sn, 'nregions', 0) and capacityc[ind, r] > 0:
                    for f in range(sn.nregions):
                        regf = sn.region[f]
                        if ist < regf.shape[0]:
                            rc = regf[ist, r]
                            gc = regf[ist, R]
                            if rc >= 0:
                                capacityc[ind, r] = min(capacityc[ind, r], rc)
                            if gc >= 0:
                                capacityc[ind, r] = min(capacityc[ind, r], gc)

            # Generate per-node state space including local variables
            # Keep infinite capacity for Source (EXT) stations
            if np.isfinite(sn.cap[ist]):
                cap_val = int(sn.cap[ist])
            else:
                sched = sn.sched[ist] if hasattr(sn, 'sched') else None
                ext_val = int(SchedStrategy.EXT.value) if hasattr(SchedStrategy.EXT, 'value') else int(SchedStrategy.EXT)
                sched_val = int(sched.value) if hasattr(sched, 'value') else int(sched) if sched is not None else -1
                if sched_val == ext_val:
                    cap_val = float('inf')  # Source: keep infinite capacity
                else:
                    cap_val = int(np.sum(capacityc[ind]))
            if not np.isfinite(cap_val):
                # Source (EXT): a Source has no finite per-class job counts, so
                # the marginal-bounds enumeration is ill-defined and collapses to
                # a wrong degenerate state when a closed class is present (mixed
                # models). Use the canonical external state instead — exactly the
                # state the global combiner looks up via fromMarginal(sn, ind, [])
                # — so the enumerated Source row matches and the lattice is
                # non-empty.
                sn.space[isf] = np.atleast_2d(fromMarginal(sn, ind, []))
            else:
                sn.space[isf] = _from_marginal_bounds(
                    sn, ind, capacityc[ind].astype(int), cap_val, options)

            # Persistent MAP server-phase augmentation. A MAP service process
            # keeps evolving its modulating phase and, on completion, restarts in
            # a phase correlated with the completing one. While a job of class r
            # is in service its phase is already encoded by the service phase
            # counts; when the server is idle the phase must still be remembered.
            # Add one phase variable per MAP class, pinned to the in-service phase
            # while busy and free (every phase) when idle, so the state count
            # matches the monolithic builder.
            _nvar_acc = 0
            _map_classes = [(r, _map_phases_at(sn, ind, r)) for r in range(R)]
            _map_classes = [(r, nph) for r, nph in _map_classes if nph > 1]
            if _map_classes:
                # Compute the in-service phase per MAP class ONCE on the clean
                # [buffer|server] state (before appending any var column, so
                # toMarginal reads the layout correctly), then append one phase
                # variable per MAP class via a cartesian product: pinned to the
                # in-service phase while busy, free across all phases when idle.
                base = np.atleast_2d(sn.space[isf])
                _, _, _, kir = toMarginal(sn, ind, base)
                kir = np.asarray(kir)
                new_rows = []
                for ri in range(base.shape[0]):
                    per_class = []
                    for r, nph in _map_classes:
                        kr = kir[ri, r, :nph] if kir.ndim == 3 else np.zeros(nph)
                        if np.sum(kr) > 0:
                            per_class.append([int(np.argmax(kr))])
                        else:
                            per_class.append(list(range(nph)))
                    for combo in iproduct(*per_class):
                        new_rows.append(np.hstack([base[ri], list(combo)]))
                sn.space[isf] = np.array(new_rows, dtype=int)
                _nvar_acc += len(_map_classes)

            # Append round-robin routing pointer local variables: a station that
            # routes via RROBIN/WRROBIN tracks its next-outlink pointer in the
            # per-node state, exactly like the non-station stateful branch below.
            _rr_var = _space_local_vars(sn, ind)
            if _rr_var is not None and np.atleast_2d(_rr_var).size > 0:
                _nvar_acc += int(np.atleast_2d(_rr_var).shape[1])
                sn.space[isf] = _cartesian_2d(sn.space[isf], _rr_var)

            # True BAS: a finite-capacity station using Blocking-After-Service holds a
            # completed job at its server (blocked) until the downstream has room. Append a
            # {0,1} blocked marker as the last local var — enumerated only for states holding
            # >=1 job (empty station -> blocked=0). Driven by after_event_station (DEP-to-full
            # -> blocked; a blocked front job departs at rate 1e7 once the dest has room) and
            # the CTMC generator become-blocked step.
            if _is_bas_station(sn, ist):
                cur = np.atleast_2d(sn.space[isf])
                new_rows = []
                for ri in range(cur.shape[0]):
                    row = cur[ri]
                    if np.sum(row) > 0:
                        new_rows.append(np.hstack([row, 0]))
                        new_rows.append(np.hstack([row, 1]))
                    else:
                        new_rows.append(np.hstack([row, 0]))
                sn.space[isf] = np.array(new_rows, dtype=int)
                _nvar_acc += 1

            # Server breakdown status: 0 = down, 1 = up. Enumerated for EVERY
            # marginal INCLUDING the empty station, because the failure clock
            # runs whenever the server is up, idle or busy. Emitting only
            # status 0 for the empty state would drop the empty-and-up state
            # and the station would behave as if it could never empty, which
            # inflates its queue length by about one job. The column is
            # exclusive with the BAS marker and the polling controller
            # (_assert_breakdown_exclusive), so it is always the trailing one.
            if _is_breakdown_station(sn, ind):
                _assert_breakdown_exclusive(sn, ind, ist)
                cur = np.atleast_2d(sn.space[isf])
                sn.space[isf] = np.vstack([
                    np.hstack([cur, np.zeros((cur.shape[0], 1), dtype=cur.dtype)]),
                    np.hstack([cur, np.ones((cur.shape[0], 1), dtype=cur.dtype)]),
                ])
                _nvar_acc += 1

            # Polling controller [pos, swk, ctr]. Appended LAST, after every
            # other local variable: the round-robin pointer is located from the
            # head of the local-variable block (see lang/sync._make_rr_prob), so
            # the polling columns must trail it. Enumerated per row, since the
            # occupiable controller configurations depend on the buffer and
            # server contents of each state (polling.polling_space).
            _pw = _polling_width(sn, ind)
            if _pw > 0:
                _Kp = np.array(sn.phasessz[ist], dtype=int)
                _Ksp = np.array(sn.phaseshift[ist], dtype=int)
                sn.space[isf] = _polling_space(sn, ind, sn.space[isf], _Kp, _Ksp)
                _nvar_acc += _pw

            if _nvar_acc > 0:
                sn.nvars[ind, 0] = _nvar_acc

            if np.isinf(sn.nservers[ist]):
                sn.nservers[ist] = int(np.sum(capacityc[ind]))

        elif sn.isstateful[ind]:
            isf = int(sn.nodeToStateful[ind])
            nt_val = _get_nodetype_val(sn, ind)
            cache_val = int(NodeType.CACHE.value) if hasattr(NodeType.CACHE, 'value') else int(NodeType.CACHE)
            router_val = int(NodeType.ROUTER.value) if hasattr(NodeType.ROUTER, 'value') else int(NodeType.ROUTER)
            transition_val = int(NodeType.TRANSITION.value) if hasattr(NodeType.TRANSITION, 'value') else int(NodeType.TRANSITION)

            if nt_val == transition_val:
                # SPN Transition: per-mode state vector
                # [idle_count_m, phase_counts_m(fK_m)] cartesian over modes,
                # then nmodes fired-counts, then local vars. Mirrors
                # matlab/src/lang/+State/spaceGeneratorNodes.m:71-118.
                capacityc[ind, :] = 0
                sn.space[isf] = _build_transition_space(sn, ind, cutoff_matrix, options)
                continue

            if nt_val == cache_val:
                capacityc[ind, :] = 1
                # Cache local state = [per-class server presence | cache contents
                # + occupancy bitmap]. The server holds the single job currently
                # being read (at most one job at a time). fromMarginal stubs
                # non-station stateful nodes, so enumerate the per-class server
                # marginals directly (empty + one job of each enabled class),
                # mirroring MATLAB fromMarginalBounds(sn,ind,capacityc,1,options).
                serv_rows = [np.zeros(R, dtype=int)]
                for rr in range(R):
                    if capacityc[ind, rr] > 0:
                        sv = np.zeros(R, dtype=int)
                        sv[rr] = 1
                        serv_rows.append(sv)
                state_bufsrv = np.array(serv_rows, dtype=int)
                state_var = _space_local_vars(sn, ind)
                if state_var is not None and np.atleast_2d(state_var).size > 0:
                    sn.nvars[ind, 0] = int(np.atleast_2d(state_var).shape[1])
                    sn.space[isf] = _cartesian_2d(state_bufsrv, state_var)
                else:
                    sn.space[isf] = state_bufsrv
                sn.space[isf] = _filter_cache_states(
                    sn, ind, sn.space[isf], state_bufsrv.shape[1])
                continue
            elif nt_val == router_val:
                # A Router holds a job only transiently (immediate pass-through),
                # a state fromMarginal stubs for non-station nodes. Build the
                # per-class occupancy directly: empty, or one job of each class
                # that visits this router. These immediate states are folded out
                # by stochastic complementation downstream.
                serv_rows = [np.zeros(R, dtype=int)]
                for r in range(R):
                    c = _find_chain(sn, r)
                    if c is not None and sn.nodevisits is not None and c in sn.nodevisits:
                        nv = sn.nodevisits[c]
                        if nv is not None and ind < nv.shape[0] and nv[ind, r] > 0:
                            capacityc[ind, r] = 1
                            sv = np.zeros(R, dtype=int)
                            sv[r] = 1
                            serv_rows.append(sv)
                state_bufsrv = np.array(serv_rows, dtype=int)
                state_var = _space_local_vars(sn, ind)
                if state_var is not None and np.atleast_2d(state_var).size > 0:
                    sn.nvars[ind, 0] = int(np.atleast_2d(state_var).shape[1])
                    sn.space[isf] = _cartesian_2d(state_bufsrv, state_var)
                else:
                    sn.space[isf] = state_bufsrv
                continue
            else:
                capacityc[ind, :] = 1

            state_bufsrv = _from_marginal_bounds(
                sn, ind, capacityc[ind].astype(int), 1, options)
            state_var = _space_local_vars(sn, ind)

            if state_var is not None and np.atleast_2d(state_var).size > 0:
                sn.nvars[ind, 0] = int(np.atleast_2d(state_var).shape[1])
                sn.space[isf] = _cartesian_2d(state_bufsrv, state_var)
            else:
                sn.space[isf] = state_bufsrv

            if nt_val == cache_val:
                sn.space[isf] = _filter_cache_states(
                    sn, ind, sn.space[isf], state_bufsrv.shape[1])

    # Restore IS nservers so the finite enumeration bound does not leak to the
    # caller's model (see snapshot note above).
    sn.nservers = _nservers_orig
    return capacityc


def _filter_cache_states(sn, ind, value, lvs):
    """
    Prune unreachable cache states for a retrieval-system cache.

    The unfiltered cartesian product enumerates states the dynamics can never
    reach -- an item recorded in a retrieval slot with no matching retrieval
    job, a retrieval-pending and -complete job coexisting for one item, or a
    miss departing a full retrieval system. Leaving them in pollutes the CTMC
    stationary distribution. Non-retrieval caches are returned unchanged.

    Mirrors the MATLAB State.spaceGeneratorNodes cache validity filter.

    Args:
        sn: NetworkStruct
        ind: cache node index (0-based)
        value: cache state space [class counts | cache contents | retrieval slots]
        lvs: local-vars start column (width of the class-counts region)
    """
    nparam = sn.nodeparam[ind] if sn.nodeparam is not None and ind in sn.nodeparam else None
    if nparam is None or value is None or value.size == 0:
        return value
    rsc = nparam.get('retrieval_system_capacity', 0) if isinstance(nparam, dict) \
        else getattr(nparam, 'retrieval_system_capacity', 0)
    if rsc <= 0:
        return value
    n_items = nparam.get('nitems', 0) if isinstance(nparam, dict) else getattr(nparam, 'nitems', 0)
    tcc = nparam.get('total_cache_capacity', 0) if isinstance(nparam, dict) \
        else getattr(nparam, 'total_cache_capacity', 0)
    rc_mat = nparam.get('retrieval_classes', None) if isinstance(nparam, dict) \
        else getattr(nparam, 'retrieval_classes', None)
    missclass = nparam.get('missclass', None) if isinstance(nparam, dict) \
        else getattr(nparam, 'missclass', None)
    if rc_mat is None:
        return value
    rc_mat = np.atleast_2d(np.asarray(rc_mat, dtype=int))
    if rc_mat.size == 0:
        return value
    missclass = np.atleast_1d(np.asarray(missclass, dtype=int)) if missclass is not None else np.array([], dtype=int)

    keep = []
    for row in range(value.shape[0]):
        st = value[row]
        valid = np.sum(st[:lvs]) <= 1  # at most one class in the buffer/server region

        # Retrieval-system occupancy bitmap: column (lvs+tcc+item) is non-zero
        # iff that item is currently being retrieved.
        items_in_rs = set()
        for col in range(lvs + tcc, st.shape[0]):
            if st[col] != 0:
                items_in_rs.add(col - (lvs + tcc))  # 0-based item id

        # If an item is in the cache, retrievals for it cannot arrive or depart.
        if valid:
            for col in range(lvs, lvs + tcc):
                item = int(st[col]) - 1  # items are 1-indexed in the cache slots
                if item < 0 or item >= rc_mat.shape[0]:
                    continue
                for jac in range(rc_mat.shape[1]):
                    rc = int(rc_mat[item, jac])
                    if rc < 0:
                        continue
                    if st[rc] > 0:
                        valid = False
                        break
                if not valid:
                    break

        # Retrieval system at full capacity -> no miss may be ready to depart.
        if valid and len(items_in_rs) == (n_items - tcc):
            for col in range(missclass.shape[0]):
                jc = int(missclass[col])
                if jc < 0:
                    continue
                if st[jc] > 0:
                    valid = False
                    break

        if valid:
            keep.append(row)

    return value[keep] if keep else np.zeros((0, value.shape[1]))


def _multichoose(n_bins: int, total: int) -> np.ndarray:
    """
    Enumerate all non-negative integer vectors of length n_bins that sum to total.

    Mirrors MATLAB multichoose used by spaceGeneratorNodes / afterGlobalEvent.
    Returns an (n_combos, n_bins) array.
    """
    if n_bins <= 0:
        return np.zeros((1 if total == 0 else 0, 0), dtype=int)
    if n_bins == 1:
        return np.array([[total]], dtype=int)
    rows = []
    for i in range(total + 1):
        sub = _multichoose(n_bins - 1, total - i)
        if sub.shape[0] == 0:
            continue
        first = np.full((sub.shape[0], 1), i, dtype=int)
        rows.append(np.hstack([first, sub]))
    if not rows:
        return np.zeros((0, n_bins), dtype=int)
    return np.vstack(rows)


def _build_transition_space(sn, ind: int, cutoff_matrix, options) -> np.ndarray:
    """
    Build per-mode SPN Transition state space.

    Layout per row (matching MATLAB spaceGeneratorNodes.m:71-117):

        [ idle_m1, phase_counts_m1(fK_1),
          idle_m2, phase_counts_m2(fK_2),
          ...
          idle_mM, phase_counts_mM(fK_M),
          fired_m1, fired_m2, ..., fired_mM,
          <local_vars> ]

    Idle counts represent servers in the disabled pool; phase counts track
    servers running through each phase of the firing process. Fired counts
    are needed by simulation (initialised to zero at SSG; CTMC ignores them).
    """
    nparam = sn.nodeparam[ind]
    nmodes = int(getattr(nparam, 'nmodes', 0) if not isinstance(nparam, dict)
                 else nparam.get('nmodes', 0))
    if nmodes <= 0:
        return np.zeros((1, 0), dtype=float)

    firing_phases = np.asarray(getattr(nparam, 'firingphases', None), dtype=float)
    if firing_phases is None or firing_phases.size != nmodes:
        firing_phases = np.full(nmodes, np.nan, dtype=float)

    # Resolve NaN firingphases from firingproc shape (matches MATLAB lines 76-85).
    firing_proc = getattr(nparam, 'firingproc', None)
    fK = np.zeros(nmodes, dtype=int)
    for m in range(nmodes):
        if not np.isnan(firing_phases[m]):
            fK[m] = int(firing_phases[m])
        elif firing_proc is not None and m < len(firing_proc) and firing_proc[m] is not None:
            D0 = np.atleast_2d(np.asarray(firing_proc[m][0]))
            fK[m] = int(D0.shape[0])
        else:
            fK[m] = 1
    fK = np.maximum(fK, 1)

    nmodeservers = np.asarray(getattr(nparam, 'nmodeservers', np.ones(nmodes)), dtype=float)
    if nmodeservers.size != nmodes:
        nmodeservers = np.broadcast_to(nmodeservers, (nmodes,)).astype(float).copy()

    # Total tokens that can ever be enabled here: closed-class population +
    # per-station cutoffs for open classes. MATLAB:
    #     max_jobs = sum(njobs(~isinf)) + (any open ? sum(cutoff(1,:)) : 0)
    R = int(sn.nclasses)
    njobs = np.asarray(sn.njobs, dtype=float).flatten() if sn.njobs is not None else np.zeros(R)
    finite_mask = np.isfinite(njobs)
    max_jobs = float(np.sum(njobs[finite_mask])) if finite_mask.any() else 0.0
    if (~finite_mask).any():
        cutoff_row = np.atleast_2d(cutoff_matrix)[0]
        max_jobs += float(np.sum(cutoff_row))
    max_jobs = max(int(max_jobs), 0)

    MAX_INT = int(GlobalConstants.MaxInt) if hasattr(GlobalConstants, 'MaxInt') else (2 ** 31 - 1)

    mode_spaces = []
    for m in range(nmodes):
        max_srv_m = nmodeservers[m]
        if not np.isfinite(max_srv_m):
            max_srv_m = max_jobs
        max_srv_m = int(min(max_srv_m, max_jobs))
        rows_idle = []
        rows_phases = []
        # buf_full = nmodeservers[m] (finite) or MaxInt (infinite); idle = buf_full - total
        buf_full_m = MAX_INT if not np.isfinite(nmodeservers[m]) else int(nmodeservers[m])
        for total in range(max_srv_m + 1):
            phase_combs = _multichoose(int(fK[m]), total)
            if phase_combs.shape[0] == 0:
                continue
            idle = buf_full_m - total
            idle_col = np.full((phase_combs.shape[0], 1), idle, dtype=int)
            rows_idle.append(idle_col)
            rows_phases.append(phase_combs)
        if rows_idle:
            mode_idle = np.vstack(rows_idle).astype(float)
            mode_phases = np.vstack(rows_phases).astype(float)
        else:
            mode_idle = np.zeros((1, 1), dtype=float)
            mode_phases = np.zeros((1, int(fK[m])), dtype=float)
        mode_spaces.append((mode_idle, mode_phases))

    # Cartesian product across modes, but separate (idle, phases) so that the
    # final layout is [idle(nmodes), phases(sum(fK)), fired(nmodes), vars].
    # This matches afterGlobalEvent.m and the per-mode slicing logic.
    n_mode = mode_spaces[0][0].shape[0]
    idle_combined = mode_spaces[0][0].copy()
    phases_combined = mode_spaces[0][1].copy()
    for m in range(1, nmodes):
        idle_m, phases_m = mode_spaces[m]
        n_a = idle_combined.shape[0]
        n_b = idle_m.shape[0]
        rep_a = np.repeat(np.arange(n_a), n_b)
        rep_b = np.tile(np.arange(n_b), n_a)
        idle_combined = np.hstack([idle_combined[rep_a], idle_m[rep_b]])
        phases_combined = np.hstack([phases_combined[rep_a], phases_m[rep_b]])

    fired_block = np.zeros((idle_combined.shape[0], nmodes), dtype=float)
    trans_space = np.hstack([idle_combined, phases_combined, fired_block])

    state_var = _space_local_vars(sn, ind)
    if state_var is not None and state_var.size > 0:
        trans_space = _cartesian_2d(trans_space, state_var)

    return trans_space


def _from_marginal_bounds(sn, ind, ub, cap, options=None):
    """
    Generate all states within marginal bounds.

    Port from MATLAB fromMarginalBounds.m.
    Iterates over all marginals n with 0 <= n[r] <= ub[r] for each class r,
    calls fromMarginal for each, and filters by capacity.

    Args:
        sn: NetworkStruct
        ind: Node index (0-based)
        ub: Upper bounds per class (array of length R)
        cap: Total capacity at this node
        options: Solver options

    Returns:
        State space matrix (n_states x n_cols)
    """
    R = sn.nclasses
    ub = np.atleast_1d(ub).astype(int)

    # Generate all marginals within bounds
    ranges = [range(int(ub[r]) + 1) for r in range(R)]
    space_list = []
    max_cols = 0

    for n_tuple in iproduct(*ranges):
        n = np.array(n_tuple, dtype=float)
        state = fromMarginal(sn, ind, n, options)
        if isinstance(state, np.ndarray) and state.size > 0:
            state = np.atleast_2d(state)
            if state.shape[1] > max_cols:
                max_cols = state.shape[1]
            space_list.append(state)

    if not space_list:
        return np.zeros((1, max(R, 1)))

    # Pad to same width (left-pad with zeros, matching MATLAB behavior)
    padded = []
    for s in space_list:
        if s.shape[1] < max_cols:
            pad = np.zeros((s.shape[0], max_cols - s.shape[1]))
            s = np.hstack([pad, s])
        padded.append(s)

    space = np.vstack(padded)
    space = np.unique(space, axis=0)

    # Filter by capacity constraints
    if sn.isstateful[ind] and space.shape[0] > 0:
        keep = []
        ist = int(sn.nodeToStation[ind]) if sn.isstation[ind] else -1
        # toMarginal expects state_i to include trailing nvars columns and
        # strips them off; we built buffer-only rows above, so pad with zeros
        # for non-station stateful nodes (e.g. Cache, Router) to match the
        # expected layout — otherwise the buf_end slice would drop the last
        # nvars buffer columns and let bogus oversize states past the cap.
        nvars_sum = int(np.sum(sn.nvars[ind])) if sn.nvars is not None else 0
        pad_needed = (not sn.isstation[ind]) and nvars_sum > 0
        for s_idx in range(space.shape[0]):
            try:
                row = space[s_idx:s_idx + 1]
                if pad_needed:
                    row = np.hstack([row, np.zeros((1, nvars_sum))])
                ni_nir = toMarginal(sn, ind, row)
                ni = ni_nir[0]
                nir = ni_nir[1]
                ni_val = float(ni.ravel()[0]) if isinstance(ni, np.ndarray) else float(ni)
                nir_vec = nir.ravel() if isinstance(nir, np.ndarray) else np.atleast_1d(nir)

                if sn.isstation[ind]:
                    if np.all(nir_vec[:R] <= sn.classcap[ist, :R]) and ni_val <= cap:
                        keep.append(s_idx)
                else:
                    if ni_val <= cap:
                        keep.append(s_idx)
            except Exception:
                keep.append(s_idx)  # Keep if toMarginal fails

        if keep:
            space = space[keep]
        else:
            return np.zeros((1, max_cols))

    # Reverse order so states with jobs in phase 1 come first
    space = space[::-1]
    return space


def _space_local_vars(sn, ind, first_only=False):
    """
    Generate state space for local state variables at a node.

    Port from MATLAB spaceLocalVars.m.

    Args:
        sn: NetworkStruct
        ind: Node index (0-based)
        first_only: when True, return only the first local-var row without
            enumerating the full cache space (used by SSA to seed one valid
            initial state). The cache component is built directly via
            _first_cache_state; the (small) round-robin pointer component is
            still enumerated and the leading row taken.

    Returns:
        State space matrix for local variables, or None if no local vars
    """
    R = sn.nclasses
    space = None

    # Cache local variables (item positions in cache lists)
    nt_val = _get_nodetype_val(sn, ind)
    cache_val = int(NodeType.CACHE.value) if hasattr(NodeType.CACHE, 'value') else int(NodeType.CACHE)

    if nt_val == cache_val:
        nparam = sn.nodeparam[ind] if sn.nodeparam is not None and ind in sn.nodeparam else None
        if nparam is not None:
            n_items = nparam.get('nitems', 0) if isinstance(nparam, dict) else getattr(nparam, 'nitems', 0)
            m = nparam.get('itemcap', []) if isinstance(nparam, dict) else getattr(nparam, 'itemcap', [])
            rsc = nparam.get('retrieval_system_capacity', 0) if isinstance(nparam, dict) else getattr(nparam, 'retrieval_system_capacity', 0)
            m = np.atleast_1d(m).astype(int)
            if n_items > 0 and len(m) > 0:
                space = (_first_cache_state(n_items, m, rsc) if first_only
                         else _space_cache(n_items, m, rsc))

    # Round-robin routing variables
    if hasattr(sn, 'routing') and sn.routing is not None:
        for r in range(R):
            routing_val = None
            if isinstance(sn.routing, dict):
                routing_val = sn.routing.get((ind, r), None)
            else:
                # ndarray or MatrixArray — both support [ind, r] indexing.
                try:
                    routing_val = sn.routing[ind, r]
                except (TypeError, KeyError, IndexError):
                    routing_val = None

            if routing_val is not None:
                rv = int(routing_val.value) if hasattr(routing_val, 'value') else int(routing_val)
                rr_val = int(RoutingStrategy.RROBIN.value) if hasattr(RoutingStrategy.RROBIN, 'value') else int(RoutingStrategy.RROBIN)
                wrr_val = int(RoutingStrategy.WRROBIN.value) if hasattr(RoutingStrategy.WRROBIN, 'value') else int(RoutingStrategy.WRROBIN)

                if rv == wrr_val:
                    # Weighted round-robin: the pointer slot holds a POSITION
                    # (1..len) in the weighted-outlink cycle, not a node index,
                    # because a destination may repeat. Enumerate positions.
                    wol = _wrr_weighted_outlinks(sn, ind, r)
                    if wol is not None and len(wol) > 0:
                        rr_space = (np.arange(len(wol)) + 1).reshape(-1, 1)
                        space = _cartesian_2d(space, rr_space) if space is not None else rr_space
                        continue
                if rv in (rr_val, wrr_val):
                    # Get outlinks for this class from nodeparam if present.
                    nparam = sn.nodeparam[ind] if sn.nodeparam is not None and ind in sn.nodeparam else None
                    outlinks_all = None
                    if nparam is not None:
                        if isinstance(nparam, dict):
                            outlinks_all = nparam.get('outlinks', None)
                        elif isinstance(nparam, (list, tuple)):
                            outlinks_all = nparam[r] if r < len(nparam) else None
                            if isinstance(outlinks_all, dict):
                                outlinks_all = outlinks_all.get('outlinks', None)
                            else:
                                outlinks_all = getattr(outlinks_all, 'outlinks', None)
                        else:
                            outlinks_all = getattr(nparam, 'outlinks', None)

                    ol = None
                    if outlinks_all is not None:
                        if isinstance(outlinks_all, (list, tuple)):
                            ol = np.atleast_1d(outlinks_all[r]) if r < len(outlinks_all) else np.array([])
                        else:
                            ol = np.atleast_1d(outlinks_all)
                    if ol is None or ol.size == 0:
                        # nodeparam.outlinks is not populated on the native path;
                        # derive the round-robin destinations from the connection
                        # matrix (matches flat _get_rrobin_outlinks).
                        if hasattr(sn, 'connmatrix') and sn.connmatrix is not None:
                            cm = np.asarray(sn.connmatrix)
                            if ind < cm.shape[0]:
                                ol = np.where(cm[ind, :] > 0)[0]

                    if ol is not None and ol.size > 0:
                        # pointer state holds the destination node index of the
                        # next outlink to be used (matches _update_rrobin_pointer,
                        # which advances current_val through the outlink list).
                        rr_space = np.asarray(ol).reshape(-1, 1)
                        if space is not None:
                            space = _cartesian_2d(space, rr_space)
                        else:
                            space = rr_space

    if first_only and space is not None:
        return space[:1]
    return space


def _space_cache(n, m, retrieval_system_capacity=0):
    """
    Generate all cache state vectors.

    Port from MATLAB spaceCache.m.
    Items are 1-indexed. Each position in the cache holds one item.
    The cache has h lists with sizes m[0], m[1], ..., m[h-1].

    When ``retrieval_system_capacity`` > 0 each emitted row has
    totalCacheCapacity + retrieval_system_capacity columns; the leading
    totalCacheCapacity columns are the cache contents and the trailing
    columns are retrieval-system slots (0 = empty, else a 1-based item id).

    Args:
        n: Total number of items
        m: Array of list capacities
        retrieval_system_capacity: extra retrieval-system slots (default 0)

    Returns:
        Matrix where each row is a cache state (item IDs in positions)
    """
    from itertools import combinations, permutations

    m = np.atleast_1d(m).astype(int)
    total_slots = int(np.sum(m))
    rsc = int(retrieval_system_capacity)
    # The retrieval system is encoded as a per-item occupancy bitmap (one column
    # per item) appended after the cache contents; bit i is set iff item i+1 is
    # being retrieved. When there is no retrieval system the bitmap is omitted.
    retrieval_width = n if rsc > 0 else 0
    n_vars = total_slots + retrieval_width

    if total_slots == 0 or n == 0:
        return np.zeros((1, max(n_vars, 1)))

    items = list(range(1, n + 1))
    states = []

    # Cache contents: every ordered placement of total_slots distinct items.
    for cache_combo in combinations(items, total_slots):
        cached = set(cache_combo)
        remaining = [it for it in items if it not in cached]
        cache_perms = list(permutations(cache_combo))
        # Retrieval-system occupancy: every subset of the remaining items of size
        # <= rsc, encoded as a one-hot bitmap (membership is a set, so each
        # physical configuration maps to exactly one bitmap row).
        for s in range(0, rsc + 1):
            for retr_combo in combinations(remaining, s):
                bitmap = [0] * retrieval_width
                for it in retr_combo:
                    bitmap[it - 1] = 1
                for perm in cache_perms:
                    states.append(list(perm) + bitmap)

    if not states:
        return np.zeros((1, n_vars))

    return np.array(states, dtype=float)


def _first_cache_state(n, m, retrieval_system_capacity=0):
    """First cache state row, built directly without enumerating the space.

    Equivalent to ``_space_cache(n, m, rsc)[0:1]`` but O(total_slots) instead of
    O(perm(n, total_slots)): the first enumerated configuration places items
    1..total_slots in the cache slots in order with an empty retrieval bitmap.
    Used to seed the SSA initial state -- a simulation needs one valid
    configuration, not the whole space (mirrors the single initial cache state
    MATLAB keeps in ``sn.state`` and uses via solver_ssa_analyzer.m's
    ``init_state = sn.state``). For large caches (e.g. 1000 items / 50 slots)
    this avoids the MemoryError that materializing the full space would cause.
    """
    m = np.atleast_1d(m).astype(int)
    total_slots = int(np.sum(m))
    rsc = int(retrieval_system_capacity)
    retrieval_width = n if rsc > 0 else 0
    n_vars = total_slots + retrieval_width
    if total_slots == 0 or n == 0 or total_slots > n:
        return np.zeros((1, max(n_vars, 1)))
    row = list(range(1, total_slots + 1)) + [0] * retrieval_width
    return np.array([row], dtype=float)


def _raise_ctmc_too_large(nstates, cap):
    """Abort CTMC enumeration when the reachable state count exceeds the cap.

    Raised as MemoryError so a solver-fallback chain (MVA->MAM->CTMC->LDES)
    treats it like the OOM it pre-empts and degrades to a simulation/approximate
    engine instead of the process being OOM-killed by the kernel.
    """
    raise MemoryError(
        "SolverCTMC aborted: reachable state count exceeded %d (cap %s). The "
        "model is too large for exact CTMC enumeration (often an unbounded open "
        "class or a finite-capacity region with a large token bound). Use a "
        "simulation/approximate solver (SSA, LDES, FLD), tighten the cutoff, or "
        "raise options.ctmc_max_states to override."
        % (int(nstates) - 1, ('%d' % int(cap)) if np.isfinite(cap) else 'inf'))


def _build_global_states(sn, Np, is_open_class, capacityc, options=None, max_states=float('inf')):
    """
    Build global state space from per-node state spaces.

    Port from MATLAB spaceGenerator.m lines 48-150.

    Args:
        sn: NetworkStruct with sn.space populated
        Np: Effective population per class (cutoff for open, actual for closed)
        is_open_class: Boolean array, True for open classes
        capacityc: Per-node capacity matrix (nnodes x R)
        options: Solver options; options['timeout'] (seconds) bounds the
            enumeration wall-clock cooperatively, matching MATLAB.

    Returns:
        Tuple of (state_space, state_space_hashed)
    """
    R = sn.nclasses
    nstateful = sn.nstateful
    is_closed_class = ~is_open_class

    # Count non-Source stateful nodes
    source_val = int(NodeType.SOURCE.value) if hasattr(NodeType.SOURCE, 'value') else int(NodeType.SOURCE)
    cache_val = int(NodeType.CACHE.value) if hasattr(NodeType.CACHE, 'value') else int(NodeType.CACHE)
    n_sources = 0
    for ind in range(sn.nnodes):
        nt_val = _get_nodetype_val(sn, ind)
        if nt_val == source_val:
            n_sources += 1
    nstateful_p = nstateful - n_sources

    # SPN models: use cartesian builder. The chain-station positioning logic
    # below assumes per-station per-class job counts can be read off the
    # per-node marginal — but Transition rows hold per-mode idle/phase/fired
    # counters that the marginal extraction does not distinguish from job
    # counts. Cartesian + population-aware filter handles SPN correctly.
    transition_val = int(NodeType.TRANSITION.value) if hasattr(NodeType.TRANSITION, 'value') else int(NodeType.TRANSITION)
    has_transition = False
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        for nt in sn.nodetype:
            ntv = int(nt.value) if hasattr(nt, 'value') else int(nt)
            if ntv == transition_val:
                has_transition = True
                break
    if has_transition:
        return _build_global_states_cartesian(sn, Np, is_open_class, max_states=max_states)

    # Generate chain-station positions
    chain_station_pos = _generate_chain_station_positions(sn, Np, is_open_class, is_closed_class, nstateful_p)

    if chain_station_pos is None or len(chain_station_pos) == 0:
        # Fallback: simple cartesian product approach
        return _build_global_states_cartesian(sn, Np, is_open_class, max_states=max_states)

    # For each chain-station position, find compatible per-node state hashes
    netstates = {}  # netstates[j][isf] = list of hash indices
    for j in range(chain_station_pos.shape[0]):
        netstates[j] = {}
        for ind in range(sn.nnodes):
            if not sn.isstateful[ind]:
                continue
            isf = int(sn.nodeToStateful[ind])
            nt_val = _get_nodetype_val(sn, ind)

            if nt_val == source_val:
                # Source: fixed state
                if isf in sn.space and sn.space[isf] is not None:
                    state_i = fromMarginal(sn, ind, [])
                    hashes = _get_hash_list(sn, ind, state_i)
                    netstates[j][isf] = hashes
                else:
                    netstates[j][isf] = [0]

            elif sn.isstation[ind]:
                # Station: extract marginal from chain_station_pos
                marg = _extract_marginal(sn, ind, j, chain_station_pos, nstateful_p, n_sources)
                if marg is not None and np.any(marg > capacityc[ind, :R]):
                    netstates[j][isf] = []
                else:
                    nvar_ind = int(np.sum(sn.nvars[ind])) if sn.nvars is not None else 0
                    if nvar_ind > 0 and isf in sn.space and sn.space[isf] is not None:
                        # Station carrying local variables (round-robin pointer):
                        # fromMarginal does not know about the pointer columns, so
                        # match per-node rows by their job columns and keep every
                        # pointer value for that marginal.
                        space_isf = np.atleast_2d(sn.space[isf])
                        if marg is not None:
                            # Phase- and var-aware per-class job counts (handles
                            # multi-phase service and the trailing pointer columns).
                            _, nir_all, _, _ = toMarginal(sn, ind, space_isf)
                            nir_all = np.atleast_2d(nir_all)
                            matching = [row_idx for row_idx in range(space_isf.shape[0])
                                        if np.allclose(nir_all[row_idx, :len(marg)], marg)]
                        else:
                            matching = list(range(space_isf.shape[0]))
                        netstates[j][isf] = matching if matching else []
                    else:
                        if marg is not None:
                            state_i = fromMarginal(sn, ind, marg)
                        else:
                            state_i = np.array([])
                        hashes = _get_hash_list(sn, ind, state_i)
                        netstates[j][isf] = hashes

            elif sn.isstateful[ind]:
                # Non-station stateful: match marginal from space
                marg = _extract_marginal(sn, ind, j, chain_station_pos, nstateful_p, n_sources)
                if marg is not None and np.any(marg > capacityc[ind, :R]):
                    netstates[j][isf] = []
                else:
                    if isf in sn.space and sn.space[isf] is not None:
                        space_isf = np.atleast_2d(sn.space[isf])
                        if marg is not None:
                            # Filter: keep rows where first R cols match marginal
                            matching = []
                            for row_idx in range(space_isf.shape[0]):
                                row_marg = space_isf[row_idx, :len(marg)]
                                if np.allclose(row_marg, marg):
                                    matching.append(row_idx)
                        else:
                            matching = list(range(space_isf.shape[0]))
                        # Cross-node (cache <-> queue) validity for retrieval-system
                        # caches: a retrieval job at a queue (read off the chain
                        # marginal) must be reflected in the cache occupancy bitmap
                        # and cannot also sit at the cache server. The local cache
                        # filter cannot see the queues, so without this the state
                        # space admits inconsistent (cache,queue) pairs and the
                        # returning-retrieval READ loops. Basic caches (no retrieval
                        # system) are returned unchanged. Mirrors MATLAB
                        # State.spaceGenerator.m.
                        if _get_nodetype_val(sn, ind) == cache_val and marg is not None and matching:
                            matching = _filter_cache_queue_consistency(
                                sn, ind, j, marg, matching, space_isf,
                                chain_station_pos, nstateful_p)
                        netstates[j][isf] = matching if matching else []
                    else:
                        netstates[j][isf] = [0]

    # Combine: cartesian product of per-node state hash sets.
    #
    # For open/mixed models the cutoff-based candidate enumeration can reach
    # several million combinations, so the inner loop must do no per-combo
    # numpy work. Precompute, once per node, the 2D local state space and its
    # rows as plain Python lists; the loop then only indexes and extends lists.
    # Duplicate global states are collapsed on the fly through a seen-set keyed
    # on the full state row, keeping first-occurrence order (equivalent to the
    # former np.unique + sort(unique_idx)) without materializing every
    # duplicate row and then sorting millions of them.
    space_rows = [None] * nstateful
    for isf in range(nstateful):
        if isf in sn.space and sn.space[isf] is not None:
            space_rows[isf] = [r.tolist() for r in np.atleast_2d(sn.space[isf])]

    # The enumeration is exhaustive: truncating the cartesian product would
    # silently change the stationary distribution (observed as RespT off by
    # 10x against MATLAB/JAR on multiserver phase-type stations). Runaway
    # enumeration is bounded cooperatively by the wall-clock budget below
    # (options['timeout']), matching the MATLAB spaceGenerator checkpoints.
    import time as _time
    _timeout = None
    if options is not None:
        _t = options.get('timeout', None) if isinstance(options, dict) \
            else getattr(options, 'timeout', None)
        if _t is not None and np.isfinite(_t) and _t > 0:
            _timeout = _time.monotonic() + float(_t)
    _tochk = 0

    SS_rows = []
    SSh_rows = []
    seen = set()

    for j in range(chain_station_pos.shape[0]):
        # Build list of hash arrays for each stateful node
        hash_lists = []
        position_feasible = True
        for isf in range(nstateful):
            if isf in netstates[j]:
                if netstates[j][isf]:
                    hash_lists.append(netstates[j][isf])
                else:
                    # Empty per-node set means this position's marginal is
                    # infeasible for node isf (exceeds capacity / no matching
                    # local state). MATLAB's spaceGenerator drops the whole
                    # position in that case (its lattice loop has no iterations);
                    # do NOT fall back to the first local state, which would
                    # fabricate an unreachable, population-violating global state.
                    position_feasible = False
                    break
            elif isf in sn.space and sn.space[isf] is not None:
                hash_lists.append([0])
            else:
                hash_lists.append([0])

        if not position_feasible:
            continue

        # Exhaustive cartesian product with a cooperative wall-clock checkpoint
        combo_source = iproduct(*hash_lists)
        for combo in combo_source:
            _tochk += 1
            if _tochk >= 4096:
                _tochk = 0
                if _timeout is not None and _time.monotonic() > _timeout:
                    raise RuntimeError(
                        'State space generation exceeded the wall-clock time '
                        'budget (options.timeout).')
            skip = False
            u_row = []
            for isf_idx in range(nstateful):
                h_idx = combo[isf_idx]
                if h_idx < 0:
                    skip = True
                    break
                rows = space_rows[isf_idx]
                if rows is not None:
                    if h_idx < len(rows):
                        u_row.extend(rows[h_idx])
                    else:
                        skip = True
                        break

            if not skip:
                key = tuple(u_row)
                if key not in seen:
                    seen.add(key)
                    SS_rows.append(u_row)
                    SSh_rows.append(list(combo))
                    if len(SS_rows) > max_states:
                        _raise_ctmc_too_large(len(SS_rows), max_states)

    if not SS_rows:
        return np.zeros((0, 0)), np.zeros((0, nstateful), dtype=int)

    SS = np.array(SS_rows, dtype=float)
    SSh = np.array(SSh_rows, dtype=int)

    return SS, SSh


def _filter_cache_queue_consistency(sn, ind, j, marg, matching, space_isf,
                                    chain_station_pos, nstateful_p):
    """
    Cross-node (cache <-> queue) validity for retrieval-system caches.

    Port of the cache global-validity block in MATLAB State.spaceGenerator.m.
    For chain position j: an item whose retrieval class is present at a
    retrieval queue (read off the chain marginal) must be recorded in the
    cache occupancy bitmap and must NOT also have a retrieval job at the cache
    server; an item in the bitmap but not at a queue must have its retrieval
    job ready at the cache server. ``matching`` is the list of candidate row
    indices into ``space_isf`` (rows are [per-class server | contents |
    bitmap]); returns the filtered list ([] if the marginal itself is invalid).
    Non-retrieval caches are returned unchanged.
    """
    nparam = sn.nodeparam[ind] if (sn.nodeparam is not None and ind in sn.nodeparam) else None
    if nparam is None:
        return matching

    def _g(k, d=None):
        return (nparam.get(k, d) if isinstance(nparam, dict) else getattr(nparam, k, d))

    rsc = int(_g('retrieval_system_capacity', 0) or 0)
    rc_mat = _g('retrieval_classes', None)
    rsqi = _g('retrieval_system_queue_indices', None)
    tcc = int(_g('total_cache_capacity', 0) or 0)
    if rsc <= 0 or rc_mat is None or rsqi is None:
        return matching

    rc_mat = np.atleast_2d(np.asarray(rc_mat, dtype=int))
    server_marg = np.atleast_1d(np.asarray(marg)).astype(int)
    lvs = len(server_marg)  # per-class server presence width

    source_val = int(NodeType.SOURCE.value) if hasattr(NodeType.SOURCE, 'value') else int(NodeType.SOURCE)

    # --- items currently at a retrieval queue (from the chain marginal) ---
    items_in_queue = set()
    valid_marginal = True
    items = rsqi.items() if hasattr(rsqi, 'items') else []
    for arrival_class, q_nodes in items:
        class_col = int(arrival_class)  # 0-based arrival class
        if class_col >= rc_mat.shape[1]:
            continue
        q_iter = q_nodes if hasattr(q_nodes, '__iter__') else [q_nodes]
        for q_idx in q_iter:
            q_idx = int(q_idx)
            n_prec = sum(1 for i in range(q_idx) if _get_nodetype_val(sn, i) == source_val)
            q_off = int(sn.nodeToStateful[q_idx]) - n_prec
            for item in range(rc_mat.shape[0]):
                r_class = int(rc_mat[item, class_col])
                if r_class < 0:
                    continue
                col = r_class * nstateful_p + q_off
                if col >= chain_station_pos.shape[1] or chain_station_pos[j, col] == 0:
                    continue
                # item is at a queue: no retrieval job may also be at the cache
                # server, and the item can be at only one queue
                if (r_class < lvs and server_marg[r_class] > 0) or (item in items_in_queue):
                    valid_marginal = False
                    break
                items_in_queue.add(item)
            if not valid_marginal:
                break
        if not valid_marginal:
            break
    if not valid_marginal:
        return []

    # --- per-state bitmap consistency ---
    keep = []
    n_cols = space_isf.shape[1]
    for ri in matching:
        st = space_isf[ri]
        valid = True
        items_in_rs = set()
        for col in range(lvs + tcc, n_cols):
            if st[col] == 0:
                continue
            item = col - (lvs + tcc)
            for jac in range(rc_mat.shape[1]):
                r_class = int(rc_mat[item, jac])
                if r_class < 0:
                    continue
                # item in the bitmap but not at a queue requires a retrieval job
                # ready to depart/arrive at the cache server
                if item not in items_in_queue and (r_class >= n_cols or st[r_class] == 0):
                    valid = False
                    break
            if not valid:
                break
            items_in_rs.add(item)
        if valid:
            for it in items_in_queue:
                if it not in items_in_rs:
                    valid = False
                    break
        if valid:
            keep.append(ri)
    return keep


def _build_global_states_cartesian(sn, Np, is_open_class, max_states=float('inf')):
    """
    Build global state space using chain-station decomposition.

    Per-node states are split into:
      * carrying nodes — Places/Queues whose marginal contributes to per-class
        job counts. We enumerate only the (per-node, per-class) marginals that
        sum to Np exactly for closed classes (and to ≤ Np for open classes).
      * free nodes — Transitions, Sources, Sinks, Caches, Routers whose state
        is decoupled from the population balance. We Cartesian-product their
        state indices.

    The result is the product of (carrying combinations passing the
    population filter) × (free Cartesian). This avoids the O(prod n_states)
    blow-up of a naive full Cartesian + filter.

    Args:
        sn: NetworkStruct with sn.space populated
        Np: Population limits per class (cutoff for open, exact pop for closed)
        is_open_class: Boolean array

    Returns:
        Tuple of (state_space, state_space_hashed)
    """
    R = sn.nclasses
    nstateful = sn.nstateful
    is_closed = ~is_open_class

    transition_v = int(NodeType.TRANSITION.value) if hasattr(NodeType.TRANSITION, 'value') else int(NodeType.TRANSITION)
    source_v = int(NodeType.SOURCE.value) if hasattr(NodeType.SOURCE, 'value') else int(NodeType.SOURCE)

    def _is_free(isf):
        if not hasattr(sn, 'statefulToNode') or sn.statefulToNode is None:
            return False
        node_ind = int(sn.statefulToNode[isf])
        if node_ind < 0 or not hasattr(sn, 'nodetype') or sn.nodetype is None:
            return False
        nt = sn.nodetype[node_ind]
        nt_int = int(nt.value) if hasattr(nt, 'value') else int(nt)
        return nt_int == transition_v or nt_int == source_v

    # Pre-compute per-node state spaces, marginals, and the carrying/free
    # partition.
    space_per_isf = {}
    n_states_per_isf = {}
    for isf in range(nstateful):
        if isf in sn.space and sn.space[isf] is not None:
            space_per_isf[isf] = np.atleast_2d(sn.space[isf])
            n_states_per_isf[isf] = space_per_isf[isf].shape[0]
        else:
            space_per_isf[isf] = None
            n_states_per_isf[isf] = 1

    carrying_isfs = []
    free_isfs = []
    for isf in range(nstateful):
        if _is_free(isf) or space_per_isf[isf] is None:
            free_isfs.append(isf)
        else:
            carrying_isfs.append(isf)

    # For each carrying node, group hashes by their per-class marginal.
    carrying_marginals = {}
    for isf in carrying_isfs:
        space_isf = space_per_isf[isf]
        node_ind = int(sn.statefulToNode[isf])
        marg_to_hashes = {}
        for h in range(space_isf.shape[0]):
            row = space_isf[h]
            try:
                result = toMarginal(sn, node_ind, row.reshape(1, -1))
                nir = result[1].ravel()[:R]
            except Exception:
                nir = np.zeros(R)
            key = tuple(int(round(float(v))) for v in nir)
            marg_to_hashes.setdefault(key, []).append(h)
        carrying_marginals[isf] = marg_to_hashes

    # Enumerate population-respecting marginal assignments across carrying
    # nodes via DFS. ``remaining`` tracks the budget per class.
    target = [int(round(float(Np[r]))) for r in range(R)]

    valid_assignments = []  # each: list of (isf, marg_tuple) in carrying order

    def recurse(idx, remaining, assignment):
        if idx == len(carrying_isfs):
            for r in range(R):
                if is_closed[r] and remaining[r] != 0:
                    return
                if (not is_closed[r]) and remaining[r] < 0:
                    return
            valid_assignments.append(list(assignment))
            return
        isf = carrying_isfs[idx]
        for marg_tuple in carrying_marginals[isf].keys():
            new_remaining = list(remaining)
            ok = True
            for r in range(R):
                new_remaining[r] -= marg_tuple[r]
                # Closed: cannot overspend; open: cannot exceed cutoff in the
                # negative direction either (already at 0 means not enabled).
                if is_closed[r] and new_remaining[r] < 0:
                    ok = False
                    break
                if (not is_closed[r]) and new_remaining[r] < 0:
                    ok = False
                    break
            if not ok:
                continue
            assignment.append((isf, marg_tuple))
            recurse(idx + 1, new_remaining, assignment)
            assignment.pop()

    recurse(0, list(target), [])

    if not valid_assignments and not carrying_isfs:
        # Pure free-node network: just take the full free Cartesian.
        valid_assignments = [[]]

    # Build (free_isfs Cartesian) once; reused per carrying assignment.
    free_index_ranges = [range(n_states_per_isf[isf]) for isf in free_isfs]

    SS_rows = []
    SSh_rows = []

    # u_row template: per-isf width offsets for fast assembly
    isf_width = {isf: (space_per_isf[isf].shape[1] if space_per_isf[isf] is not None else 0)
                 for isf in range(nstateful)}

    for assignment in valid_assignments:
        # Cartesian over hash choices for each carrying isf consistent with
        # the chosen marginal.
        carrying_hash_lists = [carrying_marginals[isf][marg] for (isf, marg) in assignment]
        for c_combo in iproduct(*carrying_hash_lists):
            carrying_hash_by_isf = {assignment[i][0]: c_combo[i] for i in range(len(assignment))}
            for f_combo in iproduct(*free_index_ranges):
                free_hash_by_isf = {free_isfs[i]: f_combo[i] for i in range(len(free_isfs))}

                # Assemble in nstateful order
                h_row = [0] * nstateful
                u_row = []
                for isf in range(nstateful):
                    if isf in carrying_hash_by_isf:
                        h_row[isf] = carrying_hash_by_isf[isf]
                    elif isf in free_hash_by_isf:
                        h_row[isf] = free_hash_by_isf[isf]
                    if space_per_isf[isf] is not None:
                        u_row.extend(space_per_isf[isf][h_row[isf]].tolist())

                SS_rows.append(u_row)
                SSh_rows.append(h_row)
                if len(SS_rows) > max_states:
                    _raise_ctmc_too_large(len(SS_rows), max_states)

    if not SS_rows:
        return np.zeros((0, 0)), np.zeros((0, nstateful), dtype=int)

    SS = np.array(SS_rows, dtype=float)
    SSh = np.array(SSh_rows, dtype=int)

    # Remove duplicates
    _, unique_idx = np.unique(SS, axis=0, return_index=True)
    unique_idx = np.sort(unique_idx)
    SS = SS[unique_idx]
    SSh = SSh[unique_idx]

    return SS, SSh


def _build_state_space_aggr(sn, state_space_hashed):
    """
    Build aggregated state space (per-station per-class job counts).

    Port from MATLAB/JAR ctmc_ssg aggregation loop.

    Args:
        sn: NetworkStruct
        state_space_hashed: Hashed state space (n_states x nstateful)

    Returns:
        state_space_aggr: Aggregated matrix (n_states x M*R)
    """
    if state_space_hashed.size == 0:
        return np.zeros((0, sn.nstations * sn.nclasses))

    M = sn.nstations
    R = sn.nclasses
    n_states = state_space_hashed.shape[0]
    aggr = np.zeros((n_states, M * R))

    for s in range(n_states):
        for ind in range(sn.nnodes):
            if not sn.isstateful[ind] or not sn.isstation[ind]:
                continue
            isf = int(sn.nodeToStateful[ind])
            ist = int(sn.nodeToStation[ind])

            if isf not in sn.space or sn.space[isf] is None:
                continue

            h_idx = int(state_space_hashed[s, isf])
            space_isf = np.atleast_2d(sn.space[isf])
            if h_idx >= space_isf.shape[0]:
                continue

            state_row = space_isf[h_idx:h_idx + 1]
            try:
                result = toMarginal(sn, ind, state_row)
                nir = result[1]
                nir_vec = nir.ravel()[:R]
                aggr[s, ist * R:(ist + 1) * R] = nir_vec
            except Exception:
                pass

    return aggr


def _generate_chain_station_positions(sn, Np, is_open, is_closed, nstateful_p):
    """
    Generate all valid job distributions across stateful nodes.

    Port from MATLAB spaceGenerator.m lines 48-82.

    Args:
        sn: NetworkStruct
        Np: Population per class
        is_open: Boolean array for open classes
        is_closed: Boolean array for closed classes
        nstateful_p: Number of stateful nodes excluding sources

    Returns:
        Matrix of chain-station positions (n_positions x nstateful_p * R)
    """
    R = sn.nclasses

    if nstateful_p == 0:
        return None

    # Generate all valid distributions of closed class populations
    positions = []
    Np_int = Np.astype(int)

    # Iterate over all possible class populations (for open classes)
    ranges = []
    for r in range(R):
        ranges.append(range(int(Np_int[r]) + 1))

    for n_tuple in iproduct(*ranges):
        n = np.array(n_tuple, dtype=int)
        # Check closed class constraint: closed classes must have exact population
        if not np.all(n[is_closed] == Np_int[is_closed]):
            continue

        # Distribute n across nstateful_p nodes
        # Use spaceClosedMultiCS equivalent
        node_dists = _distribute_across_nodes(n, nstateful_p, sn.chains if hasattr(sn, 'chains') else None)
        for dist in node_dists:
            positions.append(dist)

    if not positions:
        return None

    pos_array = np.array(positions, dtype=float)
    # Remove duplicates
    pos_array = np.unique(pos_array, axis=0)
    return pos_array


def _distribute_across_nodes(n, nstateful_p, chains):
    """
    Distribute n jobs across nstateful_p nodes respecting chain constraints.

    Simplified port of MATLAB spaceClosedMultiCS.

    Args:
        n: Total jobs per class (array of length R)
        nstateful_p: Number of non-source stateful nodes
        chains: Chain membership matrix (nchains x R)

    Returns:
        List of distribution vectors, each of length nstateful_p * R
    """
    R = len(n)
    if nstateful_p == 0:
        return []

    # Faithful port of MATLAB State.spaceClosedMultiCS: for each chain, the
    # chain's TOTAL population (summed over its classes) is first split among
    # the classes of the chain in every possible way (class switching can move
    # jobs between classes of the same chain), and only then distributed across
    # the nodes. The previous "each class independently by n[r]" approach missed
    # every class-switched state (e.g. a class with njobs==0 that jobs switch
    # into never received any jobs).
    if chains is None:
        chains_mat = np.eye(R, dtype=int)
    else:
        chains_mat = np.atleast_2d(np.asarray(chains))
    C = chains_mat.shape[0]

    chain_inchain = []   # classes belonging to each chain
    chain_splits = []    # per-chain: list of class-split vectors of its total population
    covered = np.zeros(R, dtype=bool)
    for c in range(C):
        inchain = [r for r in range(R) if c < chains_mat.shape[0] and chains_mat[c, r] != 0]
        if not inchain:
            chain_inchain.append([])
            chain_splits.append([np.array([], dtype=int)])
            continue
        for r in inchain:
            covered[r] = True
        total = int(sum(int(n[r]) for r in inchain))
        chain_inchain.append(inchain)
        chain_splits.append(_multichoose_constrained(total, len(inchain)))
    # Any class not assigned to a chain is treated as its own singleton chain.
    for r in range(R):
        if not covered[r]:
            chain_inchain.append([r])
            chain_splits.append(_multichoose_constrained(int(n[r]), 1))

    results = []
    for split_combo in iproduct(*chain_splits):
        # Reconstruct the per-class population vector subN for this class-split.
        subN = np.zeros(R, dtype=int)
        for cinchain, split in zip(chain_inchain, split_combo):
            for idx, r in enumerate(cinchain):
                subN[r] = int(split[idx])
        # spaceClosedMulti: distribute subN[r] of each class across the nodes.
        per_class_dists = [_multichoose_constrained(int(subN[r]), nstateful_p) for r in range(R)]
        for combo in iproduct(*per_class_dists):
            row = np.zeros(nstateful_p * R)
            for r in range(R):
                dist = combo[r]
                for isf in range(nstateful_p):
                    row[r * nstateful_p + isf] = dist[isf]
            results.append(row)

    return results


def _multichoose_constrained(n, k):
    """
    Generate all ways to distribute n identical items across k bins.

    Returns list of arrays, each of length k, summing to n.
    """
    if k == 0:
        return [np.array([])] if n == 0 else []
    if k == 1:
        return [np.array([n])]
    if n == 0:
        return [np.zeros(k, dtype=int)]

    results = []
    for i in range(n + 1):
        for rest in _multichoose_constrained(n - i, k - 1):
            results.append(np.concatenate([[i], rest]).astype(int))
    return results


def _extract_marginal(sn, ind, j, chain_station_pos, nstateful_p, n_sources):
    """
    Extract per-class job counts for a node from chain_station_pos.

    Port from MATLAB spaceGenerator.m line 93:
    stateMarg_i = chainStationPos(j, (isf-n_preceding_sources):nstatefulp:end)

    Args:
        sn: NetworkStruct
        ind: Node index
        j: Position index
        chain_station_pos: Position matrix
        nstateful_p: Number of non-source stateful nodes
        n_sources: Number of source nodes

    Returns:
        Marginal array (R,) or None
    """
    isf = int(sn.nodeToStateful[ind])
    R = sn.nclasses

    # Count sources before this node
    source_val = int(NodeType.SOURCE.value) if hasattr(NodeType.SOURCE, 'value') else int(NodeType.SOURCE)
    n_preceding_sources = 0
    for i in range(ind):
        nt_val = _get_nodetype_val(sn, i)
        if nt_val == source_val:
            n_preceding_sources += 1

    # MATLAB: (isf - n_preceding_sources) is 0-based offset for non-source stateful nodes
    # In the position vector, columns are [node0_class0, node1_class0, ..., node0_class1, ...]
    # We pick every nstateful_p-th element starting from (isf - n_preceding_sources)
    col_start = isf - n_preceding_sources
    if col_start < 0 or col_start >= nstateful_p:
        return None

    try:
        marg = np.zeros(R)
        for r in range(R):
            col = r * nstateful_p + col_start
            if col < chain_station_pos.shape[1]:
                marg[r] = chain_station_pos[j, col]
        return marg
    except (IndexError, ValueError):
        return None


def _space_row_index(sn, isf, raw_space, space):
    """
    Return a dict mapping each state-space row (as a rounded-tuple key) to its
    row index in sn.space[isf], memoized on sn so it is built once per node.

    The map replaces the former O(n_query * n_space) np.allclose scan in
    _get_hash_list with an O(1) lookup. State-space entries are small
    non-negative integers (job/phase counts, cache bitmaps, round-robin
    pointers), so an exact match on the 6-decimal-rounded tuple is equivalent
    to the previous allclose(atol=1e-10) test. The first occurrence wins, so a
    duplicated row resolves to the same index the linear scan would have found.

    The cache is invalidated when sn.space[isf] is replaced (identity change)
    or resized (row-count change), so a rebuilt local state space is never
    matched against a stale index.
    """
    cache = getattr(sn, '_ctmc_hashidx', None)
    if cache is None:
        cache = {}
        sn._ctmc_hashidx = cache
    token = (id(raw_space), space.shape[0])
    entry = cache.get(isf)
    if entry is not None and entry[0] == token:
        return entry[1]
    index = {}
    for j in range(space.shape[0]):
        key = tuple(np.round(space[j], 6))
        if key not in index:
            index[key] = j
    cache[isf] = (token, index)
    return index


def _get_hash_list(sn, ind, state_i):
    """
    Find hash indices for states in sn.space[isf].

    Port from MATLAB State.getHash.

    Args:
        sn: NetworkStruct
        ind: Node index
        state_i: State rows to look up

    Returns:
        List of hash (row) indices
    """
    isf = int(sn.nodeToStateful[ind])
    if isf not in sn.space or sn.space[isf] is None:
        return [0]

    if not isinstance(state_i, np.ndarray) or state_i.size == 0:
        return []

    state_i = np.atleast_2d(state_i)
    raw_space = sn.space[isf]
    space = np.atleast_2d(raw_space)
    index = _space_row_index(sn, isf, raw_space, space)

    hashes = []
    for i in range(state_i.shape[0]):
        row = state_i[i]
        # Pad to match space width
        if len(row) < space.shape[1]:
            row = np.concatenate([np.zeros(space.shape[1] - len(row)), row])
        elif len(row) > space.shape[1]:
            row = row[-space.shape[1]:]

        hashes.append(index.get(tuple(np.round(row, 6)), -1))

    return hashes


def _cartesian_2d(a, b):
    """
    Cartesian product of two 2D matrices (row-wise).

    For each pair of rows (a_row, b_row), produces [a_row, b_row].

    Args:
        a: Matrix (na x ca)
        b: Matrix (nb x cb)

    Returns:
        Matrix (na*nb x ca+cb)
    """
    a = np.atleast_2d(a)
    b = np.atleast_2d(b)
    na, ca = a.shape
    nb, cb = b.shape

    result = np.zeros((na * nb, ca + cb))
    idx = 0
    for i in range(na):
        for j in range(nb):
            result[idx, :ca] = a[i]
            result[idx, ca:] = b[j]
            idx += 1

    return result


# --- Helper functions ---

def _find_chain(sn, r):
    """Find chain index containing class r."""
    if not hasattr(sn, 'chains') or sn.chains is None:
        return None
    chains = np.atleast_2d(sn.chains)
    for c in range(chains.shape[0]):
        if chains[c, r]:
            return c
    return None


def _is_disabled_proc(sn, ist, r):
    """Check if service process at (ist, r) is disabled (NaN)."""
    if sn.proc is None:
        return False
    try:
        proc_ir = sn.proc[ist][r]
        if proc_ir is None:
            return False
        if isinstance(proc_ir, (list, tuple)) and len(proc_ir) > 0:
            D0 = np.atleast_2d(proc_ir[0])
            return bool(np.any(np.isnan(D0)))
        return False
    except (IndexError, KeyError, TypeError):
        return False


def _get_nodetype_val(sn, ind):
    """Get integer nodetype value for a node."""
    if not hasattr(sn, 'nodetype') or sn.nodetype is None or ind >= len(sn.nodetype):
        return -1
    nt = sn.nodetype[ind]
    return int(nt.value) if hasattr(nt, 'value') else int(nt)
