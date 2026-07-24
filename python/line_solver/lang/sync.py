"""
Synchronization actions for CTMC state space generation.

Each sync action pairs an active event (e.g., DEP at node i, class r)
with a passive event (e.g., ARV at node j, class s) with a routing probability.

Ported from MATLAB refreshSync.m and JAR Network.java:refreshSync().
"""

from dataclasses import dataclass
from typing import Union, Callable, Optional, List
from ..constants import EventType, RoutingStrategy, ProcessType
from .base import NodeType, SchedStrategy


@dataclass
class SyncEvent:
    """A single event in a synchronization pair."""
    event: EventType
    node: int       # 0-based node index
    job_class: int  # 0-based class index (or mode index for Transitions)
    prob: Union[float, Callable, None] = None  # routing probability or state-dep function


@dataclass
class SyncAction:
    """A pair of active + passive events that define an atomic state transition."""
    active: SyncEvent
    passive: SyncEvent


@dataclass
class ModeEvent:
    """
    A mode-keyed event for SPN global synchronization.

    Mirrors JAR ModeEvent. Used inside GlobalSync entries; weight carries
    the enabling/firing token count (PRE/POST) or 1.0 for ENABLE/FIRE actives.
    """
    event: EventType
    node: int             # 0-based node index of the participating node
    mode: int             # 0-based mode index
    weight: float = 1.0   # token weight (PRE = enabling count, POST = firing count)
    job_class: int = 0    # class index used by PRE/POST passive events


@dataclass
class GlobalSync:
    """An SPN global synchronization (one ENABLE or one FIRE per mode)."""
    active: List[ModeEvent]
    passive: List[ModeEvent]


def _make_jsq_prob(sn_ref, node_ind, dest_node):
    """Join-the-Shortest-Queue routing probability (state-dependent).

    Mirrors MATLAB refreshRoutingMatrix sub_jsq: route to dest_node iff it
    holds the minimum total job count among the linked stateful nodes; ties
    are split uniformly. With an empty local state the connection matrix
    fallback is used, as in MATLAB.
    """
    import numpy as np
    isf_i = int(sn_ref.nodeToStateful[node_ind])
    cm = np.asarray(sn_ref.connmatrix)
    linked = [int(x) for x in np.where(cm[node_ind, :] > 0)[0]]

    def prob_fn(state_before, state_after):
        from ..api.state.marginal import toMarginal
        row = state_before[isf_i] if isf_i >= 0 else None
        if row is None or np.atleast_2d(row).size == 0:
            return min(float(cm[node_ind, dest_node]), 1.0)
        n = {}
        for knd in linked:
            ksf = int(sn_ref.nodeToStateful[knd])
            if ksf < 0:
                continue
            ni = toMarginal(sn_ref, knd, np.atleast_2d(state_before[ksf])[0])[0]
            n[knd] = float(np.ravel(ni)[0])
        if int(dest_node) not in n:
            return 0.0
        nmin = min(n.values())
        if n[int(dest_node)] == nmin:
            return 1.0 / sum(1 for v in n.values() if v == nmin)
        return 0.0

    return prob_fn


def _make_sq_prob(sn_ref, node_ind, cls_r, dest_node):
    """Power-of-K-choices routing probability (state-dependent).

    Mirrors MATLAB refreshRoutingMatrix sub_sq (LDES
    selectSQDestination semantics): enumerate the m^k ordered
    candidate tuples sampled with replacement, break ties by first occurrence
    in the tuple, and return the marginal probability that dest_node is
    chosen. With memory, the previously selected destination (read from the
    node's local-variable slot) is forced as the last candidate; when the
    memory slot is absent or holds no eligible node the memory-less form is
    used, as in MATLAB.
    """
    import numpy as np
    from itertools import product
    isf_i = int(sn_ref.nodeToStateful[node_ind])
    cm = np.asarray(sn_ref.connmatrix)
    eligible = [int(x) for x in np.where(cm[node_ind, :] > 0)[0]]
    np_ir = None
    if getattr(sn_ref, 'nodeparam', None) is not None and node_ind in sn_ref.nodeparam:
        entry = sn_ref.nodeparam[node_ind]
        if isinstance(entry, dict) and isinstance(entry.get(cls_r), dict):
            np_ir = entry[cls_r]
    k_par = int(np_ir['d']) if np_ir and np_ir.get('d') else 2

    def prob_fn(state_before, state_after):
        from ..api.state.marginal import toMarginal
        row = state_before[isf_i] if isf_i >= 0 else None
        if row is None or np.atleast_2d(row).size == 0:
            return min(float(cm[node_ind, dest_node]), 1.0)
        m = len(eligible)
        if m == 0 or int(dest_node) not in eligible:
            return 0.0
        k_eff = max(1, min(k_par, m))
        n = []
        for knd in eligible:
            ksf = int(sn_ref.nodeToStateful[knd])
            ni = toMarginal(sn_ref, knd, np.atleast_2d(state_before[ksf])[0])[0]
            n.append(float(np.ravel(ni)[0]))
        pvec = np.zeros(m)
        denom = float(m) ** k_eff
        for tup in product(range(m), repeat=k_eff):
            best = tup[0]
            for i in tup[1:]:
                if n[i] < n[best]:
                    best = i
            pvec[best] += 1.0
        return float(pvec[eligible.index(int(dest_node))] / denom)

    return prob_fn


def _make_rl_prob(sn_ref, node_ind, cls_r, dest_node):
    """Reinforcement-learning routing probability (state-dependent).

    Mirrors MATLAB refreshRoutingMatrix sub_rl: with a tabular value function
    (stateSize == 0) or a quadratic linear value-function approximation
    (stateSize > 0), the destination minimizing the post-decision value is
    selected (ties split uniformly); outside the action space, or when the
    node carries no RL configuration, the policy degrades to JSQ.
    """
    import numpy as np
    isf_i = int(sn_ref.nodeToStateful[node_ind])
    cm = np.asarray(sn_ref.connmatrix)
    linked = [int(x) for x in np.where(cm[node_ind, :] > 0)[0]]
    np_ir = None
    if getattr(sn_ref, 'nodeparam', None) is not None and node_ind in sn_ref.nodeparam:
        entry = sn_ref.nodeparam[node_ind]
        if isinstance(entry, dict) and isinstance(entry.get(cls_r), dict):
            np_ir = entry[cls_r]
    valuefn = np_ir.get('valuefn') if np_ir else None
    nna = np_ir.get('nodesNeedAction', []) if np_ir else []
    state_size = int(np_ir.get('stateSize', 0)) if np_ir else 0
    jsq_fn = _make_jsq_prob(sn_ref, node_ind, dest_node)
    queue_val = int(NodeType.QUEUE.value) if hasattr(NodeType.QUEUE, 'value') else int(NodeType.QUEUE)
    nodetype = [int(nt.value) if hasattr(nt, 'value') else int(nt) for nt in sn_ref.nodetype]
    ind_queue = [i for i, nt in enumerate(nodetype) if nt == queue_val]

    def prob_fn(state_before, state_after):
        from ..api.state.marginal import toMarginal
        row = state_before[isf_i] if isf_i >= 0 else None
        if row is None or np.atleast_2d(row).size == 0:
            return min(float(cm[node_ind, dest_node]), 1.0)
        if valuefn is None or node_ind not in nna:
            return jsq_fn(state_before, state_after)
        x = np.zeros(len(ind_queue))
        for qi, knd in enumerate(ind_queue):
            ksf = int(sn_ref.nodeToStateful[knd])
            ni = toMarginal(sn_ref, knd, np.atleast_2d(state_before[ksf])[0])[0]
            x[qi] = float(np.ravel(ni)[0])
        vf = np.asarray(valuefn)
        v = {}
        for knd in linked:
            tmp = x.copy()
            for qi, q_nd in enumerate(ind_queue):
                if q_nd == knd:
                    tmp[qi] += 1.0
            if state_size == 0:
                # Tabular value function: 0-based indexing by post-decision
                # counts (MATLAB indexes 1-based at counts+1)
                idx = tmp.astype(int)
                if np.max(idx) + 1 <= vf.shape[0]:
                    v[knd] = float(vf[tuple(idx)])
            else:
                # Quadratic feature vector [1, x, {x_i*x_j, i<=j}] . coeff'
                feat = [1.0] + list(tmp)
                for i in range(len(tmp)):
                    for j in range(i, len(tmp)):
                        feat.append(tmp[i] * tmp[j])
                v[knd] = float(np.dot(np.asarray(feat), np.ravel(vf)))
        if state_size == 0:
            in_space = len(v) > 0 and float(np.max(x + 1)) < vf.shape[0]
        else:
            in_space = len(v) > 0 and float(np.max(x + 1)) < state_size
        if in_space:
            vmin = min(v.values())
            if int(dest_node) in v and v[int(dest_node)] == vmin:
                return 1.0 / sum(1 for w in v.values() if w == vmin)
            return 0.0
        return jsq_fn(state_before, state_after)

    return prob_fn


def refresh_global_sync(sn) -> "dict":
    """
    Build the per-mode SPN global synchronization map.

    Mirrors MATLAB refreshGlobalSync (called from MNetwork.refresh) and the
    JAR Network.refreshGlobalSync (Network.java:6365-6462). For every
    Transition mode emit two GlobalSync entries (in the same order as JAR/MATLAB):

      1) ENABLE  : active = ENABLE@(node,mode);
                   passive = LOCAL@(place,mode) for every place with a
                   non-zero enabling-count entry (the per-class weight is
                   carried on the firing-side as PRE in step 2).
      2) FIRE    : active = FIRE@(node,mode);
                   passive = PRE@(place,class) for every (place,class) with
                   a non-zero enabling weight, and POST@(place,class) for
                   every (place,class) with a non-zero firing weight.

    Returns:
        Dict[int, GlobalSync] keyed by emission order (0..G-1).
    """
    import numpy as np
    nclasses = int(sn.nclasses)
    nnodes = int(sn.nnodes)

    gsync = {}
    if not getattr(sn, 'nodeparam', None):
        return gsync

    transition_val = int(NodeType.TRANSITION.value) if hasattr(NodeType.TRANSITION, 'value') else int(NodeType.TRANSITION)

    for ind in range(nnodes):
        if not sn.isstateful[ind]:
            continue
        nt = sn.nodetype[ind]
        nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
        if nt_val != transition_val:
            continue
        if ind not in sn.nodeparam:
            continue
        nparam = sn.nodeparam[ind]
        nmodes = int(getattr(nparam, 'nmodes', 0) if not isinstance(nparam, dict)
                     else nparam.get('nmodes', 0))
        if nmodes <= 0:
            continue

        enabling_list = getattr(nparam, 'enabling', None) if not isinstance(nparam, dict) else nparam.get('enabling', None)
        firing_list = getattr(nparam, 'firing', None) if not isinstance(nparam, dict) else nparam.get('firing', None)
        inhibiting_list = getattr(nparam, 'inhibiting', None) if not isinstance(nparam, dict) else nparam.get('inhibiting', None)
        if enabling_list is None:
            continue

        # First pass: ENABLE events (one per mode with any enabling or inhibitor place).
        for m in range(nmodes):
            en = np.atleast_2d(np.asarray(enabling_list[m], dtype=float))
            # Reduce per (ep, class) -> per ep: sum over classes (>0 if any class enabled).
            ep_mask = np.any(en > 0, axis=1) if en.ndim == 2 else (en > 0)
            enabling_places = [int(ep) for ep in np.where(ep_mask)[0]]
            # see _kb/04-networkstruct.md (refreshGlobalSync / refreshSync) for rationale
            inhibiting_places = []
            if inhibiting_list is not None:
                inh = np.atleast_2d(np.asarray(inhibiting_list[m], dtype=float))
                ip_mask = np.any(np.isfinite(inh), axis=1) if inh.ndim == 2 else np.isfinite(inh)
                inhibiting_places = [int(ip) for ip in np.where(ip_mask)[0] if int(ip) not in enabling_places]
            if not enabling_places and not inhibiting_places:
                continue
            active = [ModeEvent(EventType.ENABLE, ind, m, 1.0)]
            passive = [ModeEvent(EventType.LOCAL, ep, m, 1.0) for ep in enabling_places]
            passive += [ModeEvent(EventType.LOCAL, ip, m, 1.0) for ip in inhibiting_places]
            gsync[len(gsync)] = GlobalSync(active=active, passive=passive)

        # Second pass: FIRE events (one per mode with any enabling or firing edge).
        for m in range(nmodes):
            en = np.atleast_2d(np.asarray(enabling_list[m], dtype=float))
            fr = np.atleast_2d(np.asarray(firing_list[m], dtype=float)) if firing_list is not None else np.zeros_like(en)
            active = [ModeEvent(EventType.FIRE, ind, m, 1.0)]
            passive = []
            # PRE events for every (place,class) with positive enabling count
            ep_idx, ec_idx = np.where(en > 0)
            for ep, ec in zip(ep_idx, ec_idx):
                passive.append(ModeEvent(EventType.PRE, int(ep), m, float(en[ep, ec]), int(ec)))
            # POST events for every (place,class) with positive firing count
            fp_idx, fc_idx = np.where(fr > 0)
            for fp, fc in zip(fp_idx, fc_idx):
                passive.append(ModeEvent(EventType.POST, int(fp), m, float(fr[fp, fc]), int(fc)))
            # see _kb/04-networkstruct.md (refreshGlobalSync / refreshSync) for rationale
            if inhibiting_list is not None:
                inh = np.atleast_2d(np.asarray(inhibiting_list[m], dtype=float))
                en_places = set(int(x) for x in np.where(np.any(en > 0, axis=1))[0])
                fr_places = set(int(x) for x in np.where(np.any(fr > 0, axis=1))[0])
                ip_mask = np.any(np.isfinite(inh), axis=1)
                for ip in np.where(ip_mask)[0]:
                    if int(ip) not in en_places and int(ip) not in fr_places:
                        passive.append(ModeEvent(EventType.LOCAL, int(ip), m, 1.0))
            if passive:
                gsync[len(gsync)] = GlobalSync(active=active, passive=passive)

    return gsync


def refresh_sync(sn) -> List[SyncAction]:
    """
    Build sync actions from network structure.

    Each sync action encodes one atomic transition:
    - PHASE/LOCAL: internal phase change at a station or transition
    - READ/READ: cache read at a cache node
    - DEP/ARV: departure from one node paired with arrival at another

    Port from MATLAB refreshSync.m / JAR Network.java:refreshSync().

    Args:
        sn: NetworkStruct

    Returns:
        List of SyncAction objects
    """
    import numpy as np

    nclasses = sn.nclasses
    nnodes = sn.nnodes
    local = nnodes  # sentinel for "no passive node" (0-based: nnodes)

    sync = []

    # Build routing mask
    has_state_dep = False
    if hasattr(sn, 'isstatedep') and sn.isstatedep is not None:
        isstatedep = np.asarray(sn.isstatedep)
        if isstatedep.any():
            has_state_dep = True

    if has_state_dep and hasattr(sn, 'rtfun') and sn.rtfun is not None:
        # State-dependent routing: evaluate rtfun with empty states
        empty_state = [np.zeros((1, 0)) for _ in range(nnodes)]
        rtmask = sn.rtfun(empty_state, empty_state)
    else:
        rt = np.asarray(sn.rt)
        rtmask = np.ceil(np.abs(rt))

    for ind in range(nnodes):
        for r in range(nclasses):
            # Phase-change actions for stations with multi-phase service
            if sn.isstation[ind]:
                ist = int(sn.nodeToStation[ind])
                if hasattr(sn, 'phasessz') and sn.phasessz is not None:
                    if sn.phasessz[ist, r] > 1:
                        sync.append(SyncAction(
                            active=SyncEvent(EventType.PHASE, ind, r, float('nan')),
                            passive=SyncEvent(EventType.LOCAL, local, r, 1.0)
                        ))
                # Reneging action (exponential patience): a waiting class-r job
                # abandons the queue at a memoryless per-job rate (passive LOCAL).
                if (getattr(sn, 'impatienceType', None) is not None
                        and int(sn.impatienceType[ist, r]) == int(ProcessType.EXP.value)):
                    sync.append(SyncAction(
                        active=SyncEvent(EventType.RENEGE, ind, r, float('nan')),
                        passive=SyncEvent(EventType.LOCAL, local, r, 1.0)
                    ))
                # see _kb/04-networkstruct.md (refreshGlobalSync / refreshSync) for rationale
                if (getattr(sn, 'retrialProc', None) is not None
                        and getattr(sn, 'retrialType', None) is not None
                        and sn.retrialProc[ist][r] is not None
                        and int(sn.retrialType[ist, r]) == int(ProcessType.EXP.value)):
                    sync.append(SyncAction(
                        active=SyncEvent(EventType.RETRY, ind, r, float('nan')),
                        passive=SyncEvent(EventType.LOCAL, local, r, 1.0)
                    ))
                # see _kb/04-networkstruct.md (refreshStruct.m field-population notes) for rationale
                _hasbd = getattr(sn, 'hasbreakdown', None)
                if (r == 0 and _hasbd is not None and ind < len(_hasbd)
                        and int(_hasbd[ind]) == 1):
                    sync.append(SyncAction(
                        active=SyncEvent(EventType.FAILURE, ind, r, float('nan')),
                        passive=SyncEvent(EventType.LOCAL, local, r, 1.0)
                    ))
                    sync.append(SyncAction(
                        active=SyncEvent(EventType.REPAIR, ind, r, float('nan')),
                        passive=SyncEvent(EventType.LOCAL, local, r, 1.0)
                    ))
                # see _kb/04-networkstruct.md (refreshGlobalSync / refreshSync) for rationale
                if sn.sched[ist] == SchedStrategy.POLLING:
                    from ..api.state.polling import polling_info
                    _pinfo = polling_info(sn, ind)
                    if _pinfo is not None and bool(_pinfo['has_sw'][r]):
                        sync.append(SyncAction(
                            active=SyncEvent(EventType.SWITCH, ind, r, float('nan')),
                            passive=SyncEvent(EventType.LOCAL, local, r, 1.0)
                        ))

            # Cache READ actions and Transition PHASE actions
            if sn.isstateful[ind]:
                nt = sn.nodetype[ind]
                nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)

                if nt_val == int(NodeType.CACHE.value if hasattr(NodeType.CACHE, 'value') else NodeType.CACHE):
                    # Cache node: check if class can read
                    if (hasattr(sn, 'nodeparam') and sn.nodeparam is not None
                            and ind in sn.nodeparam):
                        nparam = sn.nodeparam[ind]
                        pread = nparam.get('pread', None) if isinstance(nparam, dict) else getattr(nparam, 'pread', None)
                        if pread is not None and r < len(pread):
                            pr = pread[r]
                            if pr is not None and not (isinstance(pr, float) and np.isnan(pr)):
                                sync.append(SyncAction(
                                    active=SyncEvent(EventType.READ, ind, r, float('nan')),
                                    passive=SyncEvent(EventType.READ, local, r, 1.0)
                                ))

                elif nt_val == int(NodeType.TRANSITION.value if hasattr(NodeType.TRANSITION, 'value') else NodeType.TRANSITION):
                    # Transition node: phase changes per mode (only iterate once at r==0)
                    if r == 0:
                        if (hasattr(sn, 'nodeparam') and sn.nodeparam is not None
                                and ind in sn.nodeparam):
                            nparam = sn.nodeparam[ind]
                            nmodes = nparam.get('nmodes', 0) if isinstance(nparam, dict) else getattr(nparam, 'nmodes', 0)
                            for m in range(nmodes):
                                sync.append(SyncAction(
                                    active=SyncEvent(EventType.PHASE, ind, m, float('nan')),
                                    passive=SyncEvent(EventType.LOCAL, local, m, 1.0)
                                ))

                # see _kb/04-networkstruct.md (refreshGlobalSync / refreshSync) for rationale
                _is_input_read_class = False
                if (hasattr(sn, 'nodeparam') and sn.nodeparam is not None
                        and ind in sn.nodeparam):
                    _npar = sn.nodeparam[ind]
                    _pread = _npar.get('pread', None) if isinstance(_npar, dict) else getattr(_npar, 'pread', None)
                    _hitc = _npar.get('hitclass', None) if isinstance(_npar, dict) else getattr(_npar, 'hitclass', None)
                    _reads = (_pread is not None and r < len(_pread) and _pread[r] is not None
                              and not (np.isscalar(_pread[r]) and isinstance(_pread[r], float) and np.isnan(_pread[r])))
                    _has_hit = False
                    if _hitc is not None:
                        _hitc = np.atleast_1d(_hitc)
                        _has_hit = r < len(_hitc) and int(_hitc[r]) >= 0
                    _is_input_read_class = _reads and _has_hit
                if _is_input_read_class:
                    continue

                # see _kb/04-networkstruct.md (refreshGlobalSync / refreshSync) for rationale
                fork_val = int(NodeType.FORK.value) if hasattr(NodeType.FORK, 'value') else int(NodeType.FORK)
                if nt_val == fork_val:
                    continue

                # DEP/ARV routing actions
                isf = int(sn.nodeToStateful[ind])
                for jnd in range(nnodes):
                    if sn.isstateful[jnd]:
                        jsf = int(sn.nodeToStateful[jnd])
                        for s in range(nclasses):
                            row_idx = isf * nclasses + r
                            col_idx = jsf * nclasses + s
                            if row_idx < rtmask.shape[0] and col_idx < rtmask.shape[1]:
                                p = rtmask[row_idx, col_idx]
                                if p > 0:
                                    # Determine routing probability
                                    routing_s = sn.routing[ind, s] if hasattr(sn, 'routing') and sn.routing is not None else None
                                    rs_val = int(routing_s.value) if hasattr(routing_s, 'value') else int(routing_s) if routing_s is not None else -1

                                    rrobin_val = int(RoutingStrategy.RROBIN.value) if hasattr(RoutingStrategy.RROBIN, 'value') else int(RoutingStrategy.RROBIN)
                                    wrrobin_val = int(RoutingStrategy.WRROBIN.value) if hasattr(RoutingStrategy.WRROBIN, 'value') else int(RoutingStrategy.WRROBIN)
                                    jsq_val = int(RoutingStrategy.JSQ.value) if hasattr(RoutingStrategy.JSQ, 'value') else int(RoutingStrategy.JSQ)
                                    sq_val = int(RoutingStrategy.SQ.value) if hasattr(RoutingStrategy.SQ, 'value') else int(RoutingStrategy.SQ)
                                    rl_val = int(RoutingStrategy.RL.value) if hasattr(RoutingStrategy.RL, 'value') else int(RoutingStrategy.RL)

                                    if rs_val in (jsq_val, sq_val, rl_val):
                                        # see _kb/04-networkstruct.md (refreshGlobalSync / refreshSync) for rationale
                                        if r != s:
                                            continue
                                        if rs_val == jsq_val:
                                            prob = _make_jsq_prob(sn, ind, jnd)
                                        elif rs_val == sq_val:
                                            prob = _make_sq_prob(sn, ind, r, jnd)
                                        else:
                                            prob = _make_rl_prob(sn, ind, r, jnd)
                                    elif rs_val in (rrobin_val, wrrobin_val):
                                        # see _kb/04-networkstruct.md (State-space construction conventions) for rationale
                                        def _make_rr_prob(sn_ref, node_ind, cls_r, dest_node):
                                            isf_rr = int(sn_ref.nodeToStateful[node_ind])
                                            _rr = (rrobin_val, wrrobin_val)
                                            def _isrr(rt):
                                                v = int(rt.value) if hasattr(rt, 'value') else (int(rt) if rt is not None else -1)
                                                return v in _rr
                                            offset = sum(1 for rr in range(cls_r) if _isrr(sn_ref.routing[node_ind, rr]))
                                            V = int(np.sum(sn_ref.nvars[node_ind])) if sn_ref.nvars is not None else 0
                                            from ..api.state.ctmc_ssg import _map_phases_at
                                            n_map = sum(1 for _r in range(int(sn_ref.nclasses))
                                                        if _map_phases_at(sn_ref, node_ind, _r) > 1)
                                            from ..api.state.routing_pointer import wrr_weighted_outlinks
                                            wol = wrr_weighted_outlinks(sn_ref, node_ind, cls_r)
                                            def prob_fn(state_before, state_after):
                                                row = np.atleast_2d(state_after[isf_rr])
                                                if row.size == 0:
                                                    row = np.atleast_2d(state_before[isf_rr])
                                                col = row.shape[1] - V + n_map + offset
                                                if V == 0 or col < 0 or col >= row.shape[1]:
                                                    return 0.0
                                                slot = int(round(float(row[0, col])))
                                                if wol is not None and len(wol) > 0:
                                                    if slot < 1 or slot > len(wol):
                                                        return 0.0
                                                    return 1.0 if int(wol[slot - 1]) == dest_node else 0.0
                                                return 1.0 if slot == dest_node else 0.0
                                            return prob_fn
                                        prob = _make_rr_prob(sn, ind, r, jnd)
                                    else:
                                        # Static routing probability
                                        rt = np.asarray(sn.rt)
                                        prob = float(rt[row_idx, col_idx])

                                    sync.append(SyncAction(
                                        active=SyncEvent(EventType.DEP, ind, r, float('nan')),
                                        passive=SyncEvent(EventType.ARV, jnd, s, prob)
                                    ))

                # see _kb/04-networkstruct.md (refreshGlobalSync / refreshSync) for rationale
                if not sn.isstateful[ind]:
                    continue
                # see _kb/04-networkstruct.md (refreshGlobalSync / refreshSync) for rationale
                _is_cache_read_class = False
                if (hasattr(sn, 'nodeparam') and sn.nodeparam is not None
                        and ind in sn.nodeparam):
                    _np = sn.nodeparam[ind]
                    _pread = _np.get('pread', None) if isinstance(_np, dict) else getattr(_np, 'pread', None)
                    if _pread is not None and r < len(_pread):
                        _pr = _pread[r]
                        if _pr is not None and not (np.isscalar(_pr) and isinstance(_pr, float) and np.isnan(_pr)):
                            _is_cache_read_class = True
                if _is_cache_read_class:
                    continue
                rt_full = np.asarray(sn.rt)
                row_idx = isf * nclasses + r
                if row_idx < rt_full.shape[0]:
                    stateful_mass = 0.0
                    for jnd in range(nnodes):
                        if sn.isstateful[jnd]:
                            jsf = int(sn.nodeToStateful[jnd])
                            for s in range(nclasses):
                                col_idx = jsf * nclasses + s
                                if col_idx < rt_full.shape[1]:
                                    stateful_mass += float(rt_full[row_idx, col_idx])
                    exit_prob = 1.0 - stateful_mass
                    if exit_prob > 1e-9:
                        sync.append(SyncAction(
                            active=SyncEvent(EventType.DEP, ind, r, float('nan')),
                            passive=SyncEvent(EventType.LOCAL, local, r, exit_prob)
                        ))

    return sync
