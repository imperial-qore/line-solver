"""
Marginal state analysis for LINE networks (pure Python).

This module provides functions to extract and generate marginal distributions
from network states, essential for state probability computation.
"""

import numpy as np
from typing import Tuple, Optional, Union, List
from ...lang.base import (SchedStrategy, NodeType,
                          SCHED_PREEMPT as _SCHED_PREEMPT,
                          SCHED_BUFFER_CLASS_TAG as _SCHED_TAG,
                          SCHED_BUFFER_PER_CLASS_COUNT as _SCHED_COUNT)
from .multiset_perms import multiset_perms
from .reply_block import reply_width as _reply_width


def toMarginal(sn, ind: int, state_i: np.ndarray = None, phasesz: np.ndarray = None,
               phaseshift: np.ndarray = None, space_buf: np.ndarray = None,
               space_srv: np.ndarray = None, space_var: np.ndarray = None) -> Tuple[np.ndarray, ...]:
    """
    Extract marginal job distributions from global state for a specific node.

    Computes the marginal queue-length distributions and job counts for a
    specific node from the global network state, considering scheduling
    strategies and service phases.

    Args:
        sn: NetworkStruct or Network object
        ind: Node index (0-based)
        state_i: Global state vector or matrix (rows = states, columns = state components)
        phasesz: Vector of phase sizes for each class
        phaseshift: Phase shift parameters for state extraction
        space_buf: Buffer space configuration
        space_srv: Service space configuration
        space_var: Local variables configuration

    Returns:
        Tuple of:
        - ni: Total jobs in node (array: n_states)
        - nir: Total jobs per class (array: n_states x n_classes)
        - sir: Jobs in service per class (array: n_states x n_classes)
        - kir: Jobs in service per class per phase (array: n_states x n_classes x max_phases)
    """
    # Handle Network object input
    if hasattr(sn, 'get_struct'):
        sn = sn.get_struct()

    # Validate inputs
    if state_i is None:
        if hasattr(sn, 'state') and sn.state is not None:
            state_i = sn.state
        else:
            raise ValueError("state_i must be provided or available in sn.state")

    state_i = np.atleast_2d(state_i)
    n_states = state_i.shape[0]

    # Join stations on FJ-augmented structs: plain per-class count state,
    # no buffer/phase split and no service (jobs wait for synchronization)
    if getattr(sn, 'isfjaugmented', False) and sn.nodetype[ind] == NodeType.JOIN:
        R = int(sn.nclasses)
        nir = state_i[:, -R:]
        sir = np.zeros_like(nir)
        kir = np.zeros((state_i.shape[0], R, 1))
        ni = np.sum(nir, axis=1)
        return ni, nir, sir, kir

    # Handle non-station nodes (Transitions, Places, etc.)
    if hasattr(sn, 'isstation') and not sn.isstation[ind]:
        if hasattr(sn, 'nodetype') and sn.nodetype[ind] == NodeType.TRANSITION:
            # Stateful non-station node (e.g., Transition in SPN)
            _nparam = sn.nodeparam[ind] if hasattr(sn, 'nodeparam') else None
            if isinstance(_nparam, dict):
                R = _nparam.get('nmodes', 1)
            else:
                R = int(getattr(_nparam, 'nmodes', 1))
            sir = np.zeros((n_states, R))
            max_phases = max(phasesz) if phasesz is not None else 1
            kir = np.zeros((n_states, R, max_phases))

            if phasesz is not None and phaseshift is not None:
                for r in range(R):
                    for k in range(phasesz[r]):
                        kir[:, r, k] = space_srv[:, phaseshift[r] + k]
                        sir[:, r] += kir[:, r, k]

            nir = sir
            ni = np.sum(nir, axis=1)
            return ni, nir, sir, kir

    # Get network parameters
    R = sn.nclasses if hasattr(sn, 'nclasses') else 1

    # Get station-level information
    if hasattr(sn, 'nodeToStation'):
        ist = sn.nodeToStation[ind]
    else:
        ist = ind

    # Pass-and-swap / order-independent: the state is the ordered list of class
    # indices (1-based; no server/buffer split).
    if hasattr(sn, 'sched') and sn.sched[ist] == SchedStrategy.PAS:
        Vp = int(sum(sn.nvars[ind])) if hasattr(sn, 'nvars') and sn.nvars is not None else 0
        listcols = state_i[:, :state_i.shape[1] - Vp] if Vp > 0 else state_i
        nir = np.zeros((n_states, R))
        for r in range(R):
            nir[:, r] = np.sum(listcols == (r + 1), axis=1)
        # Jobs in service = positions with a nonzero rate increment Delta_mu>0.
        mu_fun = None
        if hasattr(sn, 'nodeparam') and sn.nodeparam is not None and ind in sn.nodeparam:
            _npj = sn.nodeparam[ind]
            mu_fun = _npj.get('svcRateFun') if isinstance(_npj, dict) else None
        sir = np.zeros((n_states, R))
        if mu_fun is not None:
            for row in range(n_states):
                c = listcols[row][listcols[row] > 0].astype(int)
                c0 = c - 1
                mu_prev = 0.0
                for p in range(len(c)):
                    mu_cur = float(mu_fun(c0[:p + 1]))
                    if mu_cur - mu_prev > 0:
                        sir[row, c[p] - 1] += 1
                    mu_prev = mu_cur
        else:
            sir = nir.copy()
        kir = np.zeros((n_states, R, 1))
        kir[:, :, 0] = sir
        ni = np.sum(nir, axis=1)
        return ni, nir, sir, kir

    # Set default phase parameters if not provided
    if phasesz is None:
        if hasattr(sn, 'phasessz') and sn.phasessz is not None:
            phasesz = sn.phasessz[ist]
        else:
            phasesz = np.ones(R, dtype=int)
    else:
        phasesz = np.atleast_1d(phasesz)

    if phaseshift is None:
        if hasattr(sn, 'phaseshift') and sn.phaseshift is not None:
            phaseshift = sn.phaseshift[ist]
        else:
            phaseshift = np.concatenate([[0], np.cumsum(phasesz[:-1])])
    else:
        phaseshift = np.atleast_1d(phaseshift)

    # Extract space information if not provided. Slice bounds are cast to int:
    # sn.nvars (and legacy phasessz) may carry float dtype from MATLAB-parity
    # refresh code, and numpy rejects float slice indices.
    if space_var is None or space_srv is None or space_buf is None:
        total_vars = int(sum(sn.nvars[ind])) if hasattr(sn, 'nvars') and sn.nvars is not None else 0
        total_phases = int(sum(phasesz))

        if space_var is None:
            space_var = state_i[:, -total_vars:] if total_vars > 0 else np.array([])

        if space_srv is None:
            srv_end = state_i.shape[1] - total_vars
            srv_start = srv_end - total_phases
            space_srv = state_i[:, srv_start:srv_end] if total_phases > 0 else np.array([])

        if space_buf is None:
            buf_end = state_i.shape[1] - total_phases - total_vars
            space_buf = state_i[:, :buf_end] if buf_end > 0 else np.array([])

    # Determine if distribution is exponential (single phase)
    isExponential = max(phasesz) == 1

    # Initialize output arrays
    if isExponential:
        sir = space_srv
        kir = space_srv
        nir = np.zeros((n_states, R))
    else:
        sir = np.zeros((n_states, R))
        max_phases = max(phasesz)
        kir = np.zeros((n_states, R, max_phases))
        nir = np.zeros((n_states, R))

        # Extract phase-specific information
        for r in range(R):
            for k in range(phasesz[r]):
                if phaseshift[r] + k < space_srv.shape[1]:
                    kir[:, r, k] = space_srv[:, phaseshift[r] + k]
                    sir[:, r] += kir[:, r, k]

    # Compute total jobs per class based on scheduling strategy
    if hasattr(sn, 'sched'):
        sched = sn.sched[ist]
    else:
        sched = SchedStrategy.FCFS

    # Apply scheduling strategy rules to compute nir from sir and space_buf
    if sched in [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.PSPRIO,
                 SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.LPS]:
        # These policies: jobs in service = jobs in station (no buffer distinction)
        nir = sir.copy()

    elif sched == SchedStrategy.EXT:
        # External node (Source): infinite jobs per class
        nir = np.full((n_states, R), np.inf)

    elif sched in _SCHED_TAG:
        # FCFS/LCFS: count jobs in service plus those in buffer
        # Buffer entries use 1-based class indices (MATLAB convention):
        # 0 = empty, 1 = class 0, 2 = class 1, etc.
        nir = sir.copy()
        if space_buf.size > 0:
            for r in range(R):
                # Count jobs of class r in buffer (1-based: class r stored as r+1)
                if space_buf.ndim == 1:
                    nir[:, r] += np.sum(space_buf == (r + 1))
                else:
                    nir[:, r] += np.sum(space_buf == (r + 1), axis=1)

    elif sched in _SCHED_COUNT or sched == SchedStrategy.POLLING:
        # These policies track buffer per class
        nir = sir.copy()
        if space_buf.size > 0:
            for r in range(R):
                if space_buf.ndim == 1:
                    if r < len(space_buf):
                        nir[:, r] += space_buf[r]
                else:
                    if r < space_buf.shape[1]:
                        nir[:, r] += space_buf[:, r]

    elif sched in _SCHED_PREEMPT:
        # Preempt-resume: count jobs in service plus preempted jobs held in the
        # buffer. The buffer stores interleaved [class, phase] pairs, so class
        # IDs live at the even columns (1-based: class r stored as r+1).
        nir = sir.copy()
        if space_buf.size > 0:
            buf2d = space_buf if space_buf.ndim == 2 else space_buf.reshape(1, -1)
            buf_classes = buf2d[:, 0::2]
            for r in range(R):
                nir[:, r] += np.sum(buf_classes == (r + 1), axis=1)

    else:
        # Default: jobs in station = jobs in service
        nir = sir.copy()

    # Handle disabled stations (set counts to 0). Places carry NaN rates
    # by construction since they hold tokens without service — skip them
    # so SPN buffers retain their token counts.
    is_place_station = False
    if hasattr(sn, 'nodetype') and sn.nodetype is not None:
        nt = sn.nodetype[ind]
        nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
        place_val = int(NodeType.PLACE.value) if hasattr(NodeType.PLACE, 'value') else int(NodeType.PLACE)
        is_place_station = (nt_val == place_val)

    # A Place is a token container: the [buffer, server] split in its state
    # encoding is an artifact (a transition FIRE relocates the surviving tokens
    # from the server slot into the buffer slot). For the INF-family scheduling
    # used by Places the branch above set nir = sir (server slot only), so fold
    # the buffer slot back in. This makes the marginal token count agree with
    # to_marginal_aggr (buf + srv), which drives the transition enabling degree;
    # without it the QLen measurement collapses to the oscillating server slot
    # while the dynamics stay correct.
    if is_place_station and sched in [SchedStrategy.INF, SchedStrategy.PS,
                                      SchedStrategy.PSPRIO, SchedStrategy.DPS,
                                      SchedStrategy.GPS, SchedStrategy.LPS]:
        if space_buf is not None and getattr(space_buf, 'size', 0) > 0:
            buf2d = space_buf if space_buf.ndim == 2 else space_buf.reshape(1, -1)
            if buf2d.shape[1] >= R:
                for r in range(R):
                    nir[:, r] = sir[:, r] + buf2d[:, r]
    if hasattr(sn, 'rates') and not is_place_station:
        for r in range(R):
            disabled = np.isnan(sn.rates[ist, r]) if hasattr(sn.rates[ist, r], '__iter__') else np.isnan(sn.rates[ist, r])
            if disabled:
                nir[:, r] = 0
                sir[:, r] = 0
                if kir.ndim >= 3:
                    for k in range(max(phasesz)):
                        if k < kir.shape[2]:
                            kir[:, r, k] = 0
                else:
                    kir[:, r] = 0

    # Total jobs in node
    ni = np.sum(nir, axis=1)

    return ni, nir, sir, kir


def roundMarginalPreservingChains(n, sn):
    """Round a fractional marginal queue-length matrix to integers with the
    largest remainder method, so that every closed chain keeps exactly its own
    population.

    Element-wise rounding does not: it can move a job between two classes of the
    same chain or lose one altogether, giving a state outside the state space or
    one with different chain populations (hence a different steady state). Open
    chains are rounded element-wise.

    Args:
        n: marginal queue lengths, stations by classes
        sn: NetworkStruct supplying chain membership and populations

    Returns:
        A rounded copy of n.
    """
    if n is None or sn is None or not hasattr(sn, 'chains') or not hasattr(sn, 'njobs'):
        return n
    n = np.atleast_2d(np.asarray(n, dtype=float)).copy()
    chains = np.atleast_2d(sn.chains)
    nchains = int(getattr(sn, 'nchains', chains.shape[0]))
    njobs = np.atleast_1d(sn.njobs).flatten()
    M = n.shape[0]
    for c in range(nchains):
        chain_classes = np.where(chains[c, :] > 0)[0]
        njobs_chain = sum(njobs[k] for k in chain_classes)
        if np.isinf(njobs_chain):
            for k in chain_classes:
                for i in range(M):
                    n[i, k] = round(n[i, k])
        else:
            vals = []
            indices = []
            for i in range(M):
                for k in chain_classes:
                    vals.append(n[i, k])
                    indices.append((i, k))
            vals = np.array(vals, dtype=float)
            floored = np.floor(vals)
            deficit = int(round(njobs_chain - np.sum(floored)))
            if deficit > 0:
                sort_idx = np.argsort(-(vals - floored))
                for d in range(min(deficit, len(sort_idx))):
                    floored[sort_idx[d]] += 1
            for j, (i, k) in enumerate(indices):
                n[i, k] = floored[j]
    return n


# Disciplines whose local state is phases only, with no buffer part.
_SHARE_SERVER_SCHED = [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.DPS,
                       SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.LPS]


def fromMarginal(sn, ind: int, n: Union[np.ndarray, list], options: dict = None,
                 top_row_only: bool = False) -> np.ndarray:
    """
    Generate state space with specific marginal job counts at a node.

    Creates all possible network states where the specified node has
    exactly n[r] jobs of class r. This is essential for computing
    state probabilities and performance metrics.

    Args:
        sn: NetworkStruct or Network object
        ind: Node index (0-based)
        n: Vector of job counts per class [n_classes]
        options: Optional configuration dictionary
        top_row_only: return only the first row of the space, formed in closed
            form where the branch admits it (see the share-server branch below)

    Returns:
        State space matrix where each row is a valid state with the specified
        marginal job counts at node ind.
    """
    # Handle Network object input
    if hasattr(sn, 'get_struct'):
        sn = sn.get_struct()

    if options is None:
        options = {}

    # Validate inputs
    n = np.atleast_1d(n)
    R = sn.nclasses if hasattr(sn, 'nclasses') else len(n)

    # Ensure n has correct length
    if len(n) < R:
        n = np.concatenate([n, np.zeros(R - len(n))])

    # FJ-augmented struct: the join state is the per-class count vector
    # of buffered jobs/siblings, deterministic given the marginals
    if getattr(sn, 'isfjaugmented', False) and sn.nodetype[ind] == NodeType.JOIN:
        return np.asarray(n, dtype=float).reshape(1, -1)

    # Synchronous call (REPLY signal): the node holds one server per job that
    # has left for its callee and is waiting for the reply. Those servers are
    # not derivable from the marginal n, so enumerate the held-server counts
    # here and build the rest of the state with the REMAINING servers -- with b
    # servers held, only S-b jobs can be in service, a configuration the plain
    # enumeration never produces. Recurse on a struct copy with the block
    # cleared, then append its columns, which trail the local-variable vector
    # (see api.state.reply_block).
    if _reply_width(sn, ind) > 0:
        return _from_marginal_reply(sn, ind, n, options)

    # Get station index
    if hasattr(sn, 'nodeToStation'):
        ist = sn.nodeToStation[ind]
    else:
        ist = ind

    # Handle non-station nodes (Transitions, cached variables, etc.)
    if hasattr(sn, 'isstation') and not sn.isstation[ind]:
        # For stateful non-station nodes, return empty state or fixed state
        if hasattr(sn, 'space') and len(sn.space) > 0:
            return np.array([[0]])
        else:
            return np.zeros((1, sum(n)))

    # Get scheduling strategy
    if hasattr(sn, 'sched'):
        sched = sn.sched[ist]
    else:
        sched = SchedStrategy.FCFS

    # Get number of servers. Coerce to int (the per-sched state generators use
    # it in range()), but keep inf for infinite-server (Delay/INF) stations.
    if hasattr(sn, 'nservers'):
        S = sn.nservers[ist]
        if np.isfinite(S):
            S = int(S)
    else:
        S = 1

    # Get phase information
    phases = np.zeros(R, dtype=int)
    if hasattr(sn, 'proc') and sn.proc is not None:
        for r in range(R):
            if hasattr(sn.proc[ist], '__getitem__') and r < len(sn.proc[ist]):
                proc_entry = sn.proc[ist][r]
                # sn.proc stores (D0, D1); the phase count is its order.
                from ..sn.proc_form import proc_n_phases
                if proc_entry is not None:
                    phases[r] = proc_n_phases(proc_entry)
                else:
                    phases[r] = 1
            else:
                phases[r] = 1
    else:
        phases = np.ones(R, dtype=int)

    # Check capacity constraints
    if hasattr(sn, 'classcap') and sched != SchedStrategy.EXT:
        if np.any(n > sn.classcap[ist]):
            return np.array([])

    # Source (EXT): the arrival renewal process is "always in progress" and its
    # current phase must be part of the CTMC state. Each arriving class places
    # exactly one token across its arrival phases (phases(r) one-hot rows);
    # classes that do not arrive here (disabled) contribute zeros. Mirrors
    # MATLAB State.fromMarginal.m EXT branch (spaceClosedSingle(phases(r),1)).
    # Without this, a phase-type arrival (Erlang/PH/MAP) collapses to a single
    # state and is incorrectly treated as exponential at its first-phase rate.
    if sched == SchedStrategy.EXT:
        state = np.ones((1, 0))
        for r in range(R):
            nph = int(phases[r])
            if (getattr(sn, 'markidx', None) is not None
                    and ist < sn.markidx.shape[0] and sn.markidx[ist, r] > 1):
                # Marked (MMAP) non-carrier class: the modulating chain lives
                # in the carrier's phase block; single always-zero column
                # (mirrors MATLAB fromMarginal / sn.phasessz).
                init_r = np.zeros((1, 1))
            elif nph <= 0 or _ext_class_disabled(sn, ist, r):
                init_r = np.zeros((1, max(nph, 0)))
            else:
                init_r = _generate_phase_distribution(1, nph)
            state = _cartesian_product(state, init_r)
        # Attach the infinite buffer marker before the server (phase) part.
        space = np.concatenate([np.full((state.shape[0], 1), np.inf), state], axis=1)
        if space.shape[0] > 0:
            space = np.unique(space, axis=0)[::-1]
        return space

    # Pass-and-swap / order-independent queue: the local state is the full
    # ordered list of class indices (1-based; oldest at column 0), left-aligned
    # and right zero-padded to the station capacity. No server/buffer split.
    if sched == SchedStrategy.PAS:
        W = sn.cap[ist]
        if not np.isfinite(W):
            raise ValueError('PAS stations require finite capacity for state-space generation.')
        W = int(W)
        total = int(np.sum(n))
        if total == 0:
            return np.zeros((1, W))
        if total > W:
            return np.zeros((0, W))
        vi = []
        for r in range(R):
            vi.extend([r + 1] * int(n[r]))
        mi = np.array(multiset_perms(vi), dtype=float)
        if mi.ndim == 1:
            mi = mi.reshape(1, -1)
        pad = W - mi.shape[1]
        if pad > 0:
            mi = np.concatenate([mi, np.zeros((mi.shape[0], pad))], axis=1)
        return mi

    # Generate state space based on scheduling strategy
    _is_retrial = (getattr(sn, 'retrialProc', None) is not None
                   and ist < len(sn.retrialProc)
                   and any(x is not None for x in sn.retrialProc[ist]))
    # Share-server branch: the space is the cartesian product of the per-class
    # phase compositions and is uniqued then reversed below, so its first row is
    # the lexicographic maximum, i.e. every job of every class in phase 1. Formed
    # directly here because the enumeration is binomial in the population and a
    # caller that only wants the seed row must not pay for it.
    if (top_row_only and not _is_retrial
            and sched in _SHARE_SERVER_SCHED):
        row = [np.zeros(max(int(phases[r]), 1)) for r in range(R)]
        for r in range(R):
            row[r][0] = float(n[r])
        return (np.concatenate(row).reshape(1, -1) if row else np.zeros((1, 0)))

    space = _generate_state_space_for_sched(sched, S, n, phases, R, is_retrial=_is_retrial)

    # Sort and unique the state space
    if len(space) > 0:
        space = np.unique(space, axis=0)
        # Reverse sort to put states with jobs in phase 1 earlier
        space = space[::-1]

    return space


def _from_marginal_reply(sn, ind: int, n, options) -> np.ndarray:
    """fromMarginal at a station holding servers for pending synchronous calls.

    Enumerates the held-server counts b (one counter per calling class) and,
    for each b, generates the ordinary local space with S-b servers, since the
    held servers are unavailable. The sub-spaces have different buffer widths
    (held servers push jobs into the buffer), and the FCFS buffer is
    RIGHT-aligned, so the narrow rows are widened on the LEFT before stacking.

    Mirrors the reply branch of MATLAB State/fromMarginal.m.
    """
    import copy as _copy
    from itertools import product as _iproduct

    ist = int(sn.nodeToStation[ind])
    S = int(sn.nservers[ist]) if np.isfinite(sn.nservers[ist]) else 0
    rclasses = [r for r in range(int(sn.nclasses)) if sn.replyblock[ind, r] > 0]

    snb = _copy.copy(sn)
    snb.replyblock = np.asarray(sn.replyblock).copy()
    snb.replyblock[ind, :] = 0
    snb.nservers = np.asarray(sn.nservers, dtype=float).copy()

    subspaces = []
    maxw = 0
    for b in _iproduct(*[range(S + 1) for _ in rclasses]):
        if sum(b) > S:
            continue
        snb.nservers[ist] = S - sum(b)
        subspace = fromMarginal(snb, ind, n, options)
        subspace = np.atleast_2d(np.asarray(subspace))
        if subspace.size == 0:
            continue
        subspaces.append((np.asarray(b, dtype=float), subspace))
        maxw = max(maxw, subspace.shape[1])

    rows = []
    for b, subspace in subspaces:
        if subspace.shape[1] < maxw:
            subspace = np.hstack([
                np.zeros((subspace.shape[0], maxw - subspace.shape[1])), subspace])
        rows.append(np.hstack([subspace, np.tile(b, (subspace.shape[0], 1))]))
    if not rows:
        return np.array([])
    space = np.vstack(rows)
    space = np.unique(space, axis=0)[::-1]
    return space


def _generate_state_space_for_sched(sched: int, S: int, n: np.ndarray, phases: np.ndarray, R: int, is_retrial: bool = False) -> np.ndarray:
    """
    Generate local state space for a specific scheduling strategy.

    Args:
        sched: SchedStrategy enum value
        S: Number of servers
        n: Job counts per class
        phases: Number of phases per class
        R: Number of classes

    Returns:
        State space matrix for this node
    """
    space = np.array([])

    if sched == SchedStrategy.EXT:
        # Source: one job per phase per class
        state = np.ones((1, sum(phases)))
        # Add infinite buffer marker
        space = np.concatenate([np.full((1, 1), np.inf), state], axis=1)

    elif sched in _SHARE_SERVER_SCHED:
        # Jobs only in service, no buffer
        # Generate all combinations of phase distributions
        space = _cartesian_space_for_phases(n, phases, R)

    elif is_retrial:
        # Retrial station: enumerate every (in-service csrv=0..min(n,S), orbit)
        # split. The server holds csrv jobs and the rest orbit in the buffer
        # (class-id, right-aligned). The idle-server states are essential: an
        # orbiting job re-enters service only through a RETRY event. Single
        # populated class (analyzer guard enforces). Mirrors MATLAB fromMarginal.
        rr = next((r for r in range(R) if n[r] > 0), -1)
        if rr < 0:
            space = np.zeros((1, int(sum(phases))))
        else:
            maxorbit = int(n[rr])
            Sist = int(S) if np.isfinite(S) else maxorbit
            rows = []
            for csrv in range(0, min(maxorbit, Sist) + 1):
                orbit = maxorbit - csrv
                buf = np.zeros(maxorbit)
                if orbit > 0:
                    buf[maxorbit - orbit:] = rr + 1
                svec = np.array([csrv if c == rr else 0 for c in range(R)])
                srv = _cartesian_space_for_phases(svec, phases, R)
                for srow in np.atleast_2d(srv):
                    rows.append(np.concatenate([buf, np.ravel(srow).astype(float)]))
            space = np.array(rows)

    elif sched in _SCHED_TAG:
        # Jobs in ordered buffer + service
        if sum(n) == 0:
            space = np.zeros((1, 1 + sum(phases)))
        else:
            # Generate buffer and phase states
            space = _cartesian_space_fcfs_lcfs(n, phases, S, R)

    elif sched in _SCHED_PREEMPT:
        # Preempt-resume (priority) policies: track an ordered buffer of the
        # preempted jobs together with the phase each was preempted at, plus
        # the service phases. The buffer stores interleaved [class, phase]
        # pairs (matching MATLAB State.fromMarginal and the buffer format the
        # afterEventStation preempt/resume handlers read/write).
        if sum(n) == 0:
            space = np.zeros((1, 2 + sum(phases)))
        else:
            space = _cartesian_space_preempt(n, phases, S, R)

    elif sched == SchedStrategy.POLLING:
        # Un-ordered per-class buffers and a single server, as in SIRO, but NOT
        # work-conserving over the station: the server may be idle while jobs
        # wait, because it is walking towards a buffer (a switchover) or is
        # parked. Both the empty-facility and the one-job-in-service
        # configurations are therefore enumerated; the controller columns
        # appended after the other local variables discard the combinations the
        # discipline cannot occupy (see api/state/polling.polling_space).
        if S != 1:
            raise ValueError('Polling stations must have a single server.')
        rows = [np.concatenate([np.asarray(n, dtype=float),
                                np.zeros(int(sum(phases)))])]  # facility empty
        for p in range(R):
            if n[p] <= 0:
                continue
            bufp = np.asarray(n, dtype=float).copy()
            bufp[p] -= 1
            srvp = np.zeros(int(sum(phases)))
            one = np.zeros(R)
            one[p] = 1
            srvp = _cartesian_space_for_phases(one, phases, R)
            for i in range(srvp.shape[0]):
                rows.append(np.concatenate([bufp, srvp[i]]))
        space = np.array(rows, dtype=float)

    elif sched in _SCHED_COUNT:
        # Unordered per-class-count buffer + service (SIRO/SEPT/LEPT)
        if sum(n) <= S:
            # All jobs in service
            space = _cartesian_space_for_phases(n, phases, R)
            # Add zero buffer
            space = np.concatenate([np.zeros((space.shape[0], R)), space], axis=1)
        else:
            # Some jobs in buffer, some in service
            space = _generate_siro_space(n, phases, S, R)

    else:
        # Default: just track service phases
        space = _cartesian_space_for_phases(n, phases, R)

    return space if isinstance(space, np.ndarray) and space.size > 0 else np.array([])


def _cartesian_space_for_phases(n: np.ndarray, phases: np.ndarray, R: int) -> np.ndarray:
    """
    Generate state space as cartesian product of phase distributions.

    For each class r, we generate states where class r jobs are distributed
    across phases. The total for class r equals n[r].
    """
    if sum(n) == 0:
        return np.zeros((1, sum(phases)))

    # Generate phase states for each class
    class_states = []
    for r in range(R):
        if phases[r] > 1:
            # Multiple phases: distribute n[r] jobs across phases
            phase_space = _generate_phase_distribution(int(n[r]), int(phases[r]))
        else:
            # Single phase
            phase_space = np.array([[int(n[r])]])
        class_states.append(phase_space)

    # Cartesian product of all class states
    return _cartesian_product(*class_states)


def _ext_class_disabled(sn, ist: int, r: int) -> bool:
    """True if class r does not arrive at this Source (disabled arrival).

    Mirrors MATLAB's `isnan(sn.proc{ist}{r}{1})` test in fromMarginal.m. A
    disabled class contributes no arrival-phase token to the Source state.

    MATLAB and Python spell "no arrival process" differently: MATLAB stores a
    NaN D0 (and a NaN rate), Python stores None (and a zero rate). Both spellings
    must read as disabled, otherwise every class that never arrives here gets a
    spurious one-hot arrival phase and the Source row stops matching MATLAB's.
    """
    # A zero phase count is MATLAB's own marker for "no process at (ist,r)".
    if getattr(sn, 'phases', None) is not None:
        try:
            if int(np.asarray(sn.phases)[ist, r]) == 0:
                return True
        except (IndexError, KeyError, TypeError, ValueError):
            pass
    # Prefer the (D0,D1) process representation when available.
    if getattr(sn, 'proc', None) is not None:
        try:
            proc_ir = sn.proc[ist][r]
            if proc_ir is None or len(proc_ir) == 0:
                return True
            if isinstance(proc_ir, (list, tuple)) and len(proc_ir) > 0:
                return bool(np.any(np.isnan(np.atleast_2d(proc_ir[0]))))
        except (IndexError, KeyError, TypeError):
            pass
    # Fall back to the per-class arrival rate (NaN => disabled).
    if getattr(sn, 'rates', None) is not None:
        try:
            return bool(np.isnan(sn.rates[ist, r]))
        except (IndexError, KeyError, TypeError):
            pass
    return False


def _generate_phase_distribution(n_jobs: int, n_phases: int) -> np.ndarray:
    """Generate all ways to distribute n_jobs across n_phases."""
    if n_jobs == 0:
        return np.zeros((1, n_phases), dtype=int)

    if n_phases == 1:
        return np.array([[n_jobs]])

    # Recursive: distribute across phases
    states = []
    for k in range(n_jobs + 1):
        # k jobs in phase 1, rest in remaining phases
        rest = _generate_phase_distribution(n_jobs - k, n_phases - 1)
        for r in rest:
            states.append(np.concatenate([[k], r]))

    return np.array(states)


def _cartesian_space_fcfs_lcfs(n: np.ndarray, phases: np.ndarray, S: int, R: int) -> np.ndarray:
    """
    Generate state space for FCFS/LCFS scheduling.

    Tracks ordered buffer and service phases. The state format matches MATLAB:
    - Buffer positions: class ID of job in each buffer position (0=empty)
    - Server phases: phase distribution per class for jobs in service

    Args:
        n: Array of job counts per class
        phases: Array of number of phases per class
        S: Number of servers
        R: Number of classes

    Returns:
        State space matrix where each row is [buf_positions..., phase_counts...]
    """
    total_jobs = int(sum(n))

    if total_jobs == 0:
        # Empty state: 1 buffer column (0) + phase columns
        return np.zeros((1, 1 + int(sum(phases))))

    # Build list of job classes with repetition
    # e.g., n=[2,1] -> vi=[1,1,2] (using 1-based class IDs)
    vi = []
    for r in range(R):
        vi.extend([r + 1] * int(n[r]))  # 1-based class IDs

    # Generate all unique permutations of job ordering
    # Uses an efficient multiset permutation algorithm
    # instead of set(permutations(vi)) which is O(n!) even with duplicates
    mi = multiset_perms(vi)
    mi = np.array(mi) if len(mi) > 0 else np.array([vi])

    if len(mi) == 0:
        return np.zeros((1, 1 + int(sum(phases))))

    # Build states for each permutation
    all_states = []

    for perm_idx in range(mi.shape[0]):
        perm = mi[perm_idx]

        # mi_buf: class of job in buffer position i (jobs not yet in service)
        # mi_srv: class of job in server (jobs being served)
        num_in_service = int(min(total_jobs, S))
        num_in_buffer = int(total_jobs - num_in_service)

        # Buffer: first (total_jobs - S) positions
        # Server: last S positions
        if num_in_buffer > 0:
            mi_buf = perm[:num_in_buffer]
        else:
            mi_buf = np.array([0])  # Empty buffer marker

        mi_srv = perm[max(0, num_in_buffer):]

        # Count jobs of each class in service: si[r] = number of class (r+1) jobs in server
        si = np.zeros(R, dtype=int)
        for class_id in mi_srv:
            si[int(class_id) - 1] += 1  # Convert to 0-based index

        # Generate phase distributions for jobs in service
        # kstate = Cartesian product of phase distributions per class
        kstate = _cartesian_space_for_phases(si, phases, R)

        # Build full states: [mi_buf, kstate]
        for ks in kstate:
            state = np.concatenate([mi_buf, ks])
            all_states.append(state)

    if not all_states:
        return np.zeros((1, 1 + int(sum(phases))))

    # Stack and ensure consistent column count
    space = np.array(all_states)

    # Sort and remove duplicates
    space = np.unique(space, axis=0)

    return space


def _cartesian_space_preempt(n: np.ndarray, phases: np.ndarray, S: int, R: int) -> np.ndarray:
    """
    Generate state space for the preempt family (LCFSPR, LCFSPRPRIO, LCFSPI,
    LCFSPIPRIO, FCFSPR, FCFSPRPRIO, FCFSPI, FCFSPIPRIO).

    Tracks an ordered buffer of the preempted jobs and the service phases. The
    state format matches MATLAB State.fromMarginal (the LCFSPR/FCFSPRPRIO case):
    - Buffer: interleaved [class, phase] pairs, one pair per buffer position
      (class and phase are 1-based; 0 marks an empty pair). Storing the phase is
      what distinguishes the preempt family from plain FCFS/LCFS: a PR job
      restarts in the phase it was preempted at. A PI job discards that phase
      and restarts from pie, but the pair layout is shared so that the same
      enumerated space serves both.
    - Server phases: phase distribution per class for jobs in service.

    Args:
        n: Array of job counts per class
        phases: Array of number of phases per class
        S: Number of servers
        R: Number of classes

    Returns:
        State space matrix where each row is [buf_pairs..., phase_counts...]
    """
    n = np.atleast_1d(n).astype(int)
    phases = np.atleast_1d(phases).astype(int)
    total_jobs = int(sum(n))
    sum_phases = int(sum(phases))

    if total_jobs == 0:
        # One empty [class, phase] pair + service phases
        return np.zeros((1, 2 + sum_phases))

    # Build list of job classes with repetition (1-based class IDs)
    vi = []
    for r in range(R):
        vi.extend([r + 1] * int(n[r]))

    mi = multiset_perms(vi)
    mi = np.array(mi) if len(mi) > 0 else np.array([vi])
    if mi.ndim == 1:
        mi = mi.reshape(1, -1)

    num_in_service = int(min(total_jobs, S))
    num_in_buffer = int(total_jobs - num_in_service)

    all_states = []
    for perm_idx in range(mi.shape[0]):
        perm = mi[perm_idx]

        # Buffer: first (total_jobs - S) positions; server: last S positions
        if num_in_buffer > 0:
            mi_buf = perm[:num_in_buffer]
        else:
            mi_buf = np.zeros(1, dtype=int)  # single empty pair
        mi_srv = perm[max(0, num_in_buffer):]

        # Count jobs of each class in service
        si = np.zeros(R, dtype=int)
        for class_id in mi_srv:
            if class_id > 0:
                si[int(class_id) - 1] += 1

        # Service phase distributions
        kstate = _cartesian_space_for_phases(si, phases, R)

        # Phase combinations for the buffered (preempted) jobs: each buffered
        # job of class j can have been preempted in any of its phases.
        bkstate = None
        for j in mi_buf:
            if j > 0:
                col = np.arange(1, int(phases[int(j) - 1]) + 1).reshape(-1, 1)
            else:
                col = np.array([[0]])
            bkstate = col if bkstate is None else _cartesian_product(bkstate, col)
        if bkstate is None:
            bkstate = np.zeros((1, len(mi_buf)))

        nbuf = len(mi_buf)
        for bphase_row in bkstate:
            # Interleave class/phase: even cols = class, odd cols = phase
            bufstate = np.zeros(2 * nbuf)
            bufstate[0::2] = mi_buf
            bufstate[1::2] = bphase_row
            for ks in kstate:
                all_states.append(np.concatenate([bufstate, ks]))

    if not all_states:
        return np.zeros((1, 2 + sum_phases))

    space = np.array(all_states)
    space = np.unique(space, axis=0)
    return space


def _generate_siro_space(n: np.ndarray, phases: np.ndarray, S: int, R: int) -> np.ndarray:
    """Generate state space for SIRO (Service In Random Order) scheduling.

    SIRO is work-conserving: when the total number of jobs present exceeds the
    number of servers S (the only case this helper is invoked for — the caller
    handles ``sum(n) <= S`` separately), exactly S jobs occupy the servers and
    the remainder wait in an *unordered* buffer. Only the choice of which
    classes hold the S service slots varies; the buffer is a per-class count.
    Iterating ``num_in_service`` over ``0..S`` (as an earlier version did) would
    emit non-work-conserving states with an idle server while jobs wait, which
    are spurious, reachable, and depress throughput.
    """
    space = []

    total_jobs = int(sum(n))
    num_in_service = min(total_jobs, S)
    # Distribute exactly num_in_service jobs among classes (bounded by n)
    for service_dist in _integer_partitions(num_in_service, R, n):
        buffer_dist = n - service_dist

        # Generate phase states for service portion
        phase_state = _cartesian_space_for_phases(service_dist, phases, R)

        # Add buffer info
        for ps in phase_state:
            state = np.concatenate([buffer_dist, ps])
            space.append(state)

    return np.array(space) if space else np.array([])


def _integer_partitions(total: int, num_parts: int, max_vals: np.ndarray = None):
    """Generate all ways to partition total into num_parts."""
    if num_parts == 1:
        if max_vals is not None and total <= max_vals[0]:
            yield np.array([total])
        elif max_vals is None:
            yield np.array([total])
        return

    if max_vals is None:
        max_vals = np.full(num_parts, total)

    for i in range(min(total, int(max_vals[0])) + 1):
        for rest in _integer_partitions(total - i, num_parts - 1, max_vals[1:]):
            yield np.concatenate([[i], rest])


def _cartesian_product(*arrays) -> np.ndarray:
    """
    Compute cartesian product of arrays.
    Each array is concatenated horizontally for each combination.
    """
    if not arrays:
        return np.array([])

    if len(arrays) == 1:
        return arrays[0]

    # Start with first array
    result = arrays[0]

    # Iteratively cartesian product with remaining arrays
    for arr in arrays[1:]:
        # Expand result and arr to compute cartesian product
        n_result = result.shape[0]
        n_arr = arr.shape[0]

        # Repeat result n_arr times
        result_expanded = np.repeat(result, n_arr, axis=0)

        # Tile arr n_result times
        arr_expanded = np.tile(arr, (n_result, 1))

        # Concatenate
        result = np.concatenate([result_expanded, arr_expanded], axis=1)

    return result


def fromMarginalAndRunning(sn, ind: int, n: Union[np.ndarray, list],
                           s: Union[np.ndarray, list], options: dict = None) -> np.ndarray:
    """
    Generate state space with specific marginal and running job counts.

    Creates states where node has n[r] jobs of class r total,
    with s[r] jobs of class r currently in service (running).

    Args:
        sn: NetworkStruct or Network object
        ind: Node index (0-based)
        n: Vector of total job counts per class
        s: Vector of running job counts per class
        options: Optional configuration dictionary

    Returns:
        State space matrix where each row is a valid state.
    """
    # Handle Network object input
    if hasattr(sn, 'get_struct'):
        sn = sn.get_struct()

    if options is None:
        options = {'force': False}

    n = np.atleast_1d(n).astype(int)
    s = np.atleast_1d(s).astype(int)

    R = sn.nclasses if hasattr(sn, 'nclasses') else len(n)

    # Ensure arrays have correct length
    if len(n) < R:
        n = np.concatenate([n, np.zeros(R - len(n), dtype=int)])
    if len(s) < R:
        s = np.concatenate([s, np.zeros(R - len(s), dtype=int)])

    # FJ-augmented struct: the join state is the per-class count vector of
    # buffered jobs/siblings, deterministic given the marginals (the started
    # counts s are immaterial, the Join has no service)
    if getattr(sn, 'isfjaugmented', False) and sn.nodetype[ind] == NodeType.JOIN:
        return np.asarray(n, dtype=float).reshape(1, -1)

    # Get station index
    if hasattr(sn, 'nodeToStation'):
        ist = sn.nodeToStation[ind]
    else:
        ist = ind

    # Get number of servers. Coerce to int, but keep inf for an
    # infinite-server (Delay/INF) station: int(inf) is an OverflowError, and
    # MATLAB's fromMarginalAndRunning reads sn.nservers(ist) unconverted, so
    # every model containing a Delay reached that error instead of a state.
    if hasattr(sn, 'nservers'):
        S = sn.nservers[ist]
        if np.isfinite(S):
            S = int(S)
    else:
        S = 1

    # Get phase information (same logic as fromMarginal)
    phases = np.zeros(R, dtype=int)
    if hasattr(sn, 'proc') and sn.proc is not None:
        for r in range(R):
            if hasattr(sn.proc[ist], '__getitem__') and r < len(sn.proc[ist]):
                proc_entry = sn.proc[ist][r]
                # sn.proc stores (D0, D1); the phase count is its order.
                from ..sn.proc_form import proc_n_phases
                if proc_entry is not None:
                    phases[r] = proc_n_phases(proc_entry)
                else:
                    phases[r] = 1
            else:
                phases[r] = 1
    else:
        phases = np.ones(R, dtype=int)

    # Check capacity constraints
    if hasattr(sn, 'classcap'):
        if np.any(n > sn.classcap[ist]):
            return np.array([])

    # Check running constraint: running jobs cannot exceed servers
    if S > 0 and sum(s) > S:
        return np.array([])

    # Get scheduling strategy
    if hasattr(sn, 'sched'):
        sched = sn.sched[ist]
    else:
        sched = SchedStrategy.FCFS

    # Generate state space based on scheduling strategy
    space = _generate_state_space_with_running(sched, S, n, s, phases, R)

    # Sort and unique the state space
    if len(space) > 0:
        space = np.unique(space, axis=0)
        # Reverse sort to put states with jobs in phase 1 earlier
        space = space[::-1]

    return space


def _generate_state_space_with_running(sched: int, S: int, n: np.ndarray, s: np.ndarray,
                                       phases: np.ndarray, R: int) -> np.ndarray:
    """
    Generate local state space with specific running job counts.

    Args:
        sched: SchedStrategy enum value
        S: Number of servers
        n: Total job counts per class
        s: Running job counts per class
        phases: Number of phases per class
        R: Number of classes

    Returns:
        State space matrix for this node
    """
    if sched in [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.DPS,
                 SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.LPS]:
        # For these strategies, all jobs are "in service"
        # Running = total for these schedulers
        return _cartesian_space_for_phases(s, phases, R)

    elif sched in _SCHED_TAG:
        # FCFS/LCFS: ordered buffer + service
        total_jobs = int(sum(n))
        running_jobs = int(sum(s))

        if total_jobs == 0:
            return np.zeros((1, 1 + int(sum(phases))))

        # Jobs in buffer = n - s per class
        buffer_jobs = n - s

        # Build list of job classes in buffer (with repetition)
        # Using 1-based class IDs to match MATLAB convention
        inbuf = []
        for r in range(R):
            if buffer_jobs[r] > 0:
                inbuf.extend([r + 1] * int(buffer_jobs[r]))  # 1-based class IDs

        # One buffer ordering is enough: this is an INITIAL state and the caller
        # keeps only the lexicographic maximum, whose buffer prefix is the
        # descending-sorted buffer. Enumerating every permutation first was
        # factorial in the buffer content. The phase distributions of the running
        # jobs still vary, so those rows are kept.
        if len(inbuf) > 0:
            mi_buf = np.array(sorted(inbuf, reverse=True), dtype=float)
        else:
            mi_buf = np.array([0], dtype=float)  # Empty buffer marker

        # kstate = cartesian product of phase distributions per class
        kstate = _cartesian_space_for_phases(s, phases, R)

        all_states = [np.concatenate([mi_buf, ks]) for ks in kstate]

        if not all_states:
            return np.zeros((1, 1 + int(sum(phases))))

        return np.array(all_states)

    elif sched in _SCHED_COUNT:
        # SIRO/SEPT/LEPT: unordered per-class-count buffer + service
        total_jobs = int(sum(n))

        if total_jobs <= S:
            # All jobs in service
            space = _cartesian_space_for_phases(s, phases, R)
            return np.concatenate([np.zeros((space.shape[0], R)), space], axis=1)
        else:
            # Buffer jobs = n - s
            mi_buf = n - s

            all_states = []
            # Generate phase distributions for running jobs
            kstate = _cartesian_space_for_phases(s, phases, R)

            for ks in kstate:
                state = np.concatenate([mi_buf, ks])
                all_states.append(state)

            return np.array(all_states) if all_states else np.array([])

    else:
        # Default: just track service phases
        return _cartesian_space_for_phases(s, phases, R)


def fromMarginalAndStarted(sn, ind: int, n: Union[np.ndarray, list],
                           s: Union[np.ndarray, list], options: dict = None) -> np.ndarray:
    """Wrapper appending the synchronous-call (REPLY) counter columns.

    The discipline branches below return from several places, so the counter
    columns are appended here, once, for every exit path. An initial or
    user-supplied state has no call outstanding, so the counters are zero --
    but the columns must be present, otherwise the row is narrower than the
    enumerated local space, matches no state, the unreachable-state pruning is
    silently skipped, and the enumerated-but-unreachable "counter set while
    every job is here" states remain as a second absorbing class.
    """
    if hasattr(sn, 'get_struct'):
        sn = sn.get_struct()
    space = _from_marginal_and_started(sn, ind, n, s, options)
    w = _reply_width(sn, ind)
    if w > 0:
        space = np.atleast_2d(np.asarray(space))
        if space.size > 0:
            space = np.hstack([space, np.zeros((space.shape[0], w))])
    return space


def _from_marginal_and_started(sn, ind: int, n: Union[np.ndarray, list],
                               s: Union[np.ndarray, list], options: dict = None) -> np.ndarray:
    """
    Generate state space with specific marginal and started job counts.

    Creates states where node has n[r] jobs of class r total,
    with s[r] jobs of class r that have started service.
    Started jobs are placed in phase 1 (the initial phase).

    Args:
        sn: NetworkStruct or Network object
        ind: Node index (0-based)
        n: Vector of total job counts per class
        s: Vector of started job counts per class
        options: Optional configuration dictionary

    Returns:
        State space matrix where each row is a valid state.
    """
    # Handle Network object input
    if hasattr(sn, 'get_struct'):
        sn = sn.get_struct()

    if options is None:
        options = {'force': True}

    n = np.atleast_1d(n).astype(int)
    s = np.atleast_1d(s).astype(int)

    R = sn.nclasses if hasattr(sn, 'nclasses') else len(n)

    # Ensure arrays have correct length
    if len(n) < R:
        n = np.concatenate([n, np.zeros(R - len(n), dtype=int)])
    if len(s) < R:
        s = np.concatenate([s, np.zeros(R - len(s), dtype=int)])

    # FJ-augmented struct: the join state is the per-class count vector of
    # buffered jobs/siblings, deterministic given the marginals (the started
    # counts s are immaterial, the Join has no service)
    if getattr(sn, 'isfjaugmented', False) and sn.nodetype[ind] == NodeType.JOIN:
        return np.asarray(n, dtype=float).reshape(1, -1)

    # Get station index
    if hasattr(sn, 'nodeToStation'):
        ist = sn.nodeToStation[ind]
    else:
        ist = ind

    # Get number of servers
    if hasattr(sn, 'nservers'):
        S = int(sn.nservers[ist])
    else:
        S = 1

    # Get phase information
    phases = np.zeros(R, dtype=int)
    if hasattr(sn, 'proc') and sn.proc is not None:
        for r in range(R):
            if hasattr(sn.proc[ist], '__getitem__') and r < len(sn.proc[ist]):
                proc_entry = sn.proc[ist][r]
                # sn.proc stores (D0, D1); the phase count is its order.
                from ..sn.proc_form import proc_n_phases
                if proc_entry is not None:
                    phases[r] = proc_n_phases(proc_entry)
                else:
                    phases[r] = 1
            else:
                phases[r] = 1
    else:
        phases = np.ones(R, dtype=int)

    # Check capacity constraints
    if hasattr(sn, 'classcap'):
        if np.any(n > sn.classcap[ist]):
            return np.array([])

    # Check started constraint: started jobs cannot exceed servers
    if S > 0 and sum(s) > S:
        return np.array([])

    # Get scheduling strategy
    if hasattr(sn, 'sched'):
        sched = sn.sched[ist]
    else:
        sched = SchedStrategy.FCFS

    # Pass-and-swap / order-independent: the local state is the ordered list of
    # class indices (1-based; oldest first). The started counts s are immaterial
    # (no server/buffer split). Mirrors fromMarginal.
    if sched == SchedStrategy.PAS:
        W = sn.cap[ist]
        if not np.isfinite(W):
            raise ValueError('PAS stations require finite capacity for state-space generation.')
        W = int(W)
        total = int(np.sum(n))
        if total == 0:
            return np.zeros((1, W))
        if total > W:
            return np.zeros((0, W))
        vi = []
        for r in range(R):
            vi.extend([r + 1] * int(n[r]))
        mi = np.array(multiset_perms(vi), dtype=float)
        if mi.ndim == 1:
            mi = mi.reshape(1, -1)
        pad = W - mi.shape[1]
        if pad > 0:
            mi = np.concatenate([mi, np.zeros((mi.shape[0], pad))], axis=1)
        return mi

    # Generate state space based on scheduling strategy
    space = _generate_state_space_with_started(sched, S, n, s, phases, R)

    # Sort and unique the state space
    if len(space) > 0:
        space = np.unique(space, axis=0)
        # Reverse sort to put states with jobs in phase 1 earlier
        space = space[::-1]

    return space


def _generate_state_space_with_started(sched: int, S: int, n: np.ndarray, s: np.ndarray,
                                       phases: np.ndarray, R: int) -> np.ndarray:
    """
    Generate local state space with specific started job counts.

    Started jobs are placed in phase 1 only (the initial phase).

    Args:
        sched: SchedStrategy enum value
        S: Number of servers
        n: Total job counts per class
        s: Started job counts per class
        phases: Number of phases per class
        R: Number of classes

    Returns:
        State space matrix for this node
    """
    if sched in [SchedStrategy.INF, SchedStrategy.PS, SchedStrategy.DPS,
                 SchedStrategy.GPS, SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO,
                 SchedStrategy.GPSPRIO, SchedStrategy.LPS]:
        # These policies track only the jobs in the servers, and under a shared
        # server EVERY job is in one, so the state holds the TOTAL count n[r],
        # not the started count s[r] (MATLAB fromMarginalAndStarted.m, the
        # INF/PS/DPS/GPS/LPS branch). Writing s[r] here lost the queued jobs.
        kstate = np.zeros((1, int(sum(phases))))
        col = 0
        for r in range(R):
            # Put all n[r] jobs in phase 1 (column col)
            kstate[0, col] = n[r]
            col += int(phases[r])
        return kstate

    elif sched in _SCHED_TAG:
        # FCFS/LCFS: ordered buffer + service
        total_jobs = int(sum(n))

        if total_jobs == 0:
            return np.zeros((1, 1 + int(sum(phases))))

        # Jobs in buffer = n - s per class
        buffer_jobs = n - s

        # Build list of job classes in buffer (with repetition)
        inbuf = []
        for r in range(R):
            if buffer_jobs[r] > 0:
                inbuf.extend([r + 1] * int(buffer_jobs[r]))  # 1-based class IDs

        # One buffer ordering is enough: this is an INITIAL state and the caller
        # keeps only the lexicographic maximum, i.e. the descending-sorted buffer.
        # Enumerating every permutation first was factorial in the buffer content.
        if len(inbuf) > 0:
            mi_buf = np.array([sorted(inbuf, reverse=True)], dtype=float)
        else:
            mi_buf = np.array([[0]])  # Empty buffer marker

        # For started, all jobs are in phase 1 (no phase enumeration)
        kstate = np.zeros((1, int(sum(phases))))
        col = 0
        for r in range(R):
            kstate[0, col] = s[r]  # All started jobs in phase 1
            col += int(phases[r])

        return np.array([np.concatenate([mi_buf[0], kstate[0]])])

    elif sched in _SCHED_COUNT or sched == SchedStrategy.POLLING:
        # SIRO/SEPT/LEPT/POLLING: unordered per-class-count buffer + service
        total_jobs = int(sum(n))

        # Buffer jobs = n - s
        mi_buf = n - s

        # Started jobs in phase 1
        kstate = np.zeros((1, int(sum(phases))))
        col = 0
        for r in range(R):
            kstate[0, col] = s[r]
            col += int(phases[r])

        state = np.concatenate([mi_buf, kstate[0]])
        return np.array([state])

    else:
        # Default: started jobs in phase 1
        kstate = np.zeros((1, int(sum(phases))))
        col = 0
        for r in range(R):
            kstate[0, col] = s[r]
            col += int(phases[r])
        return kstate


def toMarginalAggr(sn, ind: int, state_i: np.ndarray,
                   K: np.ndarray = None, Ks: np.ndarray = None,
                   space_buf: np.ndarray = None, space_srv: np.ndarray = None,
                   space_var: np.ndarray = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fast extraction of aggregated marginal (ni, nir) from state vector.

    Unlike toMarginal which also computes sir and kir, this only computes
    the total jobs per class (nir) and total jobs (ni), which is sufficient
    for most afterEvent rate computations.

    Port from JAR ToMarginal.java:toMarginalAggr().

    Args:
        sn: NetworkStruct
        ind: Node index (0-based)
        state_i: State vector(s), shape (n_states, n_cols)
        K: Phase counts per class, shape (R,)
        Ks: Phase shift per class, shape (R,)
        space_buf: Buffer portion of state
        space_srv: Server portion of state
        space_var: Variable portion of state

    Returns:
        Tuple of (ni, nir):
        - ni: Total jobs per state row, shape (n_states,)
        - nir: Jobs per class per state row, shape (n_states, R)
    """
    if hasattr(sn, 'get_struct'):
        sn = sn.get_struct()

    state_i = np.atleast_2d(state_i)
    n_states = state_i.shape[0]
    R = sn.nclasses

    # Join stations on FJ-augmented structs: plain per-class count state
    if getattr(sn, 'isfjaugmented', False) and sn.nodetype[ind] == NodeType.JOIN:
        nir = state_i[:, -R:]
        ni = np.sum(nir, axis=1)
        return ni, nir

    # Non-station stateful nodes (Router, Cache, Transition)
    if not sn.isstation[ind] and sn.isstateful[ind]:
        V = int(np.sum(sn.nvars[ind])) if sn.nvars is not None else 0
        n_cols = max(0, state_i.shape[1] - V)
        nir = np.zeros((n_states, max(1, n_cols)))
        for i in range(n_states):
            for j in range(min(n_cols, state_i.shape[1])):
                if j < nir.shape[1]:
                    nir[i, j] = state_i[i, j]
        ni = np.sum(nir, axis=1)
        return ni, nir

    ist = int(sn.nodeToStation[ind])

    # Pass-and-swap / order-independent: count classes from the ordered list.
    if sn.sched[ist] == SchedStrategy.PAS:
        Vp = int(np.sum(sn.nvars[ind])) if sn.nvars is not None else 0
        listcols = state_i[:, :state_i.shape[1] - Vp] if Vp > 0 else state_i
        nir = np.zeros((n_states, R))
        for r in range(R):
            nir[:, r] = np.sum(listcols == (r + 1), axis=1)
        ni = np.sum(nir, axis=1)
        return ni, nir

    # Default K, Ks from sn if not provided
    if K is None:
        K = np.array(sn.phasessz[ist], dtype=int)
    K = np.atleast_1d(K).astype(int)
    if Ks is None:
        Ks = np.array(sn.phaseshift[ist], dtype=int)
    Ks = np.atleast_1d(Ks).astype(int)

    sumK = int(np.sum(K))
    V = int(np.sum(sn.nvars[ind])) if sn.nvars is not None else 0

    # Decompose state if not provided
    if space_var is None:
        if V > 0:
            space_var = state_i[:, -V:]
        else:
            space_var = np.zeros((n_states, 0))

    if space_srv is None:
        srv_start = state_i.shape[1] - sumK - V
        srv_end = state_i.shape[1] - V
        space_srv = state_i[:, srv_start:srv_end]

    if space_buf is None:
        buf_end = state_i.shape[1] - sumK - V
        space_buf = state_i[:, :buf_end] if buf_end > 0 else np.zeros((n_states, 0))

    # Compute nir from server phases
    nir = np.zeros((n_states, R))
    for r in range(R):
        for k in range(K[r]):
            col = Ks[r] + k
            if col < space_srv.shape[1]:
                nir[:, r] += space_srv[:, col]

    # Add buffer contributions based on scheduling strategy
    sched = sn.sched[ist]

    if sched == SchedStrategy.EXT:
        nir[:] = np.inf

    elif sched in _SCHED_TAG:
        # Buffer stores 1-based class IDs
        if space_buf.size > 0 and not (space_buf.shape == (1, 1) and space_buf[0, 0] == 0):
            for r in range(R):
                if space_buf.ndim == 1:
                    nir[:, r] += np.sum(space_buf == (r + 1))
                else:
                    nir[:, r] += np.sum(space_buf == (r + 1), axis=1)

    elif sched in _SCHED_COUNT or sched == SchedStrategy.POLLING:
        # Buffer stores per-class job counts
        for r in range(R):
            if space_buf.ndim >= 2 and r < space_buf.shape[1]:
                nir[:, r] += space_buf[:, r]
            elif space_buf.ndim == 1 and r < len(space_buf):
                nir[:, r] += space_buf[r]

    elif sched in _SCHED_PREEMPT:
        # Preempt-resume: buffer stores interleaved [class, phase] pairs, so
        # the preempted-job class IDs are at the even columns.
        if space_buf.size > 0:
            buf2d = space_buf if space_buf.ndim == 2 else space_buf.reshape(1, -1)
            buf_classes = buf2d[:, 0::2]
            for r in range(R):
                nir[:, r] += np.sum(buf_classes == (r + 1), axis=1)

    # Zero out disabled classes
    if hasattr(sn, 'rates') and sn.rates is not None:
        for r in range(R):
            if np.isnan(sn.rates[ist, r]):
                nt = sn.nodetype[ind]
                nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
                place_val = int(NodeType.PLACE.value) if hasattr(NodeType.PLACE, 'value') else int(NodeType.PLACE)
                if nt_val != place_val:
                    nir[:, r] = 0

    ni = np.sum(nir, axis=1)
    return ni, nir


def fromMarg(sn, ind: int, ntot: int, options: dict = None) -> np.ndarray:
    """
    Generate the state space with a given TOTAL queue length at a node.

    This is the class-summed counterpart of fromMarginal: where fromMarginal
    fixes how many jobs of EACH class the node holds, fromMarg fixes only how
    many jobs it holds ALTOGETHER, and returns the union of fromMarginal over
    every class split of ntot the node can hold.

    A class that is disabled at the station has classcap 0 and is excluded from
    the split enumeration up front rather than after the fact. Asking
    fromMarginal for a job of such a class yields an EMPTY local space, and an
    empty factor is absorbed by the cartesian product instead of annihilating
    it, so the job would silently disappear.

    Args:
        sn: NetworkStruct or Network object
        ind: Node index (0-based)
        ntot: Total number of jobs at the node, all classes summed
        options: Optional configuration dictionary

    Returns:
        State space matrix with the requested total
    """
    from ..pfqn import multichoose

    if hasattr(sn, 'get_struct'):
        sn = sn.get_struct()
    if options is None:
        options = {'force': False}

    R = sn.nclasses
    if ntot < 0:
        return np.zeros((0, 0))

    ccap = _class_caps(sn, ind, R)

    # Enumerate the class splits of ntot that the node can hold. ntot=0 has the
    # single empty split, which fromMarginal answers with the per-discipline
    # empty state; do not re-derive that width here.
    if ntot == 0:
        nset = np.zeros((1, R), dtype=int)
    else:
        nset = np.atleast_2d(multichoose(R, ntot))
        keep = np.all(nset <= ccap.reshape(1, -1), axis=1)
        nset = nset[keep, :]

    subspaces = []
    for j in range(nset.shape[0]):
        sj = np.atleast_2d(np.asarray(fromMarginal(sn, ind, nset[j, :], options)))
        if sj.size == 0:
            continue
        subspaces.append(sj)

    return _stack_left_padded(subspaces)


def fromMargAndStarted(sn, ind: int, ntot: int, stot: int, options: dict = None) -> np.ndarray:
    """
    Generate the states with a given TOTAL queue length and a given TOTAL
    number of started jobs.

    Where fromMarginalAndStarted takes one per-class vector n and one per-class
    vector s and builds ONE row, fromMargAndStarted takes only the two totals
    and returns the union of that row over every (n,s) pair consistent with
    them: sum(n)=ntot, sum(s)=stot, and s <= n elementwise.

    Classes disabled at the station are excluded from the enumeration through
    classcap, for the reason documented in fromMarg.

    Args:
        sn: NetworkStruct or Network object
        ind: Node index (0-based)
        ntot: Total number of jobs at the node, all classes summed
        stot: Total number of jobs that have started service
        options: Optional configuration dictionary

    Returns:
        State space matrix with the requested totals
    """
    from ..pfqn import multichoose, multichoosecon

    if hasattr(sn, 'get_struct'):
        sn = sn.get_struct()
    if options is None:
        options = {'force': True}

    R = sn.nclasses
    if ntot < 0 or stot < 0 or stot > ntot:
        return np.zeros((0, 0))

    ccap = _class_caps(sn, ind, R)

    if ntot == 0:
        nset = np.zeros((1, R), dtype=int)
    else:
        nset = np.atleast_2d(multichoose(R, ntot))
        keep = np.all(nset <= ccap.reshape(1, -1), axis=1)
        nset = nset[keep, :]

    subspaces = []
    for j in range(nset.shape[0]):
        nj = nset[j, :]
        # s must be drawn from the jobs actually present, which is what
        # multichoosecon expresses; s=0 is its one uncovered base case.
        if stot == 0:
            sset = np.zeros((1, R), dtype=int)
        else:
            sset = multichoosecon(nj, stot)
        for k in range(sset.shape[0]):
            sjk = np.atleast_2d(np.asarray(fromMarginalAndStarted(sn, ind, nj, sset[k, :], options)))
            if sjk.size == 0:
                continue
            subspaces.append(sjk)

    return _stack_left_padded(subspaces)


def _class_caps(sn, ind: int, R: int) -> np.ndarray:
    """Per-class capacity of a station, or unbounded at a non-station node."""
    isstation = getattr(sn, 'isstation', None)
    classcap = getattr(sn, 'classcap', None)
    if isstation is None or classcap is None or np.size(classcap) == 0:
        return np.full(R, np.inf)
    isst = np.asarray(isstation).ravel()
    if ind >= isst.size or not isst[ind]:
        return np.full(R, np.inf)
    ist = int(np.asarray(sn.nodeToStation).ravel()[ind])
    return np.asarray(classcap, dtype=float)[ist, :].ravel()


def _stack_left_padded(subspaces) -> np.ndarray:
    """
    Stack sub-spaces of different width, then unique and reverse.

    The buffer is RIGHT-aligned, so sub-spaces of different width must be
    padded on the LEFT before they are stacked, exactly as fromMarginal does
    for the reply-block sub-spaces. The reversal after the sort puts the empty
    state first and the states with jobs in phase 1 earlier.
    """
    if not subspaces:
        return np.zeros((0, 0))
    maxw = max(s.shape[1] for s in subspaces)
    padded = []
    for s in subspaces:
        if s.shape[1] < maxw:
            s = np.hstack([np.zeros((s.shape[0], maxw - s.shape[1])), s])
        padded.append(s)
    space = np.vstack(padded)
    space = np.unique(space, axis=0)
    return space[::-1, :]
