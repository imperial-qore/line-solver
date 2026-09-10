"""
Cache node event handler for afterEvent dispatch.

Handles ARV, DEP, and READ events at Cache nodes. Cache nodes are
stateful non-station nodes that implement item caching with various
replacement policies (FIFO, LRU, SFIFO, RR).

Port from JAR AfterEventCache.java.
"""

import numpy as np
from ...constants import EventType, RoutingStrategy, GlobalConstants


def _cpos(m, i, j):
    """
    Compute absolute position in cache state vector for position j in list i.

    Args:
        m: Array of list capacities
        i: List index (0-based)
        j: Position within list (0-based)

    Returns:
        Absolute column index in spaceVar
    """
    pos = int(np.sum(m[:i])) if i > 0 else 0
    return pos + j


def _is_in_retrieval_system(var_flat, item, total_cache_capacity):
    """True iff item is currently being retrieved.

    The retrieval system is encoded as a per-item occupancy bitmap appended
    after the cache contents: column ``total_cache_capacity + item`` is
    non-zero iff that item is currently being retrieved.
    """
    col = total_cache_capacity + item
    if col >= len(var_flat):
        return False
    return var_flat[col] != 0


def _add_to_retrieval_system(var_flat, item, total_cache_capacity):
    """Mark item as in-flight by setting its occupancy bit."""
    col = total_cache_capacity + item
    if col < len(var_flat):
        var_flat[col] = 1


def _remove_from_retrieval_system(var_flat, item, total_cache_capacity):
    """Clear item's occupancy bit."""
    col = total_cache_capacity + item
    if col < len(var_flat):
        var_flat[col] = 0


def after_event_cache(sn, ind, event, job_class, R, space_buf, space_srv, space_var,
                      is_simulation=False):
    """
    Handle events at a Cache node.

    Args:
        sn: NetworkStruct
        ind: Node index (0-based)
        event: EventType (ARV, DEP, or READ)
        job_class: Job class index (0-based)
        R: Number of classes
        space_buf: Buffer state (empty for Cache)
        space_srv: Server state, shape (n_rows, R) - job counts per class
        space_var: Variable state - cache list contents
        is_simulation: sample-path mode, which has no enumerated state space and
            hence no delayed-hit truncation level

    Returns:
        Tuple of (outspace, outrate, outprob)
    """
    space_srv = space_srv.copy()
    space_var = space_var.copy()
    n_rows = space_srv.shape[0]

    if event == EventType.ARV:
        # Increment job count for arriving class
        space_srv[:, job_class] += 1
        outspace = np.hstack([space_srv, space_var]) if space_var.size > 0 else space_srv.copy()
        outrate = -1.0 * np.ones((n_rows, 1))
        outprob = np.ones((n_rows, 1))
        return outspace, outrate, outprob

    elif event == EventType.DEP:
        if space_srv[0, job_class] > 0:
            # A retrieval-class job departs the cache only to BEGIN a retrieval
            # (cache -> queue). The per-item occupancy bit is already set by the
            # READ that started the fetch, so the departure only moves the job.
            space_srv[:, job_class] -= 1

            # Update round-robin pointer if applicable
            routing_val = sn.routing[ind, job_class] if hasattr(sn, 'routing') and sn.routing is not None else None
            rs_val = int(routing_val.value) if hasattr(routing_val, 'value') else int(routing_val) if routing_val is not None else -1
            rrobin_val = int(RoutingStrategy.RROBIN.value) if hasattr(RoutingStrategy.RROBIN, 'value') else int(RoutingStrategy.RROBIN)

            if rs_val == rrobin_val:
                nvar_sum = int(np.sum(sn.nvars[ind, :R + job_class + 1]))
                space_var_idx = nvar_sum - 1
                nparam = sn.nodeparam.get(ind, None) if isinstance(sn.nodeparam, dict) else sn.nodeparam[ind] if sn.nodeparam is not None else None
                if nparam is not None:
                    outlinks = nparam.get('outlinks', None) if isinstance(nparam, dict) else getattr(nparam, 'outlinks', None)
                    if outlinks is not None and job_class < len(outlinks):
                        ol = np.atleast_1d(outlinks[job_class])
                        current_val = space_var[0, space_var_idx] if space_var.ndim >= 2 else space_var[space_var_idx]
                        idx = -1
                        for i in range(len(ol)):
                            if ol[i] == current_val:
                                idx = i
                                break
                        next_val = ol[idx + 1] if idx < len(ol) - 1 else ol[0]
                        if space_var.ndim >= 2:
                            space_var[:, space_var_idx] = next_val
                        else:
                            space_var[space_var_idx] = next_val

            outspace = np.hstack([space_srv, space_var]) if space_var.size > 0 else space_srv.copy()
            outrate = GlobalConstants.Immediate * np.ones((n_rows, 1))
            outprob = np.ones((n_rows, 1))
            return outspace, outrate, outprob

        # No job to depart
        return np.zeros((0, 0)), np.zeros((0, 0)), np.zeros((0, 0))

    elif event == EventType.READ:
        return _handle_read(sn, ind, job_class, R, space_srv, space_var, is_simulation)

    return np.zeros((0, 0)), np.zeros((0, 0)), np.zeros((0, 0))


def _handle_read(sn, ind, job_class, R, space_srv, space_var, is_simulation=False):
    """Handle READ event at a Cache node."""
    # Extract cache parameters
    nparam = sn.nodeparam[ind] if sn.nodeparam is not None and ind in sn.nodeparam else None
    if nparam is None:
        return np.zeros((0, 0)), np.zeros((0, 0)), np.zeros((0, 0))

    if isinstance(nparam, dict):
        n = nparam.get('nitems', 0)
        m = np.atleast_1d(nparam.get('itemcap', []))
        ac = nparam.get('accost', None)
        hitclass = np.atleast_1d(nparam.get('hitclass', []))
        missclass = np.atleast_1d(nparam.get('missclass', []))
        pread = nparam.get('pread', None)
        # The native struct stores the policy under 'replacestrat'; the JAR
        # wrapper exposes it as 'replacement'. Without the 'replacestrat' read
        # every policy silently defaulted to FIFO (RR coincides numerically, but
        # LRU did not), so multi-list LRU caches were solved as FIFO.
        replacement = nparam.get('replacestrat', nparam.get('replacement', 'FIFO'))
        qadm = nparam.get('qlru', 1.0)
        total_cache_capacity = nparam.get('total_cache_capacity', None)
        retrieval_classes = nparam.get('retrieval_classes', None)
        retrieval_class_indices = nparam.get('retrieval_class_indices', set())
    else:
        n = getattr(nparam, 'nitems', 0)
        m = np.atleast_1d(getattr(nparam, 'itemcap', []))
        ac = getattr(nparam, 'accost', None)
        hitclass = np.atleast_1d(getattr(nparam, 'hitclass', []))
        missclass = np.atleast_1d(getattr(nparam, 'missclass', []))
        pread = getattr(nparam, 'pread', None)
        # 'replacestrat' is the native attribute; 'replacement' is the JAR alias.
        # Defaulting to FIFO here silently mis-solved LRU caches as FIFO.
        replacement = getattr(nparam, 'replacestrat', None)
        if replacement is None:
            replacement = getattr(nparam, 'replacement', 'FIFO')
        qadm = getattr(nparam, 'qlru', 1.0)
        total_cache_capacity = getattr(nparam, 'total_cache_capacity', None)
        retrieval_classes = getattr(nparam, 'retrieval_classes', None)
        retrieval_class_indices = getattr(nparam, 'retrieval_class_indices', set())

    m = m.astype(int)
    h = len(m)  # number of cache lists
    if total_cache_capacity is None:
        total_cache_capacity = int(np.sum(m))
    tcc = int(total_cache_capacity)

    # Default access-cost matrix when none is supplied. Mirrors MATLAB
    # sanitize.m (accessProb = diag(ones(1,nLevels),1) with bottom-right=1) and
    # the MVA/NC path in cache_retrieval_inputs.py: a miss enters list 1 and a
    # hit in list i is promoted to list i+1 (the terminal list stays put), each
    # (h+1)x(h+1) row summing to 1. Without this, _get_ac falls back to 1.0 for
    # every entry, which for multi-list caches (h>1) spuriously promotes a hit
    # into *every* higher list and inserts a miss into *every* list, mis-scaling
    # the hit/miss departure rates (single-list caches are unaffected because
    # the terminal-list hit path does not consult the access-cost matrix).
    if ac is None and h > 0 and n > 0:
        Rmat = np.zeros((h + 1, h + 1))
        for j in range(h):
            Rmat[j, j + 1] = 1.0
        Rmat[h, h] = 1.0
        ac = [[Rmat for _ in range(int(n))] for _ in range(R)]
    if retrieval_classes is not None:
        retrieval_classes = np.atleast_2d(np.asarray(retrieval_classes, dtype=int))
        if retrieval_classes.size == 0:
            retrieval_classes = None
    hitclass = hitclass.astype(int)
    missclass = missclass.astype(int)

    # Get replacement strategy name
    if hasattr(replacement, 'name'):
        repl_name = replacement.name
    elif hasattr(replacement, 'value'):
        repl_name = str(replacement.value)
    else:
        repl_name = str(replacement)

    # Check precondition: exactly one job at cache, and it's of this class
    srv_sum = np.sum(space_srv[0])
    if space_srv[0, job_class] <= 0 or int(srv_sum) != 1:
        return np.zeros((0, 0)), np.zeros((0, 0)), np.zeros((0, 0))

    # Get read probabilities for this class
    p = pread[job_class] if pread is not None and job_class < len(pread) else None
    if p is None:
        return np.zeros((0, 0)), np.zeros((0, 0)), np.zeros((0, 0))
    p = list(p)

    out_srv_list = []
    out_var_list = []
    out_rate_list = []
    Immediate = GlobalConstants.Immediate

    # Block B of the local-variable vector: per-retrieval-class counts of the
    # secondary requests merged onto an in-flight fetch (see _space_cache). Its
    # width and truncation level are read off the node state space rather than
    # from nodeparam, because only sn.space is propagated back from the
    # state-space generator.
    from .ctmc_ssg import cache_retrieval_class_map
    rc_list, rc_items, rc_orig_class = cache_retrieval_class_map(sn, ind)
    block_b_offset = tcc + n
    var_width = space_var.shape[1] if space_var.ndim >= 2 else space_var.shape[0]
    width_b = len(rc_list) if (var_width - block_b_offset) == len(rc_list) else 0
    # simulation has no enumerated state space, hence no truncation: a fetch may
    # merge any number of secondary requests
    max_pending = np.inf if is_simulation else 0
    if width_b > 0 and not is_simulation:
        isfc = int(sn.nodeToStateful[ind])
        spc = sn.space[isfc] if getattr(sn, 'space', None) is not None and isfc in sn.space else None
        if spc is not None and np.asarray(spc).size > 0:
            spc = np.atleast_2d(np.asarray(spc))
            max_pending = int(np.max(np.sum(spc[:, spc.shape[1] - width_b:], axis=1)))

    # Find enabled rows
    for e in range(space_srv.shape[0]):
        if space_srv[e, job_class] <= 0:
            continue

        # A retrieval-class READ (the retrieval returning to the cache) always
        # completes the miss that started the retrieval, regardless of the
        # (possibly unreachable) enumerated state.
        is_from_retrieval = job_class in retrieval_class_indices

        # CTMC mode: iterate over ALL items
        for k in range(n):
            space_srv_e = space_srv[e].copy()
            space_srv_e[job_class] -= 1  # Remove the reading job
            var = space_var[e].copy()

            # Find item k+1 in the cache-contents region only (1..tcc); the
            # trailing retrieval-system slots are out of scope for the hit test.
            posk = -1
            for col in range(min(tcc, var.shape[0] if var.ndim == 1 else var.shape[1])):
                v = var[col] if var.ndim == 1 else var[0, col]
                if v == k + 1:  # Items are 1-indexed
                    posk = col
                    break
            if is_from_retrieval:
                posk = -1  # retrieval-complete always takes the miss branch

            if posk == -1:
                # CACHE MISS or RETRIEVAL begin/return.
                rc = -1
                if (retrieval_classes is not None
                        and k < retrieval_classes.shape[0]
                        and job_class < retrieval_classes.shape[1]):
                    rc = int(retrieval_classes[k, job_class])
                var_flat = np.atleast_1d(var).ravel()

                # A returning retrieval that is not recorded in the retrieval-system
                # bitmap is an unreachable event-loop artifact; do not continue it.
                if is_from_retrieval and not _is_in_retrieval_system(var_flat, k, tcc):
                    continue

                # Begin a retrieval: the job switches to the retrieval class for
                # item k and item k is marked as being fetched. A concurrent
                # request for an item already being fetched is a delayed hit: it
                # merges onto the in-flight fetch and is held in block B until
                # that fetch completes.
                if (not is_from_retrieval) and rc >= 0:
                    space_srv_e_b = space_srv_e.copy()
                    var_b = var.copy()
                    vb_flat = np.atleast_1d(var_b).ravel()
                    if not _is_in_retrieval_system(var_flat, k, tcc):
                        space_srv_e_b[rc] += 1
                        _add_to_retrieval_system(vb_flat, k, tcc)
                    else:
                        bslot = rc_list.index(rc) if rc in rc_list else -1
                        bcol = block_b_offset + bslot
                        if bslot < 0 or bcol >= len(vb_flat) or \
                                np.sum(vb_flat[block_b_offset:]) >= max_pending:
                            continue  # beyond the delayed-hit truncation level
                        vb_flat[bcol] += 1
                    out_srv_list.append(space_srv_e_b)
                    out_var_list.append(vb_flat.reshape(var_b.shape))
                    out_rate_list.append(p[k] * Immediate)
                    continue

                # Item has now been retrieved: mark as a miss and clear its bit.
                # Every secondary request merged onto this fetch is released in the
                # same transition and departs as a delayed hit, in the hit class of
                # the job class that issued it.
                if var.size > tcc:
                    var = var.copy()
                    vf = np.atleast_1d(var).ravel()
                    _remove_from_retrieval_system(vf, k, tcc)
                    for bslot in range(width_b):
                        if rc_items[bslot] != k + 1:
                            continue
                        bcol = block_b_offset + bslot
                        if bcol < len(vf) and vf[bcol] > 0:
                            hc = int(hitclass[rc_orig_class[bslot]])
                            if hc >= 0:
                                space_srv_e[hc] += vf[bcol]
                            vf[bcol] = 0
                    var = vf.reshape(var.shape)
                space_srv_e[int(missclass[job_class])] += 1
                _handle_miss(m, h, k, p, ac, job_class, repl_name, space_srv_e, var,
                             out_srv_list, out_var_list, out_rate_list, Immediate, qadm)

            elif posk < int(np.sum(m) - m[h - 1]):
                # CACHE HIT in list i < h (not the last list)
                space_srv_e[int(hitclass[job_class])] += 1
                # Find which list i the item is in
                mcumsum = np.cumsum(m)
                i_list = 0
                for col in range(len(mcumsum)):
                    if posk < mcumsum[col]:
                        i_list = col
                        break
                # Position within list i
                j_pos = posk - (int(np.sum(m[:i_list])) if i_list > 0 else 0)

                _handle_hit_nonterminal(m, h, k, p, ac, job_class, repl_name, i_list, j_pos,
                                        space_srv_e, var,
                                        out_srv_list, out_var_list, out_rate_list, Immediate)

            else:
                # CACHE HIT in last list h
                space_srv_e[int(hitclass[job_class])] += 1
                j_pos = posk - (int(np.sum(m[:h - 1])) if h > 1 else 0)

                _handle_hit_terminal(m, h, k, p, repl_name, j_pos,
                                     space_srv_e, var,
                                     out_srv_list, out_var_list, out_rate_list, Immediate)

    if len(out_srv_list) == 0:
        return np.zeros((0, 0)), np.zeros((0, 0)), np.zeros((0, 0))

    outspace = np.hstack([np.array(out_srv_list), np.array(out_var_list)])
    outrate = np.array(out_rate_list).reshape(-1, 1)
    outprob = np.ones((len(out_rate_list), 1))
    return outspace, outrate, outprob


def _handle_miss(m, h, k, p, ac, job_class, repl_name, space_srv_e, var,
                 out_srv_list, out_var_list, out_rate_list, Immediate, qadm=1.0):
    """Handle cache miss: insert item k+1 into cache."""
    # HLRU shares LRU dynamics.
    if repl_name == 'HLRU':
        repl_name = 'LRU'

    if repl_name == 'QLRU':
        # q-LRU: pass-through mass = structural reject + (1-q) non-admission.
        ac_reject = _get_ac(ac, job_class, k, 0, 0)
        rejw = ac_reject + (1.0 - qadm) * (1.0 - ac_reject)
        if rejw > 0:
            out_srv_list.append(space_srv_e.copy())
            out_var_list.append(var.copy())
            out_rate_list.append(rejw * p[k] * Immediate)
        for l in range(h):
            head_pos = _cpos(m, l, 0)
            tail_pos = _cpos(m, l, int(m[l]) - 1)
            varp = var.copy()
            if m[l] > 1:
                varp[head_pos + 1:tail_pos + 1] = var[head_pos:tail_pos]
            varp[head_pos] = k + 1  # LRU head insert
            out_srv_list.append(space_srv_e.copy())
            out_var_list.append(varp)
            ac_val = _get_ac(ac, job_class, k, 0, l + 1)
            out_rate_list.append(qadm * ac_val * p[k] * Immediate)
        return

    # Cache reject (accost column 0): pass through without caching
    ac_reject = _get_ac(ac, job_class, k, 0, 0)
    if ac_reject > 0:
        out_srv_list.append(space_srv_e.copy())
        out_var_list.append(var.copy())
        out_rate_list.append(ac_reject * p[k] * Immediate)
    for l in range(h):
        # l is the target list index (0-based)
        listidx = l
        head_pos = _cpos(m, listidx, 0)
        tail_pos = _cpos(m, listidx, int(m[listidx]) - 1)

        if repl_name in ('FIFO', 'LRU', 'SFIFO'):
            varp = var.copy()
            # Shift items right by one (evict tail)
            if m[listidx] > 1:
                varp[head_pos + 1:tail_pos + 1] = var[head_pos:tail_pos]
            varp[head_pos] = k + 1  # Insert at head

            out_srv_list.append(space_srv_e.copy())
            out_var_list.append(varp)
            # Rate = ac[class][k](0, l+1) * p[k] * Immediate
            # ac is indexed ac[class][item] with rows: 0=miss, 1+i=hit-in-list-i
            ac_val = _get_ac(ac, job_class, k, 0, l + 1)
            out_rate_list.append(ac_val * p[k] * Immediate)

        elif repl_name == 'RR':
            # Random replacement: place at each position
            for r_pos in range(int(m[listidx])):
                varp = var.copy()
                varp[head_pos + r_pos] = k + 1

                out_srv_list.append(space_srv_e.copy())
                out_var_list.append(varp)
                ac_val = _get_ac(ac, job_class, k, 0, l + 1)
                out_rate_list.append(ac_val * p[k] / m[listidx] * Immediate)


def _handle_hit_nonterminal(m, h, k, p, ac, job_class, repl_name, i_list, j_pos,
                             space_srv_e, var,
                             out_srv_list, out_var_list, out_rate_list, Immediate):
    """Handle cache hit in non-terminal list i < h: promote item."""
    # HLRU/QLRU share LRU hit dynamics.
    if repl_name in ('HLRU', 'QLRU'):
        repl_name = 'LRU'

    for inew in range(i_list, h):
        varp = var.copy()

        if repl_name == 'FIFO':
            # Replace position (i, j) with evicted item from target list tail
            varp[_cpos(m, i_list, j_pos)] = var[_cpos(m, inew, int(m[inew]) - 1)]
            # Shift target list right, insert item at head
            if m[inew] > 1:
                head_inew = _cpos(m, inew, 0)
                tail_inew = _cpos(m, inew, int(m[inew]) - 1)
                varp[head_inew + 1:tail_inew + 1] = var[head_inew:tail_inew]
            varp[_cpos(m, inew, 0)] = k + 1

        elif repl_name == 'RR':
            # Random replacement: generate all swap possibilities
            for r_pos in range(int(m[inew])):
                varp2 = var.copy()
                varp2[_cpos(m, i_list, j_pos)] = var[_cpos(m, inew, r_pos)]
                varp2[_cpos(m, inew, r_pos)] = k + 1
                out_srv_list.append(space_srv_e.copy())
                out_var_list.append(varp2)
                ac_val = _get_ac(ac, job_class, k, 1 + i_list, 1 + inew)
                out_rate_list.append(ac_val * p[k] / m[inew] * Immediate)
            continue  # Already appended all states for this inew

        elif repl_name in ('LRU', 'SFIFO'):
            # Shift items in current list i to free position j
            if j_pos > 0:
                head_i = _cpos(m, i_list, 0)
                varp[head_i + 1:head_i + j_pos + 1] = var[head_i:head_i + j_pos]
            # Fill freed slot with evicted item from target list tail
            varp[_cpos(m, i_list, 0)] = var[_cpos(m, inew, int(m[inew]) - 1)]
            # Shift target list right, insert item at head
            if m[inew] > 1:
                head_inew = _cpos(m, inew, 0)
                tail_inew = _cpos(m, inew, int(m[inew]) - 1)
                varp[head_inew + 1:tail_inew + 1] = var[head_inew:tail_inew]
            varp[_cpos(m, inew, 0)] = k + 1

        else:
            continue

        out_srv_list.append(space_srv_e.copy())
        out_var_list.append(varp)
        ac_val = _get_ac(ac, job_class, k, 1 + i_list, 1 + inew)
        out_rate_list.append(ac_val * p[k] * Immediate)


def _handle_hit_terminal(m, h, k, p, repl_name, j_pos,
                          space_srv_e, var,
                          out_srv_list, out_var_list, out_rate_list, Immediate):
    """Handle cache hit in terminal list h."""
    if repl_name in ('FIFO', 'SFIFO', 'RR'):
        # No cache state change
        out_srv_list.append(space_srv_e.copy())
        out_var_list.append(var.copy())
        out_rate_list.append(p[k] * Immediate)

    elif repl_name in ('LRU', 'HLRU', 'QLRU'):
        # Promote item to head of last list
        varp = var.copy()
        head_h = _cpos(m, h - 1, 0)
        if j_pos > 0:
            varp[head_h + 1:head_h + j_pos + 1] = var[head_h:head_h + j_pos]
        varp[head_h] = var[_cpos(m, h - 1, j_pos)]

        out_srv_list.append(space_srv_e.copy())
        out_var_list.append(varp)
        out_rate_list.append(p[k] * Immediate)


def _get_ac(ac, job_class, item, row, col):
    """
    Get access cost matrix entry.

    ac[class][item] is an (h+1) x (h+1) matrix.
    Row 0 = miss-to-list probs; row 1+i = hit-in-list-i-to-list probs.
    """
    if ac is None:
        return 1.0
    try:
        if isinstance(ac, (list, tuple)):
            if job_class < len(ac) and item < len(ac[job_class]):
                mat = ac[job_class][item]
                if isinstance(mat, np.ndarray):
                    return float(mat[row, col]) if mat.ndim == 2 else float(mat[col])
                elif isinstance(mat, (list, tuple)):
                    return float(mat[row][col]) if isinstance(mat[0], (list, tuple)) else float(mat[col])
        elif isinstance(ac, np.ndarray):
            return float(ac[job_class, item, row, col]) if ac.ndim == 4 else float(ac[row, col])
    except (IndexError, TypeError):
        pass
    return 1.0
