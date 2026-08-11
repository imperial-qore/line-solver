"""
Router node event handler for afterEvent dispatch.

Handles ARV and DEP events at Router nodes. Routers are stateful
non-station nodes that provide immediate pass-through with optional
round-robin routing pointer updates.

Port from JAR AfterEventRouter.java.
"""

import numpy as np
from ...constants import EventType, RoutingStrategy, GlobalConstants


def after_event_router(sn, ind, event, job_class, space_buf, space_srv, space_var):
    """
    Handle events at a Router node.

    Args:
        sn: NetworkStruct
        ind: Node index (0-based)
        event: EventType (ARV or DEP)
        job_class: Job class index (0-based)
        space_buf: Buffer state (empty for Router)
        space_srv: Server state, shape (n_rows, R) - job counts per class
        space_var: Variable state (round-robin pointers etc.)

    Returns:
        Tuple of (outspace, outrate, outprob):
        - outspace: np.ndarray of output states
        - outrate: np.ndarray of rates (-1 for ARV, Immediate for DEP)
        - outprob: np.ndarray of probabilities (all 1.0)
    """
    space_srv = space_srv.copy()
    space_var = space_var.copy()
    n_rows = space_srv.shape[0]

    if event == EventType.ARV:
        # Increment job count for arriving class
        space_srv[:, job_class] += 1
        # Compose output: [spaceSrv | spaceVar]
        outspace = np.hstack([space_srv, space_var]) if space_var.size > 0 else space_srv.copy()
        # Passive action: rate = -1
        outrate = -1.0 * np.ones((n_rows, 1))
        outprob = np.ones((n_rows, 1))
        return outspace, outrate, outprob

    elif event == EventType.DEP:
        if space_srv[0, job_class] > 0:
            # Decrement job count for departing class
            space_srv[:, job_class] -= 1

            # Update round-robin pointer if applicable. sn.routing values are
            # raw ints, so compare against the enum values (not the enum object).
            routing_val = _get_routing_strategy(sn, ind, job_class)
            _rv = int(routing_val.value) if hasattr(routing_val, 'value') else (int(routing_val) if routing_val is not None else -1)
            _rr = int(RoutingStrategy.RROBIN.value) if hasattr(RoutingStrategy.RROBIN, 'value') else int(RoutingStrategy.RROBIN)
            _wrr = int(RoutingStrategy.WRROBIN.value) if hasattr(RoutingStrategy.WRROBIN, 'value') else int(RoutingStrategy.WRROBIN)
            if _rv in (_rr, _wrr):
                _update_rrobin_pointer(sn, ind, job_class, space_var)

            # Compose output: [spaceSrv | spaceVar]
            outspace = np.hstack([space_srv, space_var]) if space_var.size > 0 else space_srv.copy()
            # Immediate action
            outrate = GlobalConstants.Immediate * np.ones((n_rows, 1))
            outprob = np.ones((n_rows, 1))
            return outspace, outrate, outprob

    # No valid transition: return empty
    outspace = np.zeros((0, space_srv.shape[1] + (space_var.shape[1] if space_var.ndim >= 2 else 0)))
    outrate = np.zeros((0, 1))
    outprob = np.zeros((0, 1))
    return outspace, outrate, outprob


def _get_routing_strategy(sn, ind, job_class):
    """Get routing strategy for a node/class pair."""
    if hasattr(sn, 'routing') and sn.routing is not None:
        val = sn.routing[ind, job_class]
        return val
    return None


def _update_rrobin_pointer(sn, ind, job_class, space_var):
    """Advance the round-robin routing pointer in space_var."""
    from ...constants import RoutingStrategy
    rr_val = int(RoutingStrategy.RROBIN.value) if hasattr(RoutingStrategy.RROBIN, 'value') else int(RoutingStrategy.RROBIN)
    wrr_val = int(RoutingStrategy.WRROBIN.value) if hasattr(RoutingStrategy.WRROBIN, 'value') else int(RoutingStrategy.WRROBIN)

    def _isrr(rt):
        v = int(rt.value) if hasattr(rt, 'value') else (int(rt) if rt is not None else -1)
        return v in (rr_val, wrr_val)

    # WRROBIN advances a POSITION pointer cyclically through the weighted-outlink
    # cycle (destinations repeated by weight), not a node index.
    from .routing_pointer import wrr_weighted_outlinks
    wol = wrr_weighted_outlinks(sn, ind, job_class)
    if wol is not None and len(wol) > 0:
        col = 0
        for rr in range(job_class):
            try:
                if _isrr(sn.routing[ind, rr]):
                    col += 1
            except (TypeError, KeyError, IndexError):
                pass
        ncol = space_var.shape[1] if space_var.ndim >= 2 else space_var.shape[0]
        if col >= ncol:
            return
        cyc = len(wol)
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
    outlinks = None
    nparam = sn.nodeparam[ind] if sn.nodeparam is not None and ind in sn.nodeparam else None
    if nparam is not None:
        ol = nparam.get('outlinks', None) if isinstance(nparam, dict) else getattr(nparam, 'outlinks', None)
        if ol is not None and job_class < len(ol):
            outlinks = np.atleast_1d(ol[job_class])
    if outlinks is None or len(outlinks) == 0:
        if hasattr(sn, 'connmatrix') and sn.connmatrix is not None:
            cm = np.asarray(sn.connmatrix)
            if ind < cm.shape[0]:
                outlinks = np.where(cm[ind, :] > 0)[0]
    if outlinks is None or len(outlinks) == 0:
        return

    # Pointer column within the variable block: one pointer per round-robin
    # class, ordered by class index.
    def _isrr(rt):
        v = int(rt.value) if hasattr(rt, 'value') else (int(rt) if rt is not None else -1)
        return v in (rr_val, wrr_val)
    space_var_idx = 0
    for rr in range(job_class):
        try:
            if _isrr(sn.routing[ind, rr]):
                space_var_idx += 1
        except (TypeError, KeyError, IndexError):
            pass
    ncol = space_var.shape[1] if space_var.ndim >= 2 else space_var.shape[0]
    if space_var_idx >= ncol:
        return

    # Find current position
    current_val = space_var[0, space_var_idx] if space_var.ndim >= 2 else space_var[space_var_idx]
    idx = -1
    for i in range(len(outlinks)):
        if outlinks[i] == current_val:
            idx = i
            break

    # Advance to next (wrap around)
    if idx < len(outlinks) - 1:
        next_val = outlinks[idx + 1]
    else:
        next_val = outlinks[0]

    if space_var.ndim >= 2:
        space_var[:, space_var_idx] = next_val
    else:
        space_var[space_var_idx] = next_val
