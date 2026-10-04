"""
Rate at which sibling tasks are discarded at each Join.

Port of matlab/src/api/fj/sn_join_droprate.m, mirrored by jline.api.sn.SnJoinDroprate.
"""

import numpy as np

from .sn_join_quorum import sn_join_quorum
from .sn_join_siblings import sn_join_siblings


def sn_join_droprate(sn, TN, AN):
    """
    A Join is the one station where the loss identity ArvR - Tput does NOT hold, because
    the two rates are in different units: AN counts the SIBLINGS offered to the join (N per
    parent job) while TN counts the PARENT jobs released by it (one per synchronisation).
    Reading ArvR - Tput there reports (N-1)/N of the offered traffic as lost at every join,
    standard joins included, when a standard join loses nothing at all.

    The siblings a join actually consumes are K per synchronisation, where K is the quorum
    (K = N on a standard join), so DropRateJoin = max(0, AN - K*TN), which is 0 for a
    standard join and (N-K)*TN for a quorum.

    This is the DERIVED value, exact given TN and AN. A solver that MEASURES the discards
    on its own sample path (SolverLDES) reports its own.

    Args:
        sn: the network structure
        TN: station throughputs, (nstations x nclasses)
        AN: station arrival rates, (nstations x nclasses)

    Returns:
        An (nstations x nclasses) array, zero away from the Join rows.
    """
    from ...lang.base import NodeType

    def _nt_val(nt):
        return int(nt.value) if hasattr(nt, 'value') else int(nt)

    M = int(sn.nstations)
    K = int(sn.nclasses)
    out = np.zeros((M, K))
    if TN is None or AN is None:
        return out
    fj = getattr(sn, 'fj', None)
    if fj is None or not np.any(np.asarray(fj)):
        return out
    TN = np.atleast_2d(np.asarray(TN, dtype=float))
    AN = np.atleast_2d(np.asarray(AN, dtype=float))
    join_val = _nt_val(NodeType.JOIN)
    nodeToStation = np.asarray(sn.nodeToStation).astype(int).flatten()
    for ind, nt in enumerate(sn.nodetype):
        if _nt_val(nt) != join_val:
            continue
        if ind >= len(nodeToStation):
            continue
        ist = int(nodeToStation[ind])
        if ist < 0 or ist >= M:
            continue
        for r in range(K):
            if ist >= AN.shape[0] or r >= AN.shape[1]:
                continue
            a = float(AN[ist, r])
            t = float(TN[ist, r]) if (ist < TN.shape[0] and r < TN.shape[1]) else 0.0
            if not np.isfinite(a) or not np.isfinite(t) or a <= 0:
                continue
            # PER CLASS: a variable forking level makes the sibling count differ
            # between classes, so it cannot be hoisted out of this loop.
            nsib = sn_join_siblings(sn, ind, r)
            if nsib <= 0:
                continue
            kreq = sn_join_quorum(sn, ind, r, nsib)
            out[ist, r] = max(0.0, a - kreq * t)
    return out
