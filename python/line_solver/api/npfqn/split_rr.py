"""
Deterministic (round-robin) traffic splitting degrees.

Port of matlab/src/api/npfqn/npfqn_traffic_split_rr.m.
"""

import numpy as np

from ...constants import RoutingStrategy, GlobalConstants


def npfqn_traffic_split_rr(sn) -> np.ndarray:
    """
    Deterministic (round-robin) split degree of the departure stream of each
    station-class.

    kRR[i, r] = k > 1 means that the class-r departures of station i are
    dispatched one-in-k by a round-robin node, so that a downstream flow
    carrying a fraction p of them is the k-fold convolution thinned with
    probability q = k*p and has SCV 1+p*(d2-k), against the Markovian
    1+p*(d2-1). kRR[i, r] = 1 marks an ordinary probabilistic split.

    Only two topologies admit the deterministic rule: the station dispatches
    round-robin itself, or it feeds with probability one a router that does
    and whose pointer no other flow advances. Anything else falls back to k=1.

    Args:
        sn: NetworkStruct

    Returns:
        (M, K) array of split degrees
    """
    M = sn.nstations
    K = sn.nclasses
    kRR = np.ones((M, K))

    routing = getattr(sn, 'routing', None)
    if routing is None or np.asarray(routing).size == 0:
        return kRR
    routing = np.asarray(routing)
    rr_val = int(RoutingStrategy.RROBIN.value) if hasattr(RoutingStrategy.RROBIN, 'value') \
        else int(RoutingStrategy.RROBIN)
    if not np.any(routing == rr_val):
        return kRR

    rtnodes = getattr(sn, 'rtnodes', None)
    rtnodes = np.asarray(rtnodes) if rtnodes is not None else np.zeros((0, 0))
    stationToNode = np.asarray(sn.stationToNode).flatten().astype(int)
    isstation = np.asarray(sn.isstation).flatten()

    for ist in range(M):
        ind = int(stationToNode[ist])
        for r in range(K):
            if int(routing[ind, r]) == rr_val:
                kRR[ist, r] = _rr_degree(sn, ind, r)
                continue
            if rtnodes.size == 0:
                continue
            # a sure transition into a router that dispatches round-robin
            dest = np.where(rtnodes[ind * K + r, :] > 0)[0]
            if dest.size != 1:
                continue
            dest = int(dest[0])
            jnd = dest // K
            s = dest - jnd * K
            if isstation[jnd] > 0 or int(routing[jnd, s]) != rr_val:
                continue
            if rtnodes[ind * K + r, dest] < 1 - GlobalConstants.FineTol:
                continue
            # the round-robin pointer must be advanced by this stream alone
            if np.count_nonzero(rtnodes[:, dest] > 0) != 1:
                continue
            kRR[ist, r] = _rr_degree(sn, jnd, s)

    return kRR


def _rr_degree(sn, ind, r) -> int:
    """Number of destinations the round-robin pointer of (ind, r) cycles through."""
    nodeparam = getattr(sn, 'nodeparam', None)
    if nodeparam is not None and ind in nodeparam:
        nparam = nodeparam[ind]
        ol = nparam.get('outlinks', None) if isinstance(nparam, dict) else getattr(nparam, 'outlinks', None)
        if ol is not None and r < len(ol) and ol[r] is not None and len(np.atleast_1d(ol[r])) > 0:
            return max(1, len(np.atleast_1d(ol[r])))
    connmatrix = getattr(sn, 'connmatrix', None)
    if connmatrix is not None:
        cm = np.asarray(connmatrix)
        if ind < cm.shape[0]:
            return max(1, int(np.count_nonzero(cm[ind, :] > 0)))
    return 1
