"""
Number of sibling tasks forked per parent job on a fork-join pair.

Port of matlab/src/api/fj/sn_join_siblings.m, mirrored by jline.api.sn.SnJoinSiblings.
"""

import numpy as np


def sn_join_siblings(sn, join_idx: int, r=None) -> int:
    """
    Siblings are counted at the FORK, as the simulation engines count them.

    THE COUNT IS PER LINK, not the out-degree times a node-wide scalar. A fork carries a
    VARIABLE FORKING LEVEL: ``set_tasks_per_link(n, jobclass[, dest])`` sets one link of
    one class, ``set_tasks_per_link_distribution`` makes the degree a draw, and
    ``set_branch_prob`` makes a link taken with probability below one. All three land in
    ``sn.nodeparam[f]['fanOutLink']`` and ``['fanOutProb']``, both (nnodes x nclasses) and
    indexed by DESTINATION NODE, with the DISTRIBUTION case storing its mean, so

        N = sum_d fanOutProb(d,r) * fanOutLink(d,r)

    which is the EXPECTED sibling count and reduces to out-degree times the node-wide
    ``fanOut`` on a model that sets none of the three. Falls back to that older product
    when ``fanOutLink`` is absent -- it is built by the local-variable refresh, which runs
    after the capacity refresh -- and to the join's in-degree when the matched fork cannot
    be identified.

    N is what a quorum is measured against: the join fires on the k-th of N siblings and
    the remaining N-k are discarded when they arrive.

    Args:
        sn: the network structure
        join_idx: node index of the Join
        r: class index, or None for the widest fork over the classes

    Returns:
        The number of siblings forked per parent job, 0 when it cannot be determined.
    """
    conn = getattr(sn, 'connmatrix', None)
    if conn is None:
        return 0
    conn = np.atleast_2d(np.asarray(conn))
    if join_idx < 0 or join_idx >= conn.shape[1]:
        return 0
    n = int(np.count_nonzero(conn[:, join_idx]))
    fj = getattr(sn, 'fj', None)
    if fj is None:
        return n
    fj = np.atleast_2d(np.asarray(fj))
    if join_idx >= fj.shape[1]:
        return n
    forks = np.nonzero(fj[:, join_idx])[0]
    if forks.size == 0:
        return n
    f = int(forks[0])
    if f >= conn.shape[0]:
        return n

    param = None
    nodeparam = getattr(sn, 'nodeparam', None)
    if nodeparam is not None:
        try:
            param = nodeparam[f]
        except (KeyError, IndexError, TypeError):
            param = None
    if not isinstance(param, dict):
        param = None

    # The per-link count, when the refresh has built it.
    if param is not None and param.get('fanOutLink') is not None:
        fol = np.atleast_2d(np.asarray(param['fanOutLink'], dtype=float))
        fop = param.get('fanOutProb')
        if fop is not None:
            fop = np.atleast_2d(np.asarray(fop, dtype=float))
            if fop.shape != fol.shape:
                fop = None
        if fop is None:
            fop = (fol > 0).astype(float)
        if fol.size > 0:
            cols = range(fol.shape[1]) if r is None else (
                [int(r)] if 0 <= int(r) < fol.shape[1] else [])
            best = 0.0
            for c in cols:
                best = max(best, float(np.sum(fol[:, c] * fop[:, c])))
            if best > 0:
                return int(round(best))

    # Fallback: out-degree times the node-wide scalar.
    w = 1
    if param is not None:
        fan_out = param.get('fanOut')
        if fan_out is not None:
            fan_out = np.atleast_1d(np.asarray(fan_out, dtype=float))
            if fan_out.size > 0 and np.isfinite(fan_out.flat[0]):
                w = max(1, int(round(float(fan_out.flat[0]))))
    return int(np.count_nonzero(conn[f, :])) * w
