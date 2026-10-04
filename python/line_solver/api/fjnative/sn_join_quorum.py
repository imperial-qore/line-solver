"""
Number of sibling tasks a Join node waits for in a given class.

Port of matlab/src/api/fj/sn_join_quorum.m, mirrored by jline.api.sn.SnJoinQuorum.
"""

from ...lang.base import JoinStrategy


def sn_join_quorum(sn, join_idx: int, r: int, nbranches: int) -> int:
    """
    A standard join, an absent declaration, a non-positive quorum and a quorum that is not
    smaller than the sibling count all return nbranches, i.e. the ordinary AND-join: those
    are the four ways a join fires only when every sibling has arrived.

    The count is the one the simulation engines apply (SolverLDES fixes it at FORK time and
    discards the stragglers when they reach the join), so an analytical solver reading it
    here charges the same synchronisation event.

    Args:
        sn: the network structure
        join_idx: node index of the Join
        r: class index the siblings are matched in
        nbranches: number of siblings the fork emits

    Returns:
        The number of siblings the join waits for.
    """
    nodeparam = getattr(sn, 'nodeparam', None)
    if nodeparam is None:
        return nbranches
    try:
        param = nodeparam[join_idx]
    except (KeyError, IndexError, TypeError):
        return nbranches
    if not isinstance(param, dict):
        return nbranches
    strategies = param.get('joinStrategy')
    required = param.get('joinRequired')
    if strategies is None or required is None:
        return nbranches
    if r >= len(strategies) or r >= len(required):
        return nbranches
    if strategies[r] == JoinStrategy.STD:
        return nbranches
    q = int(round(float(required[r])))
    if 0 < q < nbranches:
        return q
    return nbranches
