"""
Reachability set generation into a decision diagram.

After A.S. Miner, G. Ciardo, "Efficient Reachability Set Generation and Storage
Using Decision Diagrams", ICATPN 1999, LNCS 1639, pp.6-25.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

from typing import Callable, Sequence

import numpy as np

from .mdd import MDD


def mdd_reachset(domain: Sequence[int], init: Sequence[int],
                 nextfun: Callable[[tuple], object]) -> MDD:
    """Generate and store the reachability set into a quasi-reduced ordered MDD.

    Parameters
    ----------
    domain : length-K per-level local-state counts (values 0..domain[k]-1)
    init : length-K initial global state (0-based local values)
    nextfun : s -> T, an (m x K) array whose rows are the successors of s

    Returns
    -------
    MDD holding every state reachable from init.

    Notes
    -----
    This is the basic (explicit-frontier) realisation: a breadth-first search
    enumerates successors while the MDD provides the O(K) membership test that
    replaces the usual explicit visited hash. The stored set lives entirely in
    the MDD (O(#nodes) memory); only the transient BFS frontier is held
    explicitly. Symbolic image computation / saturation, which removes the
    explicit frontier too, is the natural next step but is out of scope here.
    """
    mdd = MDD(domain)
    init = tuple(int(v) for v in np.ravel(np.asarray(init)))
    mdd.insert(init)

    frontier = [init]
    head = 0
    while head < len(frontier):
        s = frontier[head]
        head += 1
        successors = nextfun(s)
        if successors is None:
            continue
        for row in successors:
            t = tuple(int(v) for v in row)
            if not mdd.member(t):
                mdd.insert(t)
                frontier.append(t)
        # drop already-expanded rows periodically to bound frontier memory
        if head > 1024 and 2 * head > len(frontier):
            frontier = frontier[head:]
            head = 0
    mdd.compact()   # reclaim the dead nodes left by the append-only inserts
    return mdd
