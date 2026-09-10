"""
MDD-rec: the normalising constant of a product-form model whose reachable set is
held in a decision diagram.

S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant for
product-form models of distributed systems with synchronisation", Future
Generation Computer Systems 111 (2020) 475-490, Sec. 4.

A product-form model has P(s) = (1/G) prod_k g_k(s_k) over its levels, and

    G = sum_{s in S} prod_k g_k(s_k).

Summing state by state is exponential and numerically unstable. MDD-rec instead
walks the diagram that already encodes S, accumulating the unnormalised mass of
each node ONCE (Def. 4.4, Algorithm 1):

    M(<l.p>) = sum_{v in S_l} g_l(v) * M(<l.p>[v]),   M(TRUE) = 1, M(FALSE) = 0

so the cost is O(sum_l |nodes_l| * |S_l|) rather than O(|S|), and G = M(root).

FORMALISM-AGNOSTIC. Nothing here knows what a level is: the paper's Appendix B
shows that on the lattice sum_k s_k = n of a closed queueing network this
collapses to Buzen's convolution, and Sec. 5 that on an S-invariant reachable
Petri net it collapses to the Coleman-Henderson-Taylor convolution (spn_conv).
Unlike either, it needs only that the reachable set be finite and encoded -- no
lattice, no S-invariant reachability.

THE MASK is how Sec. 5.3 computes measures. Restricting the sum at level l to a
subset of its local values gives the unnormalised mass of the corresponding
subset of S, so P(m_l = k) and P(e_j >= k) are the same recursion under a
different mask rather than three separate algorithms.

WHAT IS NOT HERE. The g_l themselves, and the test that the model has a product
form at all, are the caller's: the paper declares that out of scope (Sec. 3.2)
and refers to the per-formalism conditions instead. Passing g_l that do not
describe a product-form model returns a number that is not the normalising
constant of anything, and nothing here can detect it.

See also: mdd_reachset, spn_rec_enabled, spn_metrics, spn_conv.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

from typing import List, Optional, Sequence

import numpy as np

from ..io.logging import line_error
from .mdd import TERM_FALSE, TERM_TRUE


def _check(mdds, g, mask, caller):
    """Shape agreement between the diagram, the factors and the mask."""
    K = int(mdds.K)
    if len(g) != K:
        line_error(caller, 'one g_l per level is required')
    for l in range(K):
        if np.size(g[l]) != int(mdds.domain[l]):
            line_error(caller, 'g_l must have one entry per local state')
    if mask is not None:
        if len(mask) != K:
            line_error(caller, 'the mask must have one row per level')
        for l in range(K):
            if np.size(mask[l]) != int(mdds.domain[l]):
                line_error(caller, 'the mask must have one entry per local state')


def mdd_rec_masked(mdds, g: Sequence[Sequence[float]],
                   mask: Optional[Sequence[Sequence[bool]]] = None) -> float:
    """Unnormalised mass of the masked subset of the reachable set (Algorithm 1).

    Parameters
    ----------
    mdds : the reachable set, as exported by MDD.to_struct (level 0 is the root)
    g : g[l][v] is g_l(v), the per-level factor of the product form
    mask : per-level admissible local values; None admits everything, which is
        the plain MDD-rec of Algorithm 1 and returns G
    """
    _check(mdds, g, mask, 'mdd_rec')
    K = int(mdds.K)
    if int(mdds.root) == TERM_FALSE:
        return 0.0

    gl = [np.asarray(g[l], dtype=float).ravel() for l in range(K)]
    ml = None if mask is None else [np.asarray(mask[l], dtype=bool).ravel() for l in range(K)]
    memo: List[np.ndarray] = [np.zeros(int(mdds.nnodes[l])) for l in range(K)]
    done: List[np.ndarray] = [np.zeros(int(mdds.nnodes[l]), dtype=bool) for l in range(K)]

    # Explicit stack rather than recursion: a level count of a few hundred is
    # ordinary here and would otherwise meet the interpreter's recursion limit.
    stack = [(0, int(mdds.root))]
    while stack:
        l, nid = stack[-1]
        if done[l][nid - 1]:
            stack.pop()
            continue
        arcs = np.asarray(mdds.node[l][nid - 1]).ravel()
        pending = []
        acc = 0.0
        for v in range(int(mdds.domain[l])):
            if ml is not None and not ml[l][v]:
                continue
            gv = gl[l][v]
            if gv == 0.0:
                continue
            ch = int(arcs[v])
            if l + 1 == K:
                if ch != TERM_TRUE:
                    continue
                acc += gv
            else:
                if ch == TERM_FALSE:
                    continue
                if not done[l + 1][ch - 1]:
                    pending.append((l + 1, ch))
                else:
                    acc += gv * memo[l + 1][ch - 1]
        if pending:
            stack.extend(pending)
            continue
        memo[l][nid - 1] = acc
        done[l][nid - 1] = True
        stack.pop()
    return float(memo[0][int(mdds.root) - 1])


def mdd_rec(mdds, g: Sequence[Sequence[float]]) -> float:
    """The normalising constant G = sum_{s in S} prod_l g_l(s_l)."""
    return mdd_rec_masked(mdds, g, None)


def mdd_rec_marginal(mdds, g: Sequence[Sequence[float]], l: int) -> np.ndarray:
    """Unnormalised masses of {s in S : s_l = k}, one per local value k of level l.

    Divided by G these are P(m_l = k) of Sec. 5.3: the mean occupancy of a level
    is sum_k k * P(m_l = k), and its utilization 1 - P(m_l = 0).
    """
    K = int(mdds.K)
    l = int(l)
    if l < 0 or l >= K:
        line_error('mdd_rec_marginal', 'level index is out of range')
    d = int(mdds.domain[l])
    out = np.zeros(d)
    for k in range(d):
        mask = [np.ones(int(mdds.domain[j]), dtype=bool) for j in range(K)]
        mask[l][:] = False
        mask[l][k] = True
        out[k] = mdd_rec_masked(mdds, g, mask)
    return out
