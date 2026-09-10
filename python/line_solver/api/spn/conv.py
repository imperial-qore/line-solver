"""
Convolution algorithm for the normalising constant of an S-invariant reachable
product-form stochastic Petri net.

J. Coleman, W. Henderson, P. Taylor, "Product form equilibrium distributions and
a convolution algorithm for stochastic Petri nets", Performance Evaluation 26(3),
1996, 159-180; presented as the point of comparison for MDD-rec in S. Balsamo,
A. Marin, I. Stojic, FGCS 111 (2020), Sec. 5.1.

With S the minimal-support S-invariant matrix and V = S m0 the load vector, the
reachability set of an S-INVARIANT REACHABLE net is exactly {m >= 0 : S m = V},
and conditioning on the marking of one place partitions it (Lemma 5.1). Writing
G_j(W) for the mass of the markings supported on the first j places with S m = W,

    G_0(W) = [W == 0],     G_j(W) = sum_i g_j(i) G_{j-1}(W - i S_j),

and G = G_n(V). On a net whose only invariant is "the tokens are conserved" this
is Buzen's convolution for a closed queueing network, one place per station.

NO ILP IS SOLVED. The paper obtains the marking set M_p(P',W) from the
feasibility of an integer program (Prop. 5.2) so that the sum skips the terms
that contribute nothing. Here the sum simply runs over i whose residual
W - i S_j stays non-negative and the recursion returns zero on an infeasible
residual, which gives the same value: the ILP is an optimisation of the
enumeration, not part of the definition. Memoising on (j, W) keeps the walk over
the reachable residuals rather than over all of them.

S-INVARIANT REACHABILITY IS NOT CHECKED, and cannot be cheaply: no algorithm is
known that decides it without generating the reachability set (FGCS, Sec. 5.1).
On a net that fails it, {m : S m = V} is strictly larger than the reachable set
and this returns a normalising constant over unreachable markings too, which is
why mdd_rec -- which walks the reachable set itself -- is the general algorithm
and this one the special case. Compare the two on a new net before trusting this
one on it.

See also: mdd_rec, spn_sinvariants.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

from typing import Dict, Sequence

import numpy as np

from ..io.logging import line_error


def _rec(j, W, S, g, memo):
    if j == 0:
        return 1.0 if all(w == 0 for w in W) else 0.0
    key = (j, W)
    hit = memo.get(key)
    if hit is not None:
        return hit
    p = j - 1
    acc = 0.0
    for i in range(len(g[p])):
        rem = []
        feasible = True
        for r in range(len(W)):
            v = W[r] - i * S[r][p]
            if v < 0:
                feasible = False
                break
            rem.append(v)
        if not feasible:
            break              # S is non-negative, so larger i only gets worse
        if g[p][i] == 0.0:
            continue
        acc += g[p][i] * _rec(j - 1, tuple(rem), S, g, memo)
    memo[key] = acc
    return acc


def spn_conv(S, V=None, g: Sequence[Sequence[float]] = None) -> float:
    """The normalising constant by convolution over the invariant load vector.

    Parameters
    ----------
    S : S[i][p], the minimal-support S-invariants, one row per invariant; or the
        dict spn_sinvariants returned, in which case V is the factor sequence
    V : the load vector S m0, one entry per invariant
    g : g[p][i] is g_p(i), the product-form factor of i tokens in place level p;
        its length bounds the marking of that level
    """
    if isinstance(S, dict):                 # spn_conv(inv, g) off the invariant basis
        S, V, g = S['S'], S['V'], V
    if g is None:
        line_error('spn_conv', 'the per-level product-form factors g are required')
    if len(S) == 0:
        line_error('spn_conv', 'the net has no S-invariant to convolve over')
    if len(S) != len(V):
        line_error('spn_conv', 'one load-vector entry per invariant is required')
    n = len(g)
    for r in range(len(S)):
        if len(S[r]) != n:
            line_error('spn_conv',
                       'the invariant matrix and g must agree on the place-level count')
        for p in range(n):
            if S[r][p] < 0:
                line_error('spn_conv',
                           'an S-invariant has a negative weight, so the residual recursion has '
                           'no monotone bound on the marking')
    gl = [np.asarray(g[p], dtype=float).ravel().tolist() for p in range(n)]
    Sl = [[int(S[r][p]) for p in range(n)] for r in range(len(S))]
    memo: Dict = {}
    return float(_rec(n, tuple(int(v) for v in V), Sl, gl, memo))
