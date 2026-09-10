"""
Minimal-support S-invariants (P-invariants) of a stochastic Petri net, and the
load vector V = S m0.

An S-invariant is a non-negative left null vector of the incidence matrix,
U' C = 0, so U' m is conserved by every firing. The minimal-support ones form a
basis of all of them (S. Balsamo, A. Marin, I. Stojic, FGCS 111 (2020), Sec. 3.1)
and are what the convolution algorithm spn_conv decomposes the reachability set
along; spn_mdd uses a single positive invariant for a much weaker purpose, to
bound each place a priori.

FARKAS' ALGORITHM, on [C | I]: for each transition column in turn, keep the rows
that already annihilate it and add, for every pair of rows of opposite sign in
it, the positive combination that cancels it; then drop every row whose support
strictly contains another's, which is what leaves the minimal supports. Rows are
kept in integer arithmetic and divided by their gcd, so a multiplicity is never
lost to rounding and two invariants that differ only by a positive scale are the
same row.

ARC MULTIPLICITIES MUST BE INTEGRAL. A fractional arc has no Petri-net meaning
and would make the gcd normalisation and the ILP-free convolution both wrong, so
it is refused rather than rounded.

Levels are the (place, class) pairs of spn_mdd, place-major, so the invariants
come out in the coordinates the decision diagram and spn_conv both use.

See also: spn_conv, spn_mdd.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

from math import gcd
from typing import Dict, List, Optional, Sequence

import numpy as np

from ..io.logging import line_error
from ..sn.network_struct import NodeType


def _as_integer(x, what):
    """An arc multiplicity, refused unless integral."""
    r = float(np.floor(x + 0.5))
    if abs(x - r) > 1e-9:
        line_error('spn_sinvariants',
                   '%s is not integral; a fractional arc multiplicity has no Petri-net meaning '
                   'and no invariant basis over the integers' % what)
    return int(r)


def _param(nparam, name, default=None):
    if nparam is None:
        return default
    if isinstance(nparam, dict):
        return nparam.get(name, default)
    return getattr(nparam, name, default)


def _support_subset(a, b, off, n):
    """True when the support of a is contained in the support of b."""
    for k in range(n):
        if a[off + k] != 0 and b[off + k] == 0:
            return False
    return True


def _support_equal(a, b, off, n):
    for k in range(n):
        if (a[off + k] != 0) != (b[off + k] != 0):
            return False
    return True


def spn_sinvariants(sn, init: Optional[Sequence[float]] = None) -> Dict[str, object]:
    """Minimal-support S-invariants and the load vector of a net.

    Parameters
    ----------
    sn : a NetworkStruct holding Places and Transitions
    init : initial tokens per place level, place-major; None takes them from the
        reference station of each closed class, as spn_mdd does

    Returns
    -------
    dict with 'places' (0-based node indices), 'S' (one row per invariant, one
    column per place level), 'V' = S m0 and 'm0'.
    """
    nodetype = np.ravel(np.asarray(sn.nodetype, dtype=int))
    places = [int(i) for i in np.nonzero(nodetype == int(NodeType.PLACE))[0]]
    transitions = [int(i) for i in np.nonzero(nodetype == int(NodeType.TRANSITION))[0]]
    if not places or not transitions:
        line_error('spn_sinvariants', 'the model holds no Place or no Transition node')
    R = int(sn.nclasses)
    nnodes = int(sn.nnodes)
    P = len(places)
    n = P * R

    # ---- incidence matrix C[level][mode] = post - pre, one column per (transition, mode)
    cols: List[List[int]] = []
    for ind in transitions:
        nparam = sn.nodeparam[ind] if sn.nodeparam is not None else None
        nmodes = int(_param(nparam, 'nmodes', 0) or 0)
        enabling = _param(nparam, 'enabling', None) or []
        firing = _param(nparam, 'firing', None) or []
        for m in range(nmodes):
            pre = np.zeros((nnodes, R)) if m >= len(enabling) or enabling[m] is None \
                else np.asarray(enabling[m], dtype=float).reshape(nnodes, R)
            post = np.zeros((nnodes, R)) if m >= len(firing) or firing[m] is None \
                else np.asarray(firing[m], dtype=float).reshape(nnodes, R)
            col = []
            for pp in range(P):
                for k in range(R):
                    col.append(_as_integer(max(0.0, post[places[pp], k]), 'a firing arc') -
                               _as_integer(max(0.0, pre[places[pp], k]), 'an enabling arc'))
            cols.append(col)
    ncols = len(cols)

    # ---- Farkas on [C | I]: row p starts as (C[p], e_p)
    rows: List[List[int]] = []
    for p in range(n):
        row = [0] * (ncols + n)
        for c in range(ncols):
            row[c] = cols[c][p]
        row[ncols + p] = 1
        rows.append(row)

    for c in range(ncols):
        nxt = [r for r in rows if r[c] == 0]
        for a in rows:
            if a[c] <= 0:
                continue
            for b in rows:
                if b[c] >= 0:
                    continue
                pa, nb = a[c], -b[c]
                d = gcd(pa, nb)
                fa, fb = nb // d, pa // d
                combo = [fa * a[k] + fb * b[k] for k in range(ncols + n)]
                gall = 0
                for x in combo:
                    gall = gcd(gall, abs(x))
                if gall > 1:
                    combo = [x // gall for x in combo]
                if any(combo[ncols + k] != 0 for k in range(n)):
                    nxt.append(combo)
        # support-minimality filter, applied at every step so the row set cannot
        # grow combinatorially on the way to the answer
        keep = []
        for r in range(len(nxt)):
            dominated = False
            for s in range(len(nxt)):
                if s == r or dominated:
                    continue
                if not _support_subset(nxt[s], nxt[r], ncols, n):
                    continue
                same = _support_equal(nxt[s], nxt[r], ncols, n)
                if not same or s < r:      # keep the first of equal supports
                    dominated = True
            if not dominated:
                keep.append(nxt[r])
        rows = keep

    S = []
    for r in rows:
        y = [r[ncols + k] for k in range(n)]
        if any(v < 0 for v in y):
            continue                        # an S-invariant is non-negative by definition
        S.append(y)

    # ---- initial marking and the load vector V = S m0
    m0 = [0] * n
    if init is not None:
        iv = np.ravel(np.asarray(init, dtype=float))
        if iv.size != n:
            line_error('spn_sinvariants', 'init must hold one token count per place level')
        m0 = [_as_integer(float(iv[i]), 'an initial marking') for i in range(n)]
    else:
        njobs = np.ravel(np.asarray(sn.njobs, dtype=float))
        refstat = np.ravel(np.asarray(sn.refstat, dtype=int))
        node_to_station = np.ravel(np.asarray(sn.nodeToStation, dtype=int))
        for k in range(R):
            if not np.isfinite(njobs[k]):
                line_error('spn_sinvariants',
                           'class %d is open, so the net has no finite load vector' % (k + 1))
            ref = int(refstat[k]) if k < refstat.size else 0
            for pp in range(P):
                if node_to_station[places[pp]] == ref:
                    m0[pp * R + k] += _as_integer(float(njobs[k]), 'a class population')
    V = [int(sum(S[i][p] * m0[p] for p in range(n))) for i in range(len(S))]
    return {'places': places, 'S': S, 'V': V, 'm0': m0}
