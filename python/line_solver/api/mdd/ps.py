"""
Kronecker rate descriptor for shared-server stations with phase-type service.

Under processor sharing every job at a station is in service at once, each
holding its own phase, so naming a single in-service phase (what mdd_descriptor
does, which is non-preemptive semantics) cannot represent the state. The local
state here is instead the PER-PHASE COUNT vector v = (v_1,...,v_h), v_a jobs in
phase a, with n = sum(v) jobs present. That is still a per-station quantity, so
every event stays a product of per-level terms and the Kronecker form of Eq. 1
survives.

With one server shared by n jobs each job advances at rate 1/n, so from local
state v with n = sum(v):
    internal   v -> v - e_a + e_b   at v_a * D0[a,b] / n     (a != b)
    departure  v -> v - e_a         at v_a * t[a] / n * P[i,j]
    arrival    v -> v + e_b         at pie[b]
An infinite-server (delay) station is the same without the 1/n scaling. For
h = 1 the departure rate collapses to n*mu/n = mu at PS and to n*mu at IS,
reproducing the usual single-server and delay rate laws.

The local domain is the number of compositions of 0..N over h phases,
C(N+h,h), against 1+N*h for the non-preemptive encoding: the price of tracking
every job's phase rather than one.

See also: mdd_descriptor, mdd_mcd, mdd_reachset.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

from typing import Dict, List, Optional

import numpy as np
import scipy.sparse as sp

from ..io.logging import line_error
from .descriptor import _entry_law, _service_law


def mdd_ps(mu, P, servers, N, options: Optional[Dict] = None) -> Dict:
    """Kronecker rate descriptor of a closed QN whose stations are PS or IS.

    Parameters
    ----------
    mu : length-K station service rates (ignored where options['proc'] gives a law)
    P : K x K routing matrix (row-stochastic)
    servers : length-K servers per station, 1 (PS) or inf (IS); no other value
        has a per-phase-count encoding here
    N : closed population
    options : dict with keys proc and pie as in mdd_descriptor

    Returns
    -------
    dict descriptor for mdd_mcd, with K, N, domain, mu, servers, P, nphases,
    valuemap, init, nextfun, events, and the population invariant through
    valuemap.
    """
    if options is None:
        options = {}
    proc = list(options.get('proc') or [])
    pie_opt = list(options.get('pie') or [])

    mu = np.asarray(mu, dtype=float).ravel()
    servers = np.asarray(servers, dtype=float).ravel()
    P = np.asarray(P, dtype=float)
    N = int(N)
    K = mu.size

    D0: List[np.ndarray] = [None] * K
    D1: List[np.ndarray] = [None] * K
    pie: List[np.ndarray] = [None] * K
    h = np.ones(K, dtype=int)
    for i in range(K):
        if not (servers[i] == 1 or np.isinf(servers[i])):
            line_error('mdd_ps',
                       'station %d has %g servers; only processor sharing (1) and infinite server '
                       '(inf) have a per-phase-count encoding here' % (i + 1, servers[i]))
        pr = proc[i] if i < len(proc) else None
        D0[i], D1[i], _, h[i] = _service_law(pr, mu[i], i)
        if h[i] == 1 and (pr is None or (isinstance(pr, (list, tuple)) and len(pr) == 0)):
            pie[i] = np.array([1.0])
            continue
        pie[i] = _entry_law(pr, pie_opt[i] if i < len(pie_opt) else None, int(h[i]), i,
                            D1=D1[i], caller='mdd_ps')

    # ---- per-station composition state space
    comp: List[np.ndarray] = [None] * K
    lut: List[Dict] = [None] * K
    d = np.zeros(K, dtype=int)
    for i in range(K):
        comp[i], lut[i] = _comp_states(int(h[i]), N)
        d[i] = comp[i].shape[0]

    desc = {
        'K': K,
        'N': N,
        'domain': d,
        'mu': mu,
        'servers': servers,
        'P': P,
        'nphases': h,
        'valuemap': [comp[i].sum(axis=1).astype(float) for i in range(K)],
    }

    # ---- initial state: all jobs at station 1, entered in its entry phase
    init = np.zeros(K, dtype=int)
    for i in range(K):
        v = np.zeros(int(h[i]), dtype=int)
        if i == 0:
            v[int(np.nonzero(pie[0] > 0)[0][0])] = N
        init[i] = lut[i][tuple(v.tolist())]
    desc['init'] = init
    desc['nextfun'] = lambda s: _next(s, D0, D1, pie, h, comp, lut, N, P)

    # ---- events
    events = []
    for i in range(K):
        if h[i] == 1:
            continue
        Wi = _internal(D0[i], comp[i], lut[i], int(h[i]), int(d[i]), float(servers[i]))
        if Wi.nnz == 0:
            continue
        events.append({'a': i, 'b': i, 'lev': [i], 'W': [Wi]})

    aa, bb = np.nonzero(P)
    for a, b in zip(aa.tolist(), bb.tolist()):
        if a == b:
            continue
        pr_ab = float(P[a, b])
        Wa = _departure(D1[a], comp[a], lut[a], int(h[a]), int(d[a]), float(servers[a]), pr_ab)
        Wb = _arrival(pie[b], comp[b], lut[b], int(h[b]), N, int(d[b]))
        events.append({'a': a, 'b': b, 'lev': [a, b], 'W': [Wa, Wb]})
    desc['events'] = events
    return desc


# ---------------------------------------------------------------------------
def _comp_states(h, N):
    """Every per-phase count vector with 0 <= sum <= N, plus a lookup table."""
    C = _enum(h, N)
    lut = {tuple(int(x) for x in C[r, :]): r for r in range(C.shape[0])}
    return C, lut


def _enum(h, N):
    """Compositions of 0..N over h phases, in the reference's row order."""
    if h == 1:
        return np.arange(N + 1, dtype=int).reshape(-1, 1)
    sub = _enum(h - 1, N)
    subsum = sub.sum(axis=1)
    blocks = []
    for v1 in range(N + 1):
        ok = subsum <= (N - v1)
        if not np.any(ok):
            continue
        blocks.append(np.column_stack([np.full(int(ok.sum()), v1, dtype=int), sub[ok, :]]))
    return np.vstack(blocks) if blocks else np.zeros((0, h), dtype=int)


def _share(n, srv):
    """Service-rate scaling: 1/n shared by n jobs at PS, unscaled at IS."""
    if n == 0:
        return 0.0
    if np.isinf(srv):
        return 1.0
    return 1.0 / n


# ---------------------------------------------------------------------------
def _internal(D0i, C, lut, h, d, srv):
    ri, ci, vv = [], [], []
    for r in range(d):
        v = C[r, :]
        n = int(v.sum())
        if n == 0:
            continue
        sc = _share(n, srv)
        for a in range(h):
            if v[a] == 0:
                continue
            for b in range(h):
                if a == b or D0i[a, b] == 0:
                    continue
                w = v.copy()
                w[a] -= 1
                w[b] += 1
                ri.append(r)
                ci.append(lut[tuple(int(x) for x in w)])
                vv.append(float(v[a]) * D0i[a, b] * sc)
    return sp.coo_matrix((vv, (ri, ci)), shape=(d, d)).tocsr()


def _departure(D1i, C, lut, h, d, srv, pr):
    t = np.asarray(D1i).sum(axis=1).ravel()
    ri, ci, vv = [], [], []
    for r in range(d):
        v = C[r, :]
        n = int(v.sum())
        if n == 0:
            continue
        sc = _share(n, srv)
        for a in range(h):
            if v[a] == 0 or t[a] == 0:
                continue
            w = v.copy()
            w[a] -= 1
            ri.append(r)
            ci.append(lut[tuple(int(x) for x in w)])
            vv.append(float(v[a]) * t[a] * sc * pr)
    return sp.coo_matrix((vv, (ri, ci)), shape=(d, d)).tocsr()


def _arrival(pieb, C, lut, h, N, d):
    ri, ci, vv = [], [], []
    for r in range(d):
        v = C[r, :]
        if int(v.sum()) >= N:
            continue
        for b in range(h):
            if pieb[b] == 0:
                continue
            w = v.copy()
            w[b] += 1
            ri.append(r)
            ci.append(lut[tuple(int(x) for x in w)])
            vv.append(float(pieb[b]))
    return sp.coo_matrix((vv, (ri, ci)), shape=(d, d)).tocsr()


# ---------------------------------------------------------------------------
def _next(s, D0, D1, pie, h, comp, lut, N, P):
    K = len(s)
    T = []
    for i in range(K):
        v = comp[i][int(s[i]), :]
        n = int(v.sum())
        if n == 0:
            continue
        # internal phase moves
        for a in range(int(h[i])):
            if v[a] == 0:
                continue
            for b in range(int(h[i])):
                if a == b or D0[i][a, b] == 0:
                    continue
                w = v.copy()
                w[a] -= 1
                w[b] += 1
                t = list(s)
                t[i] = lut[i][tuple(int(x) for x in w)]
                T.append(t)
        # completions routed to j
        ti = np.asarray(D1[i]).sum(axis=1).ravel()
        for a in range(int(h[i])):
            if v[a] == 0 or ti[a] == 0:
                continue
            w = v.copy()
            w[a] -= 1
            for j in range(K):
                if j == i or P[i, j] <= 0:
                    continue
                vj = comp[j][int(s[j]), :]
                if int(vj.sum()) >= N:
                    continue
                for b in range(int(h[j])):
                    if pie[j][b] == 0:
                        continue
                    wj = vj.copy()
                    wj[b] += 1
                    t = list(s)
                    t[i] = lut[i][tuple(int(x) for x in w)]
                    t[j] = lut[j][tuple(int(x) for x in wj)]
                    T.append(t)
    return T
