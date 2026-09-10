"""
Position-resolved mean-field miss rates for FIFO(m) caches.

FIFO(m) and RANDOM(m) share the exact stationary distribution (Gast and Van
Houdt, SIGMETRICS 2015, Thm 1: pi_FIFO(m) = pi_RAND(m)), so their steady-state
hit ratios coincide and SolverFLD serves FIFO steady state from the cheaper
RAND(m) refined mean field (cache_miss_rmf). Their mean-field TRANSIENTS differ,
however: FIFO evicts the deterministic tail (fixed residence of m insertions)
whereas RANDOM evicts a uniformly random victim (geometric residence), so the
hit-probability trajectory H(t) from a cold cache ramps differently even though
H(inf) agrees. This module provides that dedicated FIFO transient via a
position-resolved density-dependent population process (DDPP).

FIFO(m) differs from strict FIFO(m) only in the reinsertion position on a hit:
the demoted tail of list i+1 lands at the vacated position j of list i (in
place, no within-list shift), whereas strict FIFO reinserts it at position 1.

State x[k,i,j] = P(item k in position j of list i). Dynamics (aggregate IRM):
  - miss: insert missed item at position 1 of list 1; list 1 shifts back; tail
    (position m_1) is evicted.
  - hit at position j of list i < h: promote that item to position 1 of list
    i+1 (list i+1 shifts back, tail demoted); the demoted tail replaces the
    vacated position j of list i (no other position of list i moves).
  - hit in the top list h: no change.
"""

import numpy as np
from scipy.integrate import solve_ivp

__all__ = ['cache_miss_fifo_rmf']


def _build_slotmap(m, h):
    slots = [(i, j) for i in range(1, h + 1) for j in range(1, int(m[i - 1]) + 1)]
    sidx = {s: t for t, s in enumerate(slots)}
    return slots, sidx, len(slots)


def _fifo_drift(x, p, m, n, h, slots, sidx, S):
    """Mean-field drift F(x) for FIFO(m); x is length n*S (in-cache slots)."""
    x = np.clip(x, 0.0, 1.0)

    def K(k, i, j):
        return k * S + sidx[(i, j)]

    def occ_list(k, i):
        return sum(x[K(k, i, jj)] for jj in range(1, int(m[i - 1]) + 1))

    def out(k):
        return 1.0 - sum(x[K(k, i, jj)] for (i, jj) in slots)

    Hpos = {}
    Hi = [0.0] * (h + 1)
    for (i, j) in slots:
        s = sum(p[k] * x[K(k, i, j)] for k in range(n))
        Hpos[(i, j)] = s
        Hi[i] += s
    M = sum(p[k] * out(k) for k in range(n))

    Sfull = [0.0] * (h + 1)
    Sfull[1] = M
    for i in range(2, h + 1):
        Sfull[i] = Hi[i - 1]

    d = np.zeros(n * S)
    for k in range(n):
        for (i, j) in slots:
            xk = x[K(k, i, j)]
            # outflow: full shift toward j+1 (tail leaves); promotion up if i<h
            o = Sfull[i] * xk
            if i < h:
                o += p[k] * xk
            d[K(k, i, j)] -= o
            # inflow: full shift from j-1, or front insertion at j == 1
            if j >= 2:
                d[K(k, i, j)] += Sfull[i] * x[K(k, i, j - 1)]
            else:
                if i == 1:
                    d[K(k, 1, 1)] += p[k] * out(k)             # miss inserts item k
                else:
                    d[K(k, i, 1)] += p[k] * occ_list(k, i - 1)  # promotion from i-1
            # FIFO demotion: tail of i+1 lands in place at the same position j
            if i < h:
                d[K(k, i, j)] += Hpos[(i, j)] * x[K(k, i + 1, int(m[i]))]
    return d


def _initial_state(p, m, n, h, slots, sidx, S):
    x0 = np.zeros(n * S)
    order = np.argsort(-p)
    pos = 0
    for (i, j) in slots:
        if pos < n:
            x0[order[pos] * S + sidx[(i, j)]] = 1.0
            pos += 1
    return x0


def cache_miss_fifo_rmf(gamma, m, lambd, tspan=None, x0init=None, tmax=20000.0, accost=None):
    """
    Position-resolved mean-field miss rates for FIFO(m) caches.

    Mirrors the cache_miss_rmf / cache_miss_sfifo_rmf contract. gamma is
    accepted for interface parity (used for sizing only); the popularity is
    recovered from lambd, the (u, n, h+1) per-user per-item arrival rates.

    Returns (M, MU, MI, pi0) when tspan is None, else
    (M, MU, MI, pi0, tout, pi0_t, MU_t, xtraj). The transient path is the
    intended use: at steady state FIFO(m) equals RANDOM(m) (Gast15 Thm 1), so
    prefer cache_miss_rmf when only the fixed point is needed.
    """
    lambd = np.asarray(lambd, dtype=float)
    u, n = lambd.shape[0], lambd.shape[1]
    m = np.asarray(m, dtype=float).ravel()
    h = len(m)
    slots, sidx, S = _build_slotmap(m, h)

    lam_i = np.zeros(n)
    for v in range(u):
        row = np.array(lambd[v, :, 0], dtype=float)
        row[~np.isfinite(row)] = 0.0
        lam_i += row
    tot = np.sum(lam_i)
    p = lam_i / tot if tot > 0 else np.full(n, 1.0 / n)

    def out_of(xss, k):
        return max(0.0, min(1.0, 1.0 - sum(xss[k * S + sidx[s]] for s in slots)))

    # Non-default access graph: general position-resolved drift from a cold
    # (empty) cache so non-admissible items drain; linear default keeps the
    # pre-filled path.
    from .rmf import _build_item_graphs
    from .rmf_sfifo import _pos_drift_graph
    _G = _build_item_graphs(accost, lambd, n, h)
    if _G is not None:
        _Garr = [np.asarray(g, dtype=float) for g in _G]
        drift_fn = lambda t, x: _pos_drift_graph(x, p, _Garr, m, n, h, slots, sidx, S, 'pos')
        x0 = np.zeros(n * S)  # cold start so non-admissible items drain
    else:
        drift_fn = lambda t, x: _fifo_drift(x, p, m, n, h, slots, sidx, S)
        x0 = _initial_state(p, m, n, h, slots, sidx, S)
    sol = solve_ivp(drift_fn, (0.0, tmax), x0, method='LSODA', rtol=1e-8, atol=1e-10)
    xss = np.clip(sol.y[:, -1], 0.0, 1.0)

    pi0 = np.array([out_of(xss, k) for k in range(n)])
    MI = lam_i * pi0
    MU = np.array([np.nansum(np.array(lambd[v, :, 0], dtype=float) * pi0)
                   for v in range(u)])
    M = float(np.sum(MI))

    if tspan is None:
        return M, MU, MI, pi0

    x0t = x0 if x0init is None else np.asarray(x0init, dtype=float).ravel()
    solt = solve_ivp(drift_fn, (float(tspan[0]), float(tspan[-1])), x0t, method='LSODA',
                     rtol=1e-8, atol=1e-10)
    tout = solt.t
    xtraj = np.clip(solt.y, 0.0, 1.0)
    nt = len(tout)
    pi0_t = np.zeros((n, nt))
    for k in range(n):
        for c in range(nt):
            pi0_t[k, c] = out_of(xtraj[:, c], k)
    MU_t = np.zeros((u, nt))
    for v in range(u):
        row = np.array(lambd[v, :, 0], dtype=float)
        row[~np.isfinite(row)] = 0.0
        MU_t[v, :] = row @ pi0_t
    return M, MU, MI, pi0, tout, pi0_t, MU_t, xtraj
