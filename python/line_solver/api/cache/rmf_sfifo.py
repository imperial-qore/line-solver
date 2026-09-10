"""
Position-resolved mean-field miss rates for strict FIFO(m) caches (SFIFO).

Strict FIFO(m) is NOT equivalent to RANDOM(m)/FIFO(m): while Gast and Van Houdt
(SIGMETRICS 2015, "Transient and Steady-state Regime of a Family of List-based
Cache Replacement Algorithms") prove pi_FIFO(m) = pi_RAND(m) exactly, they show
strict FIFO(m) differs (their Section 3) and provide no mean-field model for it
(it is only trace-simulated in their Section 6.2). The difference lives in the
within-list age ordering: on a hit in list i < h the demoted tail of list i+1 is
reinserted at position 1 of list i (pushing 1..j-1 back), which the per-item
per-list occupancy state of RANDOM(m) cannot represent.

This module closes that gap with a position-resolved density-dependent
population process (DDPP): the state is x[k,i,j] = P(item k occupies position j
of list i), with deterministic (age-based) demotion and eviction. The plain
mean-field fixed point (obtained by integrating the drift to steady state) is
returned; it reduces exactly to the RANDOM(m)/FIFO(m) result when
m_1 = ... = m_{h-1} = 1, matching the strict-FIFO == FIFO degeneracy.

Strict FIFO(m) dynamics (single aggregate request stream, IRM):
  - miss: insert the missed item at position 1 of list 1; list 1 shifts back one
    position; the tail (position m_1) is evicted from the cache.
  - hit at position j of list i < h: promote that item to position 1 of list
    i+1 (list i+1 shifts back, its tail demoted); the demoted tail enters
    position 1 of list i and positions 1..j-1 of list i shift back one.
  - hit in the top list h: no change.
"""

import numpy as np
from scipy.integrate import solve_ivp

__all__ = ['cache_miss_sfifo_rmf']


def _build_slotmap(m, h):
    """Return (slots, sidx, S) for lists 1..h with sizes m[0..h-1]."""
    slots = [(i, j) for i in range(1, h + 1) for j in range(1, int(m[i - 1]) + 1)]
    sidx = {s: t for t, s in enumerate(slots)}
    return slots, sidx, len(slots)


def _sfifo_drift(x, p, m, n, h, slots, sidx, S):
    """Mean-field drift F(x) for strict FIFO(m); x is length n*S (in-cache slots)."""
    x = np.clip(x, 0.0, 1.0)

    def K(k, i, j):
        return k * S + sidx[(i, j)]

    def occ_list(k, i):
        return sum(x[K(k, i, jj)] for jj in range(1, int(m[i - 1]) + 1))

    def out(k):
        return 1.0 - sum(x[K(k, i, jj)] for (i, jj) in slots)

    # per-position and per-list aggregate hit rates, and the miss rate
    Hpos = {}
    Hi = [0.0] * (h + 1)
    for (i, j) in slots:
        s = sum(p[k] * x[K(k, i, j)] for k in range(n))
        Hpos[(i, j)] = s
        Hi[i] += s
    M = sum(p[k] * out(k) for k in range(n))

    # full-shift rate of each list (list 1 on a miss; list i on a hit in i-1)
    Sfull = [0.0] * (h + 1)
    Sfull[1] = M
    for i in range(2, h + 1):
        Sfull[i] = Hi[i - 1]

    def Gi(i, jp):
        # partial-shift rate of a slot at position jp of list i: hits at deeper
        # positions of list i. The top list h never moves on a hit.
        if i == h:
            return 0.0
        return sum(Hpos[(i, jj)] for jj in range(jp + 1, int(m[i - 1]) + 1))

    d = np.zeros(n * S)
    for k in range(n):
        for (i, j) in slots:
            mi = int(m[i - 1])
            xk = x[K(k, i, j)]
            # outflow: shift toward j+1 (or leave the list at the tail), and
            # promotion up to list i+1 when this item is requested (i < h)
            o = (Sfull[i] + Gi(i, j)) * xk
            if i < h:
                o += p[k] * xk
            d[K(k, i, j)] -= o
            # inflow
            if j >= 2:
                d[K(k, i, j)] += (Sfull[i] + Gi(i, j - 1)) * x[K(k, i, j - 1)]
            else:
                if i == 1:
                    d[K(k, 1, 1)] += p[k] * out(k)          # miss inserts item k
                else:
                    d[K(k, i, 1)] += p[k] * occ_list(k, i - 1)  # promotion from i-1
                if i < h:
                    d[K(k, i, 1)] += Hi[i] * x[K(k, i + 1, int(m[i]))]  # demotion from i+1 tail
    return d


def _initial_state(p, m, n, h, slots, sidx, S):
    """Popularity-ordered warm start: most popular items fill the cache slots."""
    x0 = np.zeros(n * S)
    order = np.argsort(-p)
    pos = 0
    for (i, j) in slots:
        if pos < n:
            x0[order[pos] * S + sidx[(i, j)]] = 1.0
            pos += 1
    return x0


def _pos_drift_graph(x, p, G, m, n, h, slots, sidx, S, reinsert):
    """General position-resolved mean-field drift honouring a per-item access
    graph, for FIFO(m) (reinsert='pos') and strict FIFO(m) (reinsert='head').

    G[k] is the (h+1)x(h+1) access graph of item k (row 0 = miss admission per
    list, row 1+i = hit-in-list-i promotion target b>=i; b==i means STAY in
    place, the FIFO/SFIFO convention). A miss admits at the head of the target
    list (its tail evicted); a hit at position j of list i promotes to the head
    of target b>i (the tail of b demoted to list i -- to the vacated position j
    for FIFO, to the head with a 1..j-1 shift for SFIFO). Reduces exactly to the
    linear drift when G is the standard chain. Requires a cold (empty) initial
    state so non-admissible items drain.
    """
    x = np.clip(x, 0.0, 1.0)

    def K(k, i, j):
        return k * S + sidx[(i, j)]

    def occ(k, i):
        return sum(x[K(k, i, jj)] for jj in range(1, int(m[i - 1]) + 1))

    def out(k):
        return 1.0 - sum(x[K(k, i, jj)] for (i, jj) in slots)

    MI = np.zeros(h + 1)
    HP = np.zeros((h + 1, h + 1))
    for k in range(n):
        ok = out(k)
        gk = G[k]
        for l in range(1, h + 1):
            MI[l] += p[k] * ok * gk[0, l]
        for i in range(1, h + 1):
            oc = occ(k, i)
            for b in range(i + 1, h + 1):
                HP[i, b] += p[k] * oc * gk[i, b]
    Sin = np.zeros(h + 1)
    for l in range(1, h + 1):
        Sin[l] = MI[l] + sum(HP[s, l] for s in range(1, l))
    POp = {}
    for (i, j) in slots:
        POp[(i, j)] = sum(p[k] * x[K(k, i, j)] * (1.0 - G[k][i, i]) for k in range(n))

    def Gg(i, jp):
        return sum(POp[(i, jj)] for jj in range(jp + 1, int(m[i - 1]) + 1))

    d = np.zeros(n * S)
    for k in range(n):
        gk = G[k]
        ok = out(k)
        for (i, j) in slots:
            xk = x[K(k, i, j)]
            # outflow: active promote-out + shift
            o = p[k] * xk * (1.0 - gk[i, i])
            if reinsert == 'head':
                o += (Sin[i] + Gg(i, j)) * xk
            else:
                o += Sin[i] * xk
            d[K(k, i, j)] -= o
            # inflow via shift from j-1, or head insertion at j == 1
            if j >= 2:
                if reinsert == 'head':
                    d[K(k, i, j)] += (Sin[i] + Gg(i, j - 1)) * x[K(k, i, j - 1)]
                else:
                    d[K(k, i, j)] += Sin[i] * x[K(k, i, j - 1)]
            else:
                d[K(k, i, 1)] += p[k] * ok * gk[0, i]                  # miss admit to i
                for s in range(1, i):
                    d[K(k, i, 1)] += p[k] * occ(k, s) * gk[s, i]        # promotion-in from s<i
            # demotion inflow: tail of b -> list i (head for SFIFO, vacated pos for FIFO)
            for b in range(i + 1, h + 1):
                if reinsert == 'head':
                    if j == 1:
                        d[K(k, i, 1)] += HP[i, b] * x[K(k, b, int(m[b - 1]))]
                else:
                    poj = sum(p[kk] * x[K(kk, i, j)] * G[kk][i, b] for kk in range(n))
                    d[K(k, i, j)] += poj * x[K(k, b, int(m[b - 1]))]
    return d


def cache_miss_sfifo_rmf(gamma, m, lambd, tspan=None, x0init=None, tmax=20000.0, accost=None):
    """
    Position-resolved mean-field miss rates for strict FIFO(m) caches.

    Mirrors the cache_miss_rmf / cache_miss_fpi contract. gamma is accepted for
    interface compatibility (used for sizing only); the popularity is recovered
    from lambd, the (u, n, h+1) per-user per-item arrival rates.

    Returns (M, MU, MI, pi0) when tspan is None, else
    (M, MU, MI, pi0, tout, pi0_t, MU_t, xtraj) with tout (nt,), pi0_t (n, nt)
    per-item out-of-cache occupancy, MU_t (u, nt) per-user miss rate, and
    xtraj (n*S, nt) the full position-resolved occupancy trajectory.
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

    # A non-default access graph (accost) modulates admission/promotion per
    # item. Use the general position-resolved drift from a COLD (empty) cache so
    # non-admissible items drain; the linear default keeps the pre-filled path.
    from .rmf import _build_item_graphs
    G = _build_item_graphs(accost, lambd, n, h)
    if G is not None:
        Garr = [np.asarray(g, dtype=float) for g in G]
        drift_fn = lambda t, x: _pos_drift_graph(x, p, Garr, m, n, h, slots, sidx, S, 'head')
        x0 = np.zeros(n * S)  # cold start so non-admissible items drain
    else:
        drift_fn = lambda t, x: _sfifo_drift(x, p, m, n, h, slots, sidx, S)
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
                     rtol=1e-8, atol=1e-10, dense_output=False)
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
