"""
Refined mean field (RMF) miss rates for RANDOM(m) caches.

Computes global, per-user, and per-item miss rates for multi-list caches
with RANDOM(m) replacement using the DDPP mean field approximation with
1/N correction (refined mean field). The cache occupancy is computed from
the aggregate request stream; per-user miss rates weight the per-item miss
probabilities by each user's rates. Mirrors the cache_miss_fpi contract.

Reference:
    N. Gast, "Expected Values Estimated via Mean-Field Approximation are
    1/N-Accurate", Proc. ACM Meas. Anal. Comput. Syst., 2017.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_lyapunov, svd

__all__ = ['cache_miss_rmf']


def _idx(i, k, n):
    return i + k * n


def _hit_rate(x, p, k, n):
    return float(sum(p[i] * x[_idx(i, k, n)] for i in range(n)))


def _drift(x, p, m, n, h, dim):
    hit = [_hit_rate(x, p, k, n) for k in range(h + 1)]
    dX = np.zeros(dim)
    for i in range(n):
        for k in range(h):
            flow = p[i] * x[_idx(i, k, n)] - hit[k] * x[_idx(i, k + 1, n)] / m[k]
            dX[_idx(i, k, n)] -= flow
            dX[_idx(i, k + 1, n)] += flow
    return dX


def _drift_graph(x, p, G, m, n, h, dim):
    """General RANDOM(m) mean-field drift honouring a per-item access graph.

    G[k] is the (h+1)x(h+1) access/promotion probability matrix of item k:
    row 0 is miss admission (col 0 = reject/no-cache, col l = admit to list l);
    row a>=1 is a hit in list a (col b = promote to list b, b>=a). On a miss into
    list l or a promotion into list l a uniformly random occupant of list l is
    displaced (evicted out on a miss, swapped to the source list on a hit),
    matching the exact RR sample-path (State.afterEventCache). Reduces to the
    linear eq-8 drift when G is the standard chain (miss->1, hit a->a+1).
    """
    # aggregate insertion/promotion rate into list i from source s (s=0 miss)
    A = np.zeros((h + 1, h + 1))  # A[s, i]
    for s in range(h + 1):
        for j in range(n):
            xjs = x[_idx(j, s, n)]
            if xjs == 0.0:
                continue
            gj = G[j]
            for i in range(1, h + 1):
                A[s, i] += p[j] * xjs * gj[s, i]

    dX = np.zeros(dim)
    for k in range(n):
        outk = x[_idx(k, 0, n)]
        gk = G[k]
        for i in range(1, h + 1):
            xki = x[_idx(k, i, n)]
            # active inflow: miss admission + promotion into i from sources s<i
            infl = p[k] * outk * gk[0, i]
            for s in range(1, i):
                infl += p[k] * x[_idx(k, s, n)] * gk[s, i]
            # passive victim inflow: hit at i -> b>i displaces a victim of b to i
            for b in range(i + 1, h + 1):
                infl += A[i, b] * x[_idx(k, b, n)] / m[b - 1]
            # active outflow: k promoted out of i to some b>i
            outfl = p[k] * xki * (1.0 - gk[i, i])
            # passive victim outflow: insertions into i (miss s=0 + promo s<i)
            disp = 0.0
            for s in range(0, i):
                disp += A[s, i]
            outfl += disp * xki / m[i - 1]
            dX[_idx(k, i, n)] += infl - outfl
        # out-of-cache by conservation
        dX[_idx(k, 0, n)] = -sum(dX[_idx(k, i, n)] for i in range(1, h + 1))
    return dX


def _build_item_graphs(accost, lambd, n, h):
    """Per-item (h+1)x(h+1) access graph aggregated over users by request rate.

    Returns None when accost is absent or is the standard linear chain (so the
    caller keeps the refined linear path). accost is [user][item] matrices.
    """
    if accost is None:
        return None
    G = []
    is_linear = True
    for k in range(n):
        num = np.zeros((h + 1, h + 1))
        den = 0.0
        for v in range(len(accost)):
            wv = float(np.nansum(lambd[v, k, 0])) if v < lambd.shape[0] else 0.0
            gvk = accost[v][k] if (accost[v] is not None and k < len(accost[v])) else None
            if gvk is None:
                continue
            gvk = np.asarray(gvk, dtype=float)
            num += wv * gvk
            den += wv
        if den > 0:
            gk = num / den
        else:
            gk = np.asarray(accost[0][k], dtype=float) if accost[0] is not None else _linear_graph(h)
        # normalize rows defensively
        for a in range(h + 1):
            s = gk[a, :].sum()
            if s > 0:
                gk[a, :] = gk[a, :] / s
        G.append(gk)
        if not np.allclose(gk, _linear_graph(h), atol=1e-9):
            is_linear = False
    return None if is_linear else G


def _linear_graph(h):
    g = np.zeros((h + 1, h + 1))
    g[0, 1] = 1.0
    for a in range(1, h):
        g[a, a + 1] = 1.0
    g[h, h] = 1.0
    return g


def _fixed_point_graph(x0, p, G, m, n, h, dim, tmax=20000.0):
    sol = solve_ivp(lambda t, x: _drift_graph(x, p, G, m, n, h, dim), (0.0, tmax), x0,
                    method='LSODA', rtol=1e-8, atol=1e-10)
    return sol.y[:, -1]


def _jacobian(x, p, m, n, h, dim):
    hit = [_hit_rate(x, p, k, n) for k in range(h + 1)]
    Fp = np.zeros((dim, dim))
    for i in range(n):
        for k in range(h):
            ik = _idx(i, k, n)
            ik1 = _idx(i, k + 1, n)
            Fp[ik, ik] -= p[i]
            Fp[ik1, ik] += p[i]
            Fp[ik, ik1] += hit[k] / m[k]
            Fp[ik1, ik1] -= hit[k] / m[k]
            for j in range(n):
                jk = _idx(j, k, n)
                jk1 = _idx(j, k + 1, n)
                Fp[ik, jk1] -= p[i] * x[ik] / m[k]
                Fp[ik1, jk1] += p[i] * x[ik] / m[k]
                Fp[ik, jk] += p[j] * x[ik1] / m[k]
                Fp[ik1, jk] -= p[j] * x[ik1] / m[k]
    return Fp


def _hessian(p, m, n, h, dim):
    Fpp = np.zeros((dim, dim, dim))
    for i in range(n):
        for k in range(h):
            ik = _idx(i, k, n)
            ik1 = _idx(i, k + 1, n)
            for j in range(n):
                if j != i:
                    jk = _idx(j, k, n)
                    jk1 = _idx(j, k + 1, n)
                    Fpp[ik, jk, ik1] += p[j] / m[k]
                    Fpp[ik, ik1, jk] += p[j] / m[k]
                    Fpp[ik, jk1, ik] -= p[i] / m[k]
                    Fpp[ik, ik, jk1] -= p[i] / m[k]
                    Fpp[ik1, jk, ik1] -= p[j] / m[k]
                    Fpp[ik1, ik1, jk] -= p[j] / m[k]
                    Fpp[ik1, jk1, ik] += p[i] / m[k]
                    Fpp[ik1, ik, jk1] += p[i] / m[k]
    return Fpp


def _noise_matrix(x, p, m, n, h, dim):
    Q = np.zeros((dim, dim))
    signs = np.array([-1.0, 1.0, 1.0, -1.0])
    for i in range(n):
        for k in range(h):
            for j in range(n):
                rate = p[i] * x[_idx(i, k, n)] * x[_idx(j, k + 1, n)] / m[k]
                indices = [_idx(i, k, n), _idx(j, k, n), _idx(i, k + 1, n), _idx(j, k + 1, n)]
                for ia in range(4):
                    for ib in range(4):
                        Q[indices[ia], indices[ib]] += rate * signs[ia] * signs[ib]
    return Q


def _fixed_point(x0, p, m, n, h, dim, tmax=10000.0):
    sol = solve_ivp(lambda t, x: _drift(x, p, m, n, h, dim), (0.0, tmax), x0,
                    method='LSODA', rtol=1e-8, atol=1e-10)
    return sol.y[:, -1]


def _dimension_reduction(Fp, n, h, dim):
    rk = int(np.linalg.matrix_rank(Fp))
    C = np.zeros((dim, dim))
    d = 0
    for l_idx in range(h + 1):
        for i in range(n - 1):
            C[d, _idx(i, l_idx, n)] = 1.0
            d += 1
    U, _, _ = svd(Fp)
    C[rk:dim, :] = U[:, rk:dim].T
    Cinv = np.linalg.inv(C)
    return C, Cinv, rk


def _expansion_steady_state(x0, p, m, n, h, dim):
    pi = _fixed_point(x0, p, m, n, h, dim)
    Fp = _jacobian(pi, p, m, n, h, dim)
    Fpp = _hessian(p, m, n, h, dim)
    Q = _noise_matrix(pi, p, m, n, h, dim)

    C, Cinv, rk = _dimension_reduction(Fp, n, h, dim)
    Fp_r = (C @ Fp @ Cinv)[:rk, :rk]
    # Reduce Hessian: Fpp_r(a,b,c) = sum_{i,j,k} C(a,i) Fpp(i,j,k) Cinv(j,b) Cinv(k,c)
    Fpp_r = np.einsum('ai,ijk,jb,kc->abc', C[:rk, :], Fpp, Cinv[:, :rk], Cinv[:, :rk])
    Q_r = (C @ Q @ C.T)[:rk, :rk]

    # Solve Lyapunov equation: Fp_r W + W Fp_r' + Q_r = 0
    W_r = solve_lyapunov(Fp_r, -Q_r)

    C_r = np.einsum('abc,bc->a', Fpp_r, W_r)
    V_r = -np.linalg.solve(Fp_r, C_r / 2.0)
    V = Cinv[:, :rk] @ V_r
    return pi, V


def cache_miss_rmf(gamma, m, lambd, tspan=None, x0init=None, accost=None):
    """
    RMF (1/N-accurate) miss rates for RANDOM(m) caches.

    Mirrors the cache_miss_fpi contract. gamma is accepted for interface
    compatibility (used for sizing only); the popularity is recovered from
    lambd, the (u, n, h+1) per-user per-item arrival rates.

    Args:
        gamma: item access factors (unused beyond sizing).
        m: cache capacity vector (h,).
        lambd: arrival rates per user per item per list (u, n, >=1).
        tspan: optional [t0, t1]. When given, also integrates the plain
            mean-field drift over the window and returns the transient
            trajectory (the same drift that _fixed_point drives to steady
            state). Omitting tspan preserves the steady-state-only contract.
        x0init: optional initial occupancy (dim,) for the transient; defaults
            to the standard first-m-in-list initial state. Used to carry the
            cache mean occupancy across environment switches.

    Returns:
        Tuple (M, MU, MI, pi0) when tspan is None, else
        (M, MU, MI, pi0, tout, pi0_t, MU_t, xtraj) with tout (nt,),
        pi0_t (n, nt) per-item list-0 occupancy, MU_t (u, nt) per-user miss
        rate, and xtraj (dim, nt) full DDPP occupancy trajectory.
    """
    lambd = np.asarray(lambd, dtype=float)
    u, n = lambd.shape[0], lambd.shape[1]
    m = np.asarray(m, dtype=float).ravel()
    h = len(m)
    dim = n * (h + 1)

    # aggregate per-item request rates over users
    lam_i = np.zeros(n)
    for v in range(u):
        row = np.array(lambd[v, :, 0], dtype=float)
        row[~np.isfinite(row)] = 0.0
        lam_i += row
    p = lam_i / np.sum(lam_i)

    # initial state: first m[0] items in list 1, next m[1] in list 2, etc.;
    # remaining items outside the cache (list 0)
    x0 = np.zeros(dim)
    obj_idx = 0
    for k in range(1, h + 1):
        for _ in range(int(m[k - 1])):
            if obj_idx < n:
                x0[_idx(obj_idx, k, n)] = 1.0
                obj_idx += 1
    for i in range(obj_idx, n):
        x0[_idx(i, 0, n)] = 1.0

    # mean field fixed point, then 1/N refinement when computable. The
    # refinement solves a reduced Jacobian/Lyapunov system that turns singular
    # for small or non-hyperbolic fixed points; numpy returns non-finite
    # entries (inf/nan) rather than raising, so a bare except does not catch
    # it. Accept the 1/N correction only when finite, otherwise keep the plain
    # mean-field fixed point.
    # A non-default access graph (accost) modulates admission/promotion per item
    # (row 0 = miss admission per list, row a = hit-in-list-a promotion target).
    # The refined mean-field machinery is specialized to the linear chain, so a
    # non-linear graph uses the general drift with the plain mean-field fixed
    # point; the linear default keeps the 1/N-refined path unchanged.
    G = _build_item_graphs(accost, lambd, n, h)
    if G is not None:
        xss = _fixed_point_graph(x0, p, G, m, n, h, dim)
    else:
        xss = _fixed_point(x0, p, m, n, h, dim)
        try:
            pi_mf, V = _expansion_steady_state(x0, p, m, n, h, dim)
            xref = pi_mf + V / n
            if np.all(np.isfinite(xref)):
                xss = xref
        except Exception:
            pass  # fall back to plain mean field

    # per-item miss probability (occupancy of list 0), clipped to [0,1]
    pi0 = np.array([min(1.0, max(0.0, xss[_idx(i, 0, n)])) for i in range(n)])

    MI = lam_i * pi0
    MU = np.zeros(u)
    for v in range(u):
        row = np.array(lambd[v, :, 0], dtype=float)
        row[~np.isfinite(row)] = 0.0
        MU[v] = float(np.dot(row, pi0))
    M = float(np.sum(MI))

    if tspan is None:
        return M, MU, MI, pi0

    # Transient mean-field trajectory: integrate the plain drift over the
    # window from the supplied (or default) initial occupancy. xtraj is the
    # full DDPP occupancy trajectory used to carry the cache mean occupancy
    # across environment switches.
    x0t = np.asarray(x0, dtype=float) if x0init is None else np.asarray(x0init, dtype=float).ravel()
    t0, t1 = float(tspan[0]), float(tspan[1])
    t_eval = np.linspace(t0, t1, 200)
    sol = solve_ivp(lambda t, x: _drift(x, p, m, n, h, dim), (t0, t1), x0t,
                    method='LSODA', rtol=1e-8, atol=1e-10, t_eval=t_eval)
    tout = sol.t
    xtraj = sol.y  # dim x nt
    nt = tout.shape[0]
    pi0_t = np.clip(xtraj[[_idx(i, 0, n) for i in range(n)], :], 0.0, 1.0)
    MU_t = np.zeros((u, nt))
    for v in range(u):
        row = np.array(lambd[v, :, 0], dtype=float)
        row[~np.isfinite(row)] = 0.0
        MU_t[v, :] = row @ pi0_t
    return M, MU, MI, pi0, tout, pi0_t, MU_t, xtraj
