"""MAMAP(2,m) fitting: marks a second-order acyclic MAP so that the class
probabilities are matched exactly and the forward moments, backward moments and
one-step class transition probabilities as closely as the form allows.

Port of m3a ``mamap2m_fit_fb_multiclass.m``, ``mamap2m_fit_gamma_fb.m``,
``mamap2m_fit.m``, ``mamap2m_fit_mmap.m`` and ``mamap2m_fit_trace.m``. An MMAP
is a list ``[D0, D1, D11, ..., D1m]``.
"""
from typing import List, Optional, Sequence, Tuple

import numpy as np

from .amap2 import amap2_fit_gamma, amap2_fitall_gamma, _map_repair
from .maph2m import maph2m_fit, maph2m_fit_multiclass
from .qp import solve_qp

DEGENTOL = 1e-8
GAMMATOL = 1e-4


def _marked_poisson(D0: np.ndarray, D1: np.ndarray, p: np.ndarray) -> List[np.ndarray]:
    out = [np.asarray(D0, float).copy(), np.asarray(D1, float).copy()]
    for c in range(p.size):
        out.append(np.asarray(D1, float) * p[c])
    return out


def _fwd_bwd(mmap: Sequence[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    from line_solver.api.mam.mmap_ops import mmap_backward_moment, mmap_forward_moment
    fF = np.asarray(mmap_forward_moment(list(mmap), [1]), dtype=float).ravel()
    fB = np.asarray(mmap_backward_moment(list(mmap), [1]), dtype=float).ravel()
    return fF, fB


def mamap2m_fit_fb_multiclass(map_: Sequence[np.ndarray], p: Sequence[float],
                              F: Sequence[float], B: Sequence[float],
                              classWeights: Optional[Sequence[float]] = None,
                              fbWeights: Optional[Sequence[float]] = None
                              ) -> Tuple[List[np.ndarray], np.ndarray, np.ndarray]:
    """Mark a canonical acyclic MAP(2) on its forward and backward moments."""
    D0 = np.asarray(map_[0], dtype=float)
    D1 = np.asarray(map_[1], dtype=float)
    if D0.shape[0] != 2:
        raise ValueError('Underlying MAP must be of second-order.')
    if D0[1, 0] != 0:
        raise ValueError('Underlying MAP must be acyclic')
    if D1[0, 1] == 0:
        form = 1
    elif D1[0, 0] == 0:
        form = 2
    else:
        raise ValueError('Underlying MAP must be in canonical acyclic form')

    p = np.asarray(p, dtype=float).ravel()
    F = np.asarray(F, dtype=float).ravel()
    B = np.asarray(B, dtype=float).ravel()
    k = p.size
    cw = np.ones(k) if classWeights is None else np.asarray(classWeights, float).ravel()
    fbw = np.ones(2) if fbWeights is None else np.asarray(fbWeights, float).ravel()

    mmap: List[np.ndarray] = [D0.copy(), D1.copy()]

    h1 = -1.0 / D0[0, 0]
    h2 = -1.0 / D0[1, 1]
    r1 = D0[0, 1] * h1
    r2 = D1[1, 1] * h2

    poisson = ((form == 1 and (r1 < DEGENTOL or r2 > 1 - DEGENTOL
                               or abs(h2 - h1 * r2) < DEGENTOL
                               or abs(h1 - h2 + h2 * r1) < DEGENTOL))
               or (form == 2 and (r2 > 1 - DEGENTOL
                                  or abs(h1 - h2 + h2 * r1) < DEGENTOL
                                  or abs(h1 - h2 - h1 * r1 + h1 * r1 * r2) < DEGENTOL)))

    if poisson:
        from line_solver.api.mam.map_analysis import map_mean
        h = map_mean(D0, D1)
        out = [np.array([[-1.0 / h]]), np.array([[1.0 / h]])]
        for c in range(k):
            out.append(out[1] * p[c])
        fF, fB = _fwd_bwd(out)
        return out, fF, fB

    q = np.zeros((3, k))

    if form == 2 and r2 < DEGENTOL and abs(1 - r1) < DEGENTOL:
        # degenerate phase-type
        for c in range(k):
            q[0, c] = p[c]
            q[1, c] = p[c]
            q[2, c] = p[c]
    elif form == 1 and r2 < DEGENTOL:
        # canonical phase-type: fall back to the MAPH fit
        aph0 = D0.copy()
        aph1 = D1.copy()
        aph1[1, 1] = 0.0
        aph = _map_repair(aph0, aph1)
        out, _ = maph2m_fit_multiclass(aph, p, B, cw)
        fF, fB = _fwd_bwd(out)
        return out, fF, fB
    elif abs(1 - r1) < DEGENTOL:
        # non-canonical phase-type: only the forward moments can be fitted
        q_f = np.zeros((2, k))
        q_0 = np.zeros((2, k))
        for c in range(k):
            q_f[0, c] = p[c] * (-1.0 / ((h1 + h2 * (r1 - 1)) * (r2 - 1) * (r1 + r2 - r1 * r2)))
            q_0[0, c] = p[c] * (h2 / ((r2 - 1) * (r1 + r2 - r1 * r2) * (h1 - h2 + h2 * r1)))
            q_f[1, c] = p[c] * (-1.0 / (r2 * (h1 + h2 * (r1 - 1)) * (r1 + r2 - r1 * r2)))
            q_0[1, c] = p[c] * ((h1 + h2 * r1) / (r2 * (r1 + r2 - r1 * r2) * (h1 - h2 + h2 * r1)))
        fF = _solve_single(q_f, q_0, F, cw * fbw[0], k)
        for c in range(k):
            q[0, c] = 1.0 / k
            q[1, c] = fF[c] * q_f[0, c] + q_0[0, c]
            q[2, c] = fF[c] * q_f[1, c] + q_0[1, c]
    elif form == 2 and r2 < DEGENTOL:
        # degenerate form for gamma < 0: fit either forward or backward
        if fbw[0] >= fbw[1]:
            q_f = np.zeros((2, k))
            q_0 = np.zeros((2, k))
            for c in range(k):
                q_f[0, c] = p[c] * (-(r1 - 2) / ((h1 + h2 * (r1 - 1)) * (r1 - 1)))
                q_0[0, c] = p[c] * (1 - (h1 + h2) / ((r1 - 1) * (h1 - h2 + h2 * r1)))
                q_f[1, c] = p[c] * (-(r1 - 2) / (h1 + h2 * (r1 - 1)))
                q_0[1, c] = p[c] * ((h2 * (r1 - 2)) / (h1 - h2 + h2 * r1))
            x = _solve_single(q_f, q_0, F, cw * fbw[0], k)
            for c in range(k):
                for j in range(2):
                    q[j, c] = x[c] * q_f[j, c] + q_0[j, c]
                q[2, c] = 1.0 / k
        else:
            q_b = np.zeros((2, k))
            q_0 = np.zeros((2, k))
            for c in range(k):
                q_b[0, c] = p[c] * (-(r1 - 2) / ((h2 + h1 * (r1 - 1)) * (r1 - 1)))
                q_0[0, c] = p[c] * (1 - (h1 + h2) / ((r1 - 1) * (h2 - h1 + h1 * r1)))
                q_b[1, c] = p[c] * (-(r1 - 2) / (h2 + h1 * (r1 - 1)))
                q_0[1, c] = p[c] * ((h1 * (r1 - 2)) / (h2 - h1 + h1 * r1))
            x = _solve_single(q_b, q_0, B, cw * fbw[1], k)
            for c in range(k):
                for j in range(2):
                    q[j, c] = x[c] * q_b[j, c] + q_0[j, c]
                q[2, c] = 1.0 / k
    else:
        # general case: q(j,c) = F(c) q_f(j,c) + B(c) q_b(j,c) + q_0(j,c)
        q_f = np.zeros((3, k))
        q_b = np.zeros((3, k))
        q_0 = np.zeros((3, k))
        for c in range(k):
            pc = p[c]
            if form == 1:
                q_f[0, c] = 0.0
                q_b[0, c] = -(pc * (r1 * r2 - r2 + 1)) / ((h2 - h1 * r2) * (r1 - 1) * (r2 - 1))
                q_0[0, c] = (pc * (h1 + h2 - h1 * r2) * (r1 * r2 - r2 + 1)) / ((h2 - h1 * r2) * (r1 - 1) * (r2 - 1))
                q_f[1, c] = -(pc * (r1 * r2 - r2 + 1)) / (r1 * (h1 + h2 * (r1 - 1)) * (r2 - 1))
                q_b[1, c] = -(pc * (r1 * r2 - r2 + 1)) / (r1 * (h2 - h1 * r2) * (r2 - 1))
                q_0[1, c] = ((pc * (r1 * r2 - r2 + 1)) / ((r1 - 1) * (r2 - 1))
                             + (h1 * pc * (r1 * r2 - r2 + 1)) / (r1 * (h2 - h1 * r2) * (r2 - 1))
                             - (h1 * pc * (r1 * r2 - r2 + 1)) / (r1 * (h1 + h2 * (r1 - 1)) * (r1 - 1) * (r2 - 1)))
                q_f[2, c] = -(pc * (r1 * r2 - r2 + 1)) / (r1 * r2 * (h1 - h2 + h2 * r1))
                q_b[2, c] = 0.0
                q_0[2, c] = (pc * (h1 + h2 * r1) * (r1 * r2 - r2 + 1)) / (r1 * r2 * (h1 - h2 + h2 * r1))
            else:
                q_f[0, c] = 0.0
                q_b[0, c] = -(pc * (r1 + r2 - r1 * r2 - 2)) / ((r1 - 1) * (r2 - 1) * (h1 - h2 - h1 * r1 + h1 * r1 * r2))
                q_0[0, c] = (pc * (h2 + h1 * r1 - h1 * r1 * r2) * (r1 + r2 - r1 * r2 - 2)) / ((r1 - 1) * (r2 - 1) * (h1 - h2 - h1 * r1 + h1 * r1 * r2))
                q_f[1, c] = (pc * (r1 + r2 - r1 * r2 - 2)) / ((r2 - 1) * (h1 - h2 + h2 * r1))
                q_b[1, c] = 0.0
                q_0[1, c] = -(h2 * pc * (r1 + r2 - r1 * r2 - 2)) / ((r2 - 1) * (h1 - h2 + h2 * r1))
                q_f[2, c] = (pc * (r1 + r2 - r1 * r2 - 2)) / (r2 * (h1 + h2 * (r1 - 1)))
                q_b[2, c] = (pc * (r1 + r2 - r1 * r2 - 2)) / (r2 * (h1 - h2 - h1 * r1 + h1 * r1 * r2))
                q_0[2, c] = ((h1 * pc * (r1 + r2 - r1 * r2 - 2)) / (r2 * (h1 + h2 * (r1 - 1)) * (r1 - 1))
                             - (h1 * pc * (r1 + r2 - r1 * r2 - 2)) / (r2 * (h1 - h2 - h1 * r1 + h1 * r1 * r2))
                             - (pc * (r1 + r2 - r1 * r2 - 2)) / (r2 * (r1 - 1)))

        n = 2 * k
        A = np.zeros((6 * k, n))
        b = np.zeros(6 * k)
        for c in range(k):
            for j in range(3):
                row = c * 6 + j * 2
                col = c * 2
                A[row, col] = q_f[j, c]
                A[row, col + 1] = q_b[j, c]
                b[row] = 1 - q_0[j, c]
                A[row + 1, col] = -q_f[j, c]
                A[row + 1, col + 1] = -q_b[j, c]
                b[row + 1] = q_0[j, c]

        Aeq = np.zeros((3, n))
        beq = np.ones(3)
        for c in range(k):
            for j in range(3):
                Aeq[j, c * 2] = q_f[j, c]
                Aeq[j, c * 2 + 1] = q_b[j, c]
                beq[j] -= q_0[j, c]

        H = np.zeros((n, n))
        hv = np.zeros(n)
        for c in range(k):
            base = c * 2
            fw = cw[c] * fbw[0]
            bw = cw[c] * fbw[1]
            H[base, base] = 2.0 / F[c] ** 2 * fw
            H[base + 1, base + 1] = 2.0 / B[c] ** 2 * bw
            hv[base] = -2.0 / F[c] * fw
            hv[base + 1] = -2.0 / B[c] * bw

        x0 = np.empty(n)
        x0[0::2] = F
        x0[1::2] = B
        x = solve_qp(H, hv, A, b, Aeq, beq, lb=1e-6, ub=1e6, x0=x0)
        for c in range(k):
            for j in range(3):
                q[j, c] = x[c * 2] * q_f[j, c] + x[c * 2 + 1] * q_b[j, c] + q_0[j, c]

    for c in range(k):
        if form == 1:
            mask = np.array([[q[0, c], 0.0], [q[1, c], q[2, c]]])
        else:
            mask = np.array([[0.0, q[0, c]], [q[1, c], q[2, c]]])
        mmap.append(D1 * mask)

    fF, fB = _fwd_bwd(mmap)
    return mmap, fF, fB


def _solve_single(q_x: np.ndarray, q_0: np.ndarray, target: np.ndarray,
                  weights: np.ndarray, k: int) -> np.ndarray:
    """QP over one moment vector, shared by the degenerate branches."""
    A = np.zeros((4 * k, k))
    b = np.zeros(4 * k)
    for c in range(k):
        for j in range(2):
            row = c * 4 + j * 2
            A[row, c] = q_x[j, c]
            b[row] = 1 - q_0[j, c]
            A[row + 1, c] = -q_x[j, c]
            b[row + 1] = q_0[j, c]
    Aeq = np.zeros((2, k))
    beq = np.ones(2)
    for c in range(k):
        for j in range(2):
            Aeq[j, c] = q_x[j, c]
            beq[j] -= q_0[j, c]
    H = np.zeros((k, k))
    h = np.zeros(k)
    for c in range(k):
        H[c, c] = 2.0 / target[c] ** 2 * weights[c]
        h[c] = -2.0 / target[c] * weights[c]
    return solve_qp(H, h, A, b, Aeq, beq, lb=1e-6, ub=1e6, x0=np.asarray(target, float).copy())


def mamap2m_fit_gamma_fb(M1: float, M2: float, M3: float, GAMMA: float,
                         P: Sequence[float], F: Sequence[float], B: Sequence[float]
                         ) -> List[np.ndarray]:
    """MAMAP(2,m) fitted on the forward and backward moments only."""
    P = np.asarray(P, float).ravel()
    F = np.asarray(F, float).ravel()
    B = np.asarray(B, float).ravel()

    _, maps = amap2_fit_gamma(M1, M2, M3, GAMMA)
    if len(maps) == 1 and np.asarray(maps[0][0]).shape[0] == 1:
        return _marked_poisson(maps[0][0], maps[0][1], P)

    best, best_err = None, np.inf
    for m in maps:
        mmap, fF, fB = mamap2m_fit_fb_multiclass(m, P, F, B)
        err = float(np.sum((fF / F - 1) ** 2) + np.sum((fB / B - 1) ** 2))
        if err < best_err:
            best, best_err = mmap, err
    return best


def mamap2m_fit(M1: float, M2: float, M3: float, GAMMA: float,
                P: Sequence[float], F: Sequence[float], B: Sequence[float],
                S, fbsWeights: Optional[Sequence[float]] = None) -> List[np.ndarray]:
    """MAPH(2,m) or MAMAP(2,m) matching the inter-arrival moments and decay
    rate, the class probabilities (exactly) and, as far as the underlying form
    allows, the forward moments, backward moments and one-step class
    transition probabilities.

    Args:
        M1, M2, M3: moments of the inter-arrival times
        GAMMA: autocorrelation decay rate
        P: class probabilities
        F: first-order forward moments
        B: first-order backward moments
        S: one-step class transition probabilities (k x k)
        fbsWeights: weights of forward moments, backward moments and sigma
    """
    from line_solver.api.mam.mmap_ops import mmap_sigma

    P = np.asarray(P, float).ravel()
    F = np.asarray(F, float).ravel()
    B = np.asarray(B, float).ravel()
    S = np.atleast_2d(np.asarray(S, float))
    w = np.ones(3) if fbsWeights is None else np.asarray(fbsWeights, float).ravel()
    fbW = np.array([w[0], w[1]])
    m = P.size

    if m > 2:
        return mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B)

    if abs(GAMMA) < GAMMATOL:
        return maph2m_fit(M1, M2, M3, P, B)

    _, maps = amap2_fit_gamma(M1, M2, M3, GAMMA)

    if len(maps) == 1 and np.asarray(maps[0][0]).shape[0] == 1:
        # Poisson underlying process: perturb just above the exponential to
        # recover a second-order form, else return a marked Poisson process
        M2a = M2 * (1 + 1e-4)
        M3a = M3 * (M2a / M2) ** 1.5
        maps2 = amap2_fitall_gamma(M1, M2a, M3a, GAMMA)
        if maps2:
            maps = [_map_repair(x[0], x[1]) for x in maps2]
        else:
            return _marked_poisson(maps[0][0], maps[0][1], P)

    best, best_err = None, np.inf
    for mp in maps:
        D0 = np.asarray(mp[0], float)
        D1 = np.asarray(mp[1], float)
        h1 = -1.0 / D0[0, 0]
        h2 = -1.0 / D0[1, 1]
        r1 = h1 * D0[0, 1]
        r2 = h2 * D1[1, 1]

        fitted = None
        degen = True
        if GAMMA > 0:
            if r1 < DEGENTOL or (1 - r2) < DEGENTOL:
                raise RuntimeError('Fitting MAMAP(2,m): should not happen')
            elif abs(h2 - h1 * r2) < DEGENTOL:
                fitted = _fs_unavailable('h2 = h1*r2')
            elif abs(h1 - h2 + h2 * r1) < DEGENTOL:
                fitted = _bs_unavailable('h1 - h2 + h2*r1 = 0')
            elif (1 - r1) < DEGENTOL:
                fitted = _fs_unavailable('r1 = 1 (non-canonical APH(2))')
            elif r2 < DEGENTOL:
                fitted, _ = maph2m_fit_multiclass(mp, P, B)
            else:
                degen = False
        else:
            if (1 - r2) < DEGENTOL:
                raise RuntimeError('Fitting MAMAP(2,m): should not happen')
            elif abs(h1 - h2 - h1 * r1 + h1 * r1 * r2) < DEGENTOL:
                fitted = _fs_unavailable('h1 - h2 - h1*r1 + h1*r1*r2 = 0')
            elif abs(h1 - h2 + h2 * r1) < DEGENTOL:
                fitted = _bs_unavailable('h1 - h2 + h2*r1 = 0')
            elif r2 < DEGENTOL and (1 - r1) < DEGENTOL:
                fitted, _ = maph2m_fit_multiclass(mp, P, B)
            elif r2 < DEGENTOL:
                if w[0] >= w[1]:
                    fitted = _fs_unavailable('r2 = 0, forward preferred')
                else:
                    fitted = _bs_unavailable('r2 = 0, backward preferred')
            else:
                degen = False

        if not degen:
            if w[0] >= w[2] and w[1] >= w[2]:
                fitted, _, _ = mamap2m_fit_fb_multiclass(mp, P, F, B, None, fbW)
            elif w[0] >= w[1]:
                fitted = _fs_unavailable('forward and sigma preferred')
            else:
                fitted = _bs_unavailable('backward and sigma preferred')

        fF, fB = _fwd_bwd(fitted)
        fS = np.atleast_2d(np.asarray(mmap_sigma(list(fitted)), dtype=float))
        err = (w[0] * (F[0] / fF[0] - 1) ** 2 + w[1] * (B[0] / fB[0] - 1) ** 2
               + w[2] * (S[0, 0] / fS[0, 0] - 1) ** 2)
        if err < best_err:
            best, best_err = fitted, err

    return best


def _fs_unavailable(reason: str):
    raise NotImplementedError(
        'mamap2m_fit: this underlying form (%s) needs mamap22_fit_fs_multiclass, '
        'which is not yet ported to native Python (see _kb/03-api-layer.md)' % reason)


def _bs_unavailable(reason: str):
    raise NotImplementedError(
        'mamap2m_fit: this underlying form (%s) needs mamap22_fit_bs_multiclass, '
        'which is not yet ported to native Python (see _kb/03-api-layer.md)' % reason)


def mamap2m_fit_mmap(mmap: Sequence[np.ndarray]) -> List[np.ndarray]:
    """MAMAP(2,m) fitting the characteristics of a given MMAP."""
    from line_solver.api.mam.map_analysis import map_gamma, map_moment
    from line_solver.api.mam.mmap_ops import (mmap_backward_moment, mmap_forward_moment,
                                              mmap_pc, mmap_sigma)
    D0 = np.asarray(mmap[0], float)
    D1 = np.asarray(mmap[1], float)
    M1 = map_moment(D0, D1, 1)
    M2 = map_moment(D0, D1, 2)
    M3 = map_moment(D0, D1, 3)
    GAMMA = map_gamma(D0, D1)
    P = np.asarray(mmap_pc(list(mmap)), float).ravel()
    F = np.asarray(mmap_forward_moment(list(mmap), [1]), float).ravel()
    B = np.asarray(mmap_backward_moment(list(mmap), [1]), float).ravel()
    S = np.atleast_2d(np.asarray(mmap_sigma(list(mmap)), float))
    return mamap2m_fit(M1, M2, M3, GAMMA, P, F, B, S)


def mamap2m_fit_trace(T: Sequence[float], A: Sequence[int]) -> List[np.ndarray]:
    """MAMAP(2,m) fitting the characteristics of a marked trace."""
    from line_solver.api.trace.trace_analysis import (mtrace_backward_moment,
                                                      mtrace_forward_moment,
                                                      mtrace_sigma)
    T = np.asarray(T, dtype=float).ravel()
    A = np.asarray(A).ravel()
    classes = np.unique(A)
    M1 = float(np.mean(T))
    M2 = float(np.mean(T ** 2))
    M3 = float(np.mean(T ** 3))
    Tc = T - M1
    denom = float(np.sum(Tc * Tc))
    rho1 = float(np.sum(Tc[:-1] * Tc[1:]) / denom) if denom > 0 else 0.0
    rho2 = float(np.sum(Tc[:-2] * Tc[2:]) / denom) if denom > 0 else 0.0
    GAMMA = rho2 / rho1 if abs(rho1) > 1e-12 else 0.0
    P = np.array([np.mean(A == c) for c in classes], dtype=float)
    F = np.asarray(mtrace_forward_moment(T, A, [1]), dtype=float).ravel()
    B = np.asarray(mtrace_backward_moment(T, A, [1]), dtype=float).ravel()
    S = np.atleast_2d(np.asarray(mtrace_sigma(T, A), dtype=float))
    return mamap2m_fit(M1, M2, M3, GAMMA, P, F, B, S)
