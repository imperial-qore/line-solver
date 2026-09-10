"""AMAP(2) fitting: acyclic MAPs of second order matching three moments and the
autocorrelation decay rate.

Port of m3a ``amap2_assemble.m``, ``amap2_fitall_gamma.m``,
``amap2_adjust_gamma.m`` and ``amap2_fit_gamma.m``. MAPs are represented as a
list ``[D0, D1]`` of numpy arrays, as elsewhere in the native Python m3a code.
"""
from typing import List, Optional, Tuple

import numpy as np

DEGENTOL = 1e-8
R12TOL = 1e-6
TOL = 1e-2


def amap2_assemble(l1: float, l2: float, p1: float, p2: float, form: int) -> List[np.ndarray]:
    """AMAP(2) with the given phase means and branching probabilities."""
    if form == 1:
        D0 = np.array([[-1.0 / l1, p1 / l1], [0.0, -1.0 / l2]])
        D1 = np.array([[(1 - p1) / l1, 0.0], [(1 - p2) / l2, p2 / l2]])
    elif form == 2:
        D0 = np.array([[-1.0 / l1, p1 / l1], [0.0, -1.0 / l2]])
        D1 = np.array([[0.0, (1 - p1) / l1], [(1 - p2) / l2, p2 / l2]])
    else:
        raise ValueError('Invalid form: should be either 1 (gamma > 0) or 2 (gamma < 0)')
    return [D0, D1]


def _map_repair(D0: np.ndarray, D1: np.ndarray) -> List[np.ndarray]:
    """Generator repair of MATLAB map_normalize.m: drop imaginary parts, clamp
    negatives to zero and close each row of D0 + D1.

    NOTE: this is NOT what ``line_solver.api.mam.map_normalize`` does -- that
    one rescales a MAP to unit mean, which would destroy the fitted moments.
    """
    D0 = np.real(np.asarray(D0, dtype=complex)).astype(float).copy()
    D1 = np.real(np.asarray(D1, dtype=complex)).astype(float).copy()
    D0[D0 < 0] = 0.0
    D1[D1 < 0] = 0.0
    n = D0.shape[0]
    for i in range(n):
        D0[i, i] = 0.0
        D0[i, i] = -(D0[i, :].sum() + D1[i, :].sum())
    return [D0, D1]


def _feasible(r1: float, r2: float) -> bool:
    return (np.isreal(r1) and np.isreal(r2)
            and r1 >= -R12TOL and r1 <= 1 + R12TOL
            and r2 >= -R12TOL and r2 <= 1 + R12TOL)


def _fix(q: float) -> float:
    return max(min(q, 1.0), 0.0)


def amap2_fitall_gamma(M1: float, M2: float, M3: float, GAMMA: float) -> List[List[np.ndarray]]:
    """Every AMAP(2) matching (M1, M2, M3) and the decay rate GAMMA exactly."""
    SCV = (M2 - M1 ** 2) / M1 ** 2
    M3lb = 3 * M1 ** 3 * (3 * SCV - 1 + np.sqrt(2) * (1 - SCV) ** 1.5) if SCV <= 1 else np.nan

    if SCV <= 1 and abs(M3 - M3lb) < DEGENTOL:
        # at the lower bound of M3 the argument of the square root is zero
        tmp0 = 0.0
    else:
        tmp0 = (M3 ** 2 / 9 + ((8 * M1 ** 3) / 3 - 2 * M2 * M1) * M3
                - 3 * M1 ** 2 * M2 ** 2 + 2 * M2 ** 3)
        if tmp0 < 0:
            return []

    tmp1 = 3 * np.sqrt(tmp0)
    tmp2 = M3 - 3 * M1 * M2
    tmp3 = 6 * M2 - 12 * M1 ** 2

    if tmp0 == 0:
        n = 1
        h2v = [tmp2 / tmp3]
        h1v = [tmp2 / tmp3]
    else:
        n = 2
        h2v = [(tmp2 + tmp1) / tmp3, (tmp2 - tmp1) / tmp3]
        h1v = [h2v[1], h2v[0]]

    if min(h2v) <= 0:
        return []

    amaps: List[List[np.ndarray]] = []
    for j in range(n):
        h1 = h1v[j]
        h2 = h2v[j]
        if GAMMA >= 0:
            # first canonical form
            z = (M1 ** 2 * GAMMA ** 2
                 + (2 * M1 * h1 + 2 * M1 * h2 - 4 * h1 * h2 - 2 * M1 ** 2) * GAMMA
                 + M1 ** 2 - 2 * M1 * h1 - 2 * M1 * h2 + h1 ** 2 + 2 * h1 * h2 + h2 ** 2)
            if abs(z) < DEGENTOL:
                r2 = (h1 - M1 + h2 + GAMMA * M1) / (2 * h1)
                r1 = (M1 - h1 - M1 * r2 + h1 * r2) / (h2 - M1 * r2)
                if _feasible(r1, r2):
                    amaps.append(amap2_assemble(h1, h2, _fix(r1), _fix(r2), 1))
            elif z > 0:
                r2v = [(h1 - M1 + h2 - np.sqrt(z) + GAMMA * M1) / (2 * h1),
                       (h1 - M1 + h2 + np.sqrt(z) + GAMMA * M1) / (2 * h1)]
                for r2 in r2v:
                    r1 = (M1 - h1 - M1 * r2 + h1 * r2) / (h2 - M1 * r2)
                    if _feasible(r1, r2):
                        amaps.append(amap2_assemble(h1, h2, _fix(r1), _fix(r2), 1))
        else:
            # second canonical form
            r2 = (h1 - M1 + h2 + GAMMA * M1) / h1
            r1 = (r2 + (h1 + h2 - h1 * r2) / M1 - 2) / (r2 - 1)
            if _feasible(r1, r2):
                amaps.append(amap2_assemble(h1, h2, _fix(r1), _fix(r2), 2))
    return amaps


def _gamma_bounds(M1: float, M2a: float, M3a: float) -> Tuple[float, float]:
    """Feasible interval of the decay rate at the given normalized moments."""
    N2a = M2a / M1 ** 2
    N3a = M3a / (M2a * M1)
    if N2a < 2:
        lb = -(N2a * (N3a - 6) + 6) / (3 * N2a - 6)
        ub = -(2 * (0.5 * (N2a - 2) + 0.5 * np.sqrt(N2a ** 2 - 2 * N2a * N3a / 3)) ** 2) / (N2a - 2)
        ub = ub * (1 - TOL)
    elif N3a < 9 - 12 / N2a:
        lb = -(N2a * (N3a - 6) + 6) / (3 * N2a - 6)
        ub = 1 - TOL
    else:
        x1 = np.sqrt(N2a * (N2a * (18 * N2a + N3a * (N3a - 18) - 27) + 24 * N3a))
        x2 = N2a * (N3a - 9)
        lb = (x2 - x1 + 12) / (x2 + x1 + 12)
        ub = 1 - TOL
    return float(lb), float(ub)


_ADJUST_TOL = 1e-2


def _nonlcon_theoretical(M1, xM2, xM3, xGAMMA):
    """The five feasibility residuals c(x) <= 0 of amap2_adjust_gamma.m."""
    n2 = xM2 / M1 ** 2
    n3 = xM3 / (M1 * xM2)

    # The moment-bound formulas are only defined on 3/2 <= n2 < 2, so they are
    # evaluated at a CLAMPED n2: an iterate below 3/2 (which constraint 1 pushes
    # back) or at exactly 2 would leave the third-moment rows at zero and M3
    # unconstrained while constraint 1 is satisfied to solver tolerance.
    n2b = min(max(n2, 1.5), 2 - 1e-12)
    p2 = 3 * (n2b - 2) / (3 * n2b) * (-2 * np.sqrt(3) / np.sqrt(12 - 6 * n2b) - 1)
    a2 = (n2b - 2) / (p2 * (1 - n2b) + np.sqrt(p2 ** 2 + (2 * p2 * (n2b - 2))))
    l2 = 3 * (a2 + 1) / (a2 * p2 + 1) - (6 * a2) / (2 + a2 * p2 * (2 * a2 + 2))
    u2 = 6 * (n2b - 1) / n2b

    eps = 1e-8
    c = np.zeros(5)
    c[0] = 1.5 - n2                       # CONSTRAINT 1: 3/2 <= n2
    if n2 <= 2:
        c[1] = l2 - n3                    # CONSTRAINT 2: l2 <= n3
        c[2] = n3 - u2                    # CONSTRAINT 3: n3 <= u2
    else:
        c[1] = 1.5 * n2 - n3 + eps        # CONSTRAINT 2: 3/2 n2 < n3
        c[2] = 0.0                        # CONSTRAINT 3: disabled

    with np.errstate(invalid='ignore', divide='ignore'):
        lb1 = -(n2 * (n3 - 6) + 6) / (3 * n2 - 6)
        ub1 = -(2 * (0.5 * (n2 - 2) + 0.5 * np.sqrt(n2 ** 2 - (2 * n2 * n3) / 3)) ** 2) \
            / (n2 - 2)
        tmp1 = n2 * (n3 - 9)
        tmp2 = np.sqrt(n2 * (n2 * (18 * n2 + n3 * (n3 - 18) - 27) + 24 * n3))
        lb2 = (tmp1 - tmp2 + 12) / (tmp1 + tmp2 + 12)
    if n2 < 2:
        c[3] = lb1 - xGAMMA
        c[4] = xGAMMA - ub1
    elif n2 > 2 and n3 < 9 - 12 / n2:
        c[3] = lb1 - xGAMMA
        c[4] = xGAMMA - 1
    elif n2 > 2:
        c[3] = lb2 - xGAMMA
        c[4] = xGAMMA - 1
    return np.nan_to_num(c, nan=1.0, posinf=1.0, neginf=-1.0)


def _constrained_min(fun, x0, bounds, cons_fun):
    """Minimize a smooth objective under c(x) <= 0 on a box.

    MATLAB drives these two searches with patternsearch, a derivative-free
    pattern search with a Nelder-Mead inner search. scipy has no patternsearch;
    SLSQP from the same feasible start, restarted from the box mid-point if it
    lands infeasible, solves the same program. The result is the nearest
    feasible point under the same objective and the same constraints, not the
    same iterate sequence.
    """
    from scipy.optimize import minimize, NonlinearConstraint

    nlc = NonlinearConstraint(cons_fun, -np.inf, 0.0)
    best, best_obj = None, np.inf
    starts = [np.asarray(x0, dtype=float)]
    mid = np.array([0.5 * (lo + (hi if np.isfinite(hi) else lo + 2 * abs(lo) + 1))
                    for lo, hi in bounds])
    starts.append(mid)
    for s in starts:
        try:
            res = minimize(fun, s, method='SLSQP', bounds=bounds,
                           constraints=[{'type': 'ineq',
                                         'fun': lambda z: -np.asarray(cons_fun(z))}],
                           options={'maxiter': 500, 'ftol': 1e-12})
        except Exception:
            continue
        if not np.all(np.isfinite(res.x)):
            continue
        if np.max(cons_fun(res.x)) > 1e-6:
            continue
        if res.fun < best_obj:
            best, best_obj = res.x, res.fun
    if best is None:
        # No feasible point was reached; the start is feasible by construction
        # in both callers, so returning it is the nearest point actually found.
        return np.asarray(x0, dtype=float)
    return best


def amap2_adjust_gamma(M1: float, M2: float, M3: float, GAMMA: float,
                       weights: Optional[np.ndarray] = None, method: int = 3
                       ) -> Tuple[float, float, float]:
    """Nearest feasible (M2, M3, GAMMA) for an AMAP(2).

    Args:
        M1, M2, M3: moments of the marginal distribution
        GAMMA: autocorrelation decay rate
        weights: three-element vector weighting M2, M3, GAMMA; default
            [10, 1, 10]. Ignored by methods 3 and 4's M2 handling.
        method: 1) search on (M2, M3, GAMMA);
            2) fit M2 as closely as possible, then search on (M3, GAMMA), the
               weight of M2 being ignored;
            3) (default) prioritize M2, then M3, then GAMMA;
            4) prioritize M2, then GAMMA, then M3, by a bounded 1-D search on
               M3 minimizing w2*(M3a/M3-1)^2 + w3*(GAMMAa/GAMMA-1)^2.

    Returns:
        (M2a, M3a, GAMMAa), a feasible triple; M1 is always feasible.

    Note:
        Methods 1, 2 and 4 are optimization searches. MATLAB runs them with
        patternsearch and PSwarm; here they run with the scipy equivalents
        (SLSQP under the same constraints, bounded scalar search plus
        differential evolution). The program solved is the same, so the results
        agree in feasibility and objective value but not iterate for iterate.
    """
    from line_solver.api.mam.aph2_fitting import aph2_adjust

    if weights is None:
        weights = np.array([10.0, 1.0, 10.0])
    weights = np.asarray(weights, dtype=float).ravel()
    tol = _ADJUST_TOL

    if method == 3:
        # priorities are M2 > M3 > GAMMA
        M2a, M3a = aph2_adjust(M1, M2, M3)
        lb, ub = _gamma_bounds(M1, M2a, M3a)
        GAMMAa = max(lb, min(GAMMA, ub))
        return float(M2a), float(M3a), float(GAMMAa)

    if method == 1:
        target = np.array([M2, M3, GAMMA], dtype=float)

        def fun(x):
            return float(np.linalg.norm((x - target) / target * weights))

        # Feasible start, as MATLAB builds it from amap2_assemble
        from line_solver.api.mam.map_analysis import map_scale, map_moment, map_acf
        if GAMMA > 0:
            seed = amap2_assemble(1, 1 / 3, 1 / 2, 2 / 3, 1)
        else:
            seed = amap2_assemble(1, 2 / 3, 1 / 2, 2 / 3, 2)
        f0, f1 = map_scale(seed[0], seed[1], M1)
        fm2 = map_moment(f0, f1, 2)
        fm3 = map_moment(f0, f1, 3)
        fgamma = float(np.asarray(map_acf(f0, f1, 4)).ravel()[0]
                       / np.asarray(map_acf(f0, f1, 3)).ravel()[0])
        x0 = np.array([fm2, fm3, fgamma], dtype=float)
        bounds = [(0.0, np.inf), (0.0, np.inf), (-1.0, 1 - tol)]
        x = _constrained_min(fun, x0, bounds,
                             lambda z: _nonlcon_theoretical(M1, z[0], z[1], z[2]))
        return float(x[0]), float(x[1]), float(x[2])

    if method == 2:
        # force feasibility of the second moment
        M2a = max(1.5 * M1 ** 2, M2)
        # The M3 interval must come from the ADJUSTED second moment: with the
        # unadjusted M2 the branch test and the bound formulas refer to a
        # normalized moment the search is not constrained to, which can return
        # an inverted interval.
        n2 = M2a / M1 ** 2
        p2 = 3 * (n2 - 2) / (3 * n2) * (-2 * np.sqrt(3) / np.sqrt(12 - 6 * n2) - 1) \
            if n2 < 2 else np.nan
        if 1.5 <= n2 < 2:
            a2 = (n2 - 2) / (p2 * (1 - n2) + np.sqrt(p2 ** 2 + (2 * p2 * (n2 - 2))))
            l2 = 3 * (a2 + 1) / (a2 * p2 + 1) - (6 * a2) / (2 + a2 * p2 * (2 * a2 + 2))
            u2 = 6 * (n2 - 1) / n2
            fm3 = (l2 + (u2 - l2) / 2) * M1 * M2a
            m3_lb = l2 * M1 * M2a
            m3_ub = u2 * M1 * M2a
        else:
            fm3 = 1.5 * M2a ** 2 / M1 + tol
            m3_lb = fm3
            m3_ub = np.inf
        target = np.array([M3, GAMMA], dtype=float)

        def fun2(x):
            return float(np.linalg.norm((x - target) / target * weights[1:3]))

        # a null decay rate is always feasible
        x0 = np.array([fm3, 0.0], dtype=float)
        bounds = [(m3_lb, m3_ub), (-1.0, 1 - tol)]
        x = _constrained_min(fun2, x0, bounds,
                             lambda z: _nonlcon_theoretical(M1, M2a, z[0], z[1]))
        return float(M2a), float(x[0]), float(x[1])

    if method == 4:
        # priorities are M2 > GAMMA > M3
        M1sq = M1 ** 2
        scv = (M2 - M1sq) / M1sq
        if scv < 0.5:
            M2a = 1.5 * M1sq
            scva = (M2a - M1sq) / M1sq
        else:
            M2a = M2
            scva = scv
        if scva <= 1:
            m3_lb = 3 * M1 ** 3 * (3 * scva - 1 + np.sqrt(2) * (1 - scva) ** 1.5)
            m3_ub = 6 * M1 ** 3 * scva
        else:
            m3_lb = 1.5 * M1 ** 3 * (1 + scva) ** 2
            m3_ub = np.inf

        if abs(m3_lb - m3_ub) < tol:
            # exponential
            M3a = 0.5 * (m3_lb + m3_ub)
            return float(M2a), float(M3a), 0.0

        def gamma_objective(M3a_val):
            M3a_val = float(np.asarray(M3a_val).ravel()[0])
            lb, ub = _gamma_bounds(M1, M2a, M3a_val)
            gam = max(lb, min(GAMMA, ub))
            return (weights[1] * (M3a_val / M3 - 1) ** 2
                    + weights[2] * (gam / GAMMA - 1) ** 2)

        lo = m3_lb + tol
        hi = m3_ub if np.isfinite(m3_ub) else max(lo * 10.0, 10.0 * M3)
        # PSwarm is a bounded global search on this one variable; scipy's
        # differential_evolution is, refined by a bounded scalar search.
        from scipy.optimize import differential_evolution, minimize_scalar
        with np.errstate(invalid='ignore', divide='ignore'):
            res = differential_evolution(lambda z: gamma_objective(z[0]),
                                         [(lo, hi)], tol=1e-10, seed=0,
                                         polish=False)
            ref = minimize_scalar(gamma_objective, bounds=(lo, hi),
                                  method='bounded',
                                  options={'xatol': 1e-12})
        M3a = float(res.x[0]) if res.fun <= ref.fun else float(ref.x)
        lb, ub = _gamma_bounds(M1, M2a, M3a)
        GAMMAa = max(lb, min(GAMMA, ub))
        return float(M2a), float(M3a), float(GAMMAa)

    raise ValueError('amap2_adjust_gamma: invalid method %r for adjusting '
                     'AMAP(2) characteristics' % (method,))


def amap2_fit_gamma(M1: float, M2: float, M3: float, GAMMA: float
                    ) -> Tuple[List[np.ndarray], List[List[np.ndarray]]]:
    """AMAP(2) fitting the given moments and decay rate.

    Returns (AMAP, AMAPS): the preferred fit and every candidate form. Falls
    back to the nearest feasible characteristics, and finally to a Poisson
    process, exactly as ``amap2_fit_gamma.m`` does.
    """
    if abs(M2 - 2 * M1 ** 2) < 1e-6:
        # coefficient of variation one: marked Poisson process
        amap = [np.array([[-1.0 / M1]]), np.array([[1.0 / M1]])]
        return amap, [amap]

    amaps = amap2_fitall_gamma(M1, M2, M3, GAMMA)
    amaps = [_map_repair(m[0], m[1]) for m in amaps]

    if not amaps:
        M2a, M3a, GAMMAa = amap2_adjust_gamma(M1, M2, M3, GAMMA)
        amaps = amap2_fitall_gamma(M1, M2a, M3a, GAMMAa)

    if not amaps:
        amap = [np.array([[-1.0 / M1]]), np.array([[1.0 / M1]])]
        return amap, [amap]

    return amaps[0], amaps


def amap2_fit_gamma_map(MAP) -> Tuple[List[np.ndarray], List[List[np.ndarray]]]:
    """AMAP(2) fitting the characteristics of another MAP."""
    from line_solver.api.mam.map_analysis import map_gamma, map_moment
    D0, D1 = (np.asarray(MAP[0], float), np.asarray(MAP[1], float))
    return amap2_fit_gamma(map_moment(D0, D1, 1), map_moment(D0, D1, 2),
                           map_moment(D0, D1, 3), map_gamma(D0, D1))


def amap2_fit_gamma_trace(T) -> Tuple[List[np.ndarray], List[List[np.ndarray]]]:
    """AMAP(2) fitting the characteristics of a trace of inter-arrival times."""
    T = np.asarray(T, dtype=float).ravel()
    M1 = float(np.mean(T))
    M2 = float(np.mean(T ** 2))
    M3 = float(np.mean(T ** 3))
    # amap2_fit_gamma_trace.m reads the decay rate off trace_gamma, the robust
    # geometric fit of the whole ACF, not off a lag-2 / lag-1 ratio.
    from line_solver.api.trace.trace_analysis import trace_gamma
    GAMMA = float(np.asarray(trace_gamma(T), dtype=float).ravel()[0])
    return amap2_fit_gamma(M1, M2, M3, GAMMA)
