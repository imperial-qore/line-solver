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


def amap2_adjust_gamma(M1: float, M2: float, M3: float, GAMMA: float,
                       weights: Optional[np.ndarray] = None, method: int = 3
                       ) -> Tuple[float, float, float]:
    """Nearest feasible (M2, M3, GAMMA) for an AMAP(2).

    Only the default method 3 is ported: adjust (M2, M3) with ``aph2_adjust``
    (priority M2 > M3 > GAMMA), then clamp GAMMA to the interval feasible at the
    adjusted moments. Methods 1, 2 and 4 of ``amap2_adjust_gamma.m`` drive
    MATLAB global optimizers (patternsearch, PSwarm) that have no counterpart
    here, so they are rejected rather than approximated.
    """
    if method != 3:
        raise NotImplementedError(
            'amap2_adjust_gamma: only method 3 is available in native Python; '
            'methods 1, 2 and 4 require the MATLAB global optimizers')
    from line_solver.api.mam.aph2_fitting import aph2_adjust
    M2a, M3a = aph2_adjust(M1, M2, M3)
    lb, ub = _gamma_bounds(M1, M2a, M3a)
    GAMMAa = max(lb, min(GAMMA, ub))
    return float(M2a), float(M3a), float(GAMMAa)


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
    acf = np.correlate(T - M1, T - M1, mode='full') / (len(T) * np.var(T))
    mid = len(T) - 1
    rho1 = float(acf[mid + 1])
    rho2 = float(acf[mid + 2])
    GAMMA = rho2 / rho1 if abs(rho1) > 1e-12 else 0.0
    return amap2_fit_gamma(M1, M2, M3, GAMMA)
