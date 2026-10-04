"""
M3PP: Marked Markov Modulated Poisson Processes.

Native Python implementations for M3PP (second-order Marked MMPP) fitting
and manipulation algorithms.

References:
    A. Horvath, G. Horvath, M. Telek, "A traffic based decomposition of two-class
    queueing networks with priority service", Computer Networks 2013.
"""

import numpy as np
from typing import Union, Tuple, Optional, List
from scipy.optimize import minimize, fsolve

ArrayLike = Union[np.ndarray, list]


def m3pp_rand(order: int, classes: int) -> List[np.ndarray]:
    """
    Generate a random M3PP (Marked Markov Modulated Poisson Process).

    Args:
        order: Number of phases
        classes: Number of arrival classes

    Returns:
        MMAP representation [D0, D1, D1_1, D1_2, ..., D1_classes]
    """
    from line_solver.lib.kpctoolbox import mmpp_rand
    from line_solver.api.mam import map_issym

    MAP = mmpp_rand(order)
    MMAP = [MAP[0], MAP[1]]

    for _ in range(classes):
        p = np.random.rand(classes)
        p = p / np.sum(p)
        for c in range(classes):
            MMAP.append(MAP[1] * p[c])

    return MMAP


def m3pp2m_interleave(m3pps: List[List[np.ndarray]]) -> List[np.ndarray]:
    """
    Compute the interleaved MMAP from multiple M3PP(2,m).

    Args:
        m3pps: List of M3PP processes to interleave

    Returns:
        Interleaved M3PP as [D0, D1, D1_1, ..., D1_M]
    """
    L = len(m3pps)

    # Compute off-diagonal rates
    r = np.zeros((2, L))
    r[0, L - 1] = m3pps[L - 1][0][0, 1]
    for i in range(L - 2, -1, -1):
        r[0, i] = m3pps[i][0][0, 1] - np.sum(r[0, i + 1:])

    r[1, 0] = m3pps[0][0][1, 0]
    for i in range(1, L):
        r[1, i] = m3pps[i][0][1, 0] - np.sum(r[1, :i])

    # Total number of classes
    M = sum(len(m3pp) - 2 for m3pp in m3pps)

    # State space size
    n = 2 + (L - 1)

    # Initialize result
    s = [None] * (2 + M)

    # Build D0 (off-diagonal part)
    D0 = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if j > i:
                D0[i, j] = r[0, j - 1]
            elif j < i:
                D0[i, j] = r[1, j]

    s[0] = D0

    # Build D1_c matrices
    c = 0
    for i in range(L):
        m = len(m3pps[i]) - 2  # Number of classes in i-th M3PP
        for j in range(m):
            D1c = np.zeros((n, n))
            for h in range(n):
                if h <= i:
                    D1c[h, h] = m3pps[i][2 + j][0, 0]
                else:
                    D1c[h, h] = m3pps[i][2 + j][1, 1]
            s[2 + c] = D1c
            c += 1

    # Build D1 (sum of all D1_c)
    D1 = np.zeros((n, n))
    for i in range(M):
        D1 = D1 + s[2 + i]
    s[1] = D1

    # Set diagonal of D0 to make rows sum to zero
    for h in range(n):
        s[0][h, h] = -np.sum(s[0][h, :] + s[1][h, :])

    return s


def _m3pp_split_coeffs_dv(l1: float, l2: float, r1: float, r2: float,
                          t: float) -> Tuple[float, float, float, float, float, float]:
    """
    Affine coefficients of the per-class marking probabilities for the
    variance-difference split of ``m3pp2m_fitc_approx.m``.

    Returns ``(q1_a, q1_x, q1_c, q2_a, q2_x, q2_c)`` with
    ``q1i = q1_a*ai + q1_x*dvi + q1_c`` and likewise for q2i.
    """
    # sinh(u)*exp(-u) with u = (r1+r2)*t/2, written without the cancellation
    sh = (1.0 - np.exp(-(r1 + r2) * t)) / 2.0

    br = 2 * l1 * sh - 2 * l2 * sh - l1 * r1 * t - l1 * r2 * t + l2 * r1 * t + l2 * r2 * t
    den1 = l1 * r2 * (r1 + r2) * br
    br2 = (l2 ** 2 * r1 ** 2 * t - 2 * l2 ** 2 * r1 * sh - l1 * l2 * r1 ** 2 * t
           + l2 ** 2 * r1 * r2 * t + 2 * l1 * l2 * r1 * sh - l1 * l2 * r1 * r2 * t)
    den2 = (r1 + r2) * br2

    p4 = r1 ** 4 + 4 * r1 ** 3 * r2 + 6 * r1 ** 2 * r2 ** 2 + 4 * r1 * r2 ** 3 + r2 ** 4
    cn = (l1 * r2 ** 4 * t + l2 * r1 ** 4 * t + 3 * l1 * r1 * r2 ** 3 * t
          + l1 * r1 ** 3 * r2 * t + l2 * r1 * r2 ** 3 * t + 3 * l2 * r1 ** 3 * r2 * t
          + 3 * l1 * r1 ** 2 * r2 ** 2 * t + 2 * l1 ** 2 * r1 * r2 ** 2 * t
          + 2 * l1 ** 2 * r1 ** 2 * r2 * t + 3 * l2 * r1 ** 2 * r2 ** 2 * t
          + 2 * l2 ** 2 * r1 * r2 ** 2 * t + 2 * l2 ** 2 * r1 ** 2 * r2 * t
          - 4 * l1 ** 2 * r1 * r2 * sh - 4 * l2 ** 2 * r1 * r2 * sh
          - 4 * l1 * l2 * r1 * r2 ** 2 * t - 4 * l1 * l2 * r1 ** 2 * r2 * t
          + 8 * l1 * l2 * r1 * r2 * sh)

    q1_a = (r1 ** 4 * t / 2 + r2 ** 4 * t / 2 - l1 * r2 ** 3 * t + l2 * r2 ** 3 * t
            + 2 * r1 * r2 ** 3 * t + 2 * r1 ** 3 * r2 * t + 3 * r1 ** 2 * r2 ** 2 * t
            + 2 * l1 * r2 ** 2 * sh - 2 * l2 * r2 ** 2 * sh - 2 * l1 * r1 * r2 ** 2 * t
            - l1 * r1 ** 2 * r2 * t + 2 * l2 * r1 * r2 ** 2 * t + l2 * r1 ** 2 * r2 * t
            + 2 * l1 * r1 * r2 * sh - 2 * l2 * r1 * r2 * sh) / den1
    q1_x = -p4 / (4 * den1)
    q1_c = -cn / (4 * den1)
    q2_a = -(r1 ** 4 * t / 2 + r2 ** 4 * t / 2 + l1 * r1 ** 3 * t - l2 * r1 ** 3 * t
             + 2 * r1 * r2 ** 3 * t + 2 * r1 ** 3 * r2 * t + 3 * r1 ** 2 * r2 ** 2 * t
             - 2 * l1 * r1 ** 2 * sh + 2 * l2 * r1 ** 2 * sh + l1 * r1 * r2 ** 2 * t
             + 2 * l1 * r1 ** 2 * r2 * t - l2 * r1 * r2 ** 2 * t - 2 * l2 * r1 ** 2 * r2 * t
             - 2 * l1 * r1 * r2 * sh + 2 * l2 * r1 * r2 * sh) / den2
    q2_x = p4 / (4 * den2)
    q2_c = cn / (4 * den2)
    return q1_a, q1_x, q1_c, q2_a, q2_x, q2_c


def _m3pp_split_coeffs_ag(l1: float, l2: float, r1: float, r2: float,
                          t: float) -> Tuple[float, float, float, float, float, float]:
    """
    Affine coefficients of the per-class marking probabilities for the
    variance-plus-covariance split of ``m3pp2m_fitc_approx_ag_multiclass.m``.
    Both constant terms are zero there.
    """
    e = np.exp(-r1 * t - r2 * t)
    f1num = l1 * r2 * (2 * l2 * r1 - 2 * l1 * r1 + r1 ** 3 * t + r2 ** 3 * t
                       + 2 * l1 * r1 ** 2 * t - 2 * l2 * r1 ** 2 * t
                       + 3 * r1 * r2 ** 2 * t + 3 * r1 ** 2 * r2 * t
                       + 2 * l1 * r1 * e - 2 * l2 * r1 * e
                       + 2 * l1 * r1 * r2 * t - 2 * l2 * r1 * r2 * t)
    f1 = f1num / (r1 + r2) ** 4
    f2num = l2 * r1 * (2 * l1 * r2 - 2 * l2 * r2 + r1 ** 3 * t + r2 ** 3 * t
                       - 2 * l1 * r2 ** 2 * t + 2 * l2 * r2 ** 2 * t
                       + 3 * r1 * r2 ** 2 * t + 3 * r1 ** 2 * r2 * t
                       - 2 * l1 * r2 * e + 2 * l2 * r2 * e
                       - 2 * l1 * r1 * r2 * t + 2 * l2 * r1 * r2 * t)
    f2 = f2num / (r1 + r2) ** 4
    tmp = f1 * l2 * r1 - f2 * l1 * r2
    if tmp == 0.0:
        raise ValueError('m3pp2m_fitc_approx_ag: degenerate per-class split')
    return (-(f2 * (r1 + r2)) / tmp, (l2 * r1) / tmp, 0.0,
            (f1 * (r1 + r2)) / tmp, -(l1 * r2) / tmp, 0.0)


def _m3pp_split_solve(coeffs, ai: np.ndarray, target: np.ndarray,
                      a: float) -> np.ndarray:
    """
    Least-squares compromise over the per-class targets, subject to the marking
    probabilities being non-negative and summing to one in each phase:

        min  sum_i (x_i/target_i - 1)^2
        s.t. q1i(x_i) >= 0, q2i(x_i) >= 0,  sum_i q1i(x_i) = sum_i q2i(x_i) = 1

    which is the reference's quadprog program written in relative form (the
    reference objective 1/2 x'Hx + f'x with H = diag(2/target^2), f = -2/target
    is this sum minus the constant m).
    """
    q1_a, q1_x, q1_c, q2_a, q2_x, q2_c = coeffs
    ai = np.asarray(ai, dtype=float).ravel()
    target = np.asarray(target, dtype=float).ravel()
    m = len(ai)
    if np.any(target == 0.0):
        raise ValueError('m3pp2m_fitc_approx: a per-class target is zero, so the '
                         'relative objective of the reference is undefined')

    # Neither q1_x nor q2_x depends on i, so both equality rows constrain sum_i
    # x_i alone; the reference's algebra makes them consistent, and the program
    # reduces to one equality plus a box.
    if q1_x == 0.0 or q2_x == 0.0:
        raise ValueError('m3pp2m_fitc_approx: the per-class split does not depend '
                         'on the free variables')
    S1 = (1.0 - m * q1_c - q1_a * a) / q1_x
    S2 = (1.0 - m * q2_c - q2_a * a) / q2_x
    if abs(S1 - S2) > 1e-6 * max(1.0, abs(S1), abs(S2)):
        raise ValueError('m3pp2m_fitc_approx: the two phase-normalization rows '
                         'disagree, so the split has no solution')
    S = 0.5 * (S1 + S2)

    lo = np.full(m, -np.inf)
    hi = np.full(m, np.inf)
    for qa, qx, qc in ((q1_a, q1_x, q1_c), (q2_a, q2_x, q2_c)):
        bound = -(qa * ai + qc) / qx
        if qx > 0:
            lo = np.maximum(lo, bound)
        else:
            hi = np.minimum(hi, bound)
    if np.sum(np.maximum(lo, -1e300)) > S + 1e-9 or np.sum(np.minimum(hi, 1e300)) < S - 1e-9:
        raise ValueError('m3pp2m_fitc_approx: empty feasibility region for the '
                         'per-class split')

    # KKT of min sum (x_i/t_i - 1)^2 s.t. sum x_i = S and lo <= x <= hi:
    # x_i(mu) = t_i (1 + mu t_i / 2) clipped to the box, increasing in mu.
    def xof(mu):
        return np.clip(target * (1.0 + mu * target / 2.0), lo, hi)

    mu_lo, mu_hi = -1.0, 1.0
    while np.sum(xof(mu_lo)) > S:
        mu_lo *= 2.0
        if mu_lo < -1e18:
            break
    while np.sum(xof(mu_hi)) < S:
        mu_hi *= 2.0
        if mu_hi > 1e18:
            break
    for _ in range(200):
        mu = 0.5 * (mu_lo + mu_hi)
        if np.sum(xof(mu)) < S:
            mu_lo = mu
        else:
            mu_hi = mu
    return xof(0.5 * (mu_lo + mu_hi))


def _m3pp2m_trivial_split(D0: np.ndarray, D1: np.ndarray, ai: np.ndarray,
                          a: float) -> Optional[List[np.ndarray]]:
    """Poisson and single-class short-circuits shared by the split entry points."""
    m = len(ai)
    if D0.shape[0] == 1:
        pi = ai / a
        return [D0, D1] + [pi[i] * D1 for i in range(m)]
    if m == 1:
        return [D0, D1, D1.copy()]
    return None


def m3pp2m_fitc_approx_ag_multiclass(mmpp: List[np.ndarray],
                                      ac: np.ndarray,
                                      gtc: np.ndarray,
                                      t: float) -> List[np.ndarray]:
    """
    Fit an M3PP(2,m) given the underlying MMPP(2), matching the per-class rates
    exactly and the per-class variance-plus-covariance ``gtc`` at scale ``t`` in
    the least-squares sense. Port of MATLAB
    ``m3a/m3pp/m3pp2m_fitc_approx_ag_multiclass.m``.

    Args:
        mmpp: MMPP(2) as [D0, D1]
        ac: Per-class arrival rates, which must sum to the rate of mmpp
        gtc: Per-class variance + covariance with all other classes, at scale t
        t: Time scale

    Returns:
        M3PP(2,m) as [D0, D1, D1_1, ..., D1_m]
    """
    from line_solver.api.mam import map_count_mean

    ac = np.asarray(ac, dtype=float).ravel()
    gtc = np.asarray(gtc, dtype=float).ravel()
    m = len(ac)

    D0 = np.atleast_2d(np.asarray(mmpp[0], dtype=float))
    D1 = np.atleast_2d(np.asarray(mmpp[1], dtype=float))

    a = float(np.ravel(map_count_mean(D0, D1, 1.0))[0])
    if abs(a - np.sum(ac)) > 1e-8:
        raise ValueError('Inconsistent per-class arrival rates.')

    trivial = _m3pp2m_trivial_split(D0, D1, ac, a)
    if trivial is not None:
        return trivial

    coeffs = _m3pp_split_coeffs_ag(D1[0, 0], D1[1, 1], D0[0, 1], D0[1, 0], t)
    x = _m3pp_split_solve(coeffs, ac, gtc, a)
    q1_a, q1_x, q1_c, q2_a, q2_x, q2_c = coeffs

    FIT = [D0, D1]
    for i in range(m):
        q1 = q1_a * ac[i] + q1_x * x[i] + q1_c
        q2 = q2_a * ac[i] + q2_x * x[i] + q2_c
        FIT.append(D1 * np.array([[q1, 0.0], [0.0, q2]]))
    return FIT


def m3pp2m_fitc_approx(a: float, bt1: float, bt2: float, binf: float,
                        m3t2: float, t1: float, t2: float,
                        ai: np.ndarray, dvt3: np.ndarray, t3: float
                        ) -> List[np.ndarray]:
    """
    Fit a second-order Marked MMPP, splitting the classes on their count
    variance DIFFERENCES at t3. Port of MATLAB ``m3a/m3pp/m3pp2m_fitc_approx.m``.

    The free variables are variance differences, which are routinely negative (a
    minority class has less variance than all the others combined) and whose sum
    is fixed by the equality rows, so no lower bound is imposed on them; this
    follows the sibling ``m3pp2m_fitc_approx_ag_multiclass``, whose quadprog call
    passes no box.

    Args:
        a: Total arrival rate
        bt1: IDC at scale t1
        bt2: IDC at scale t2
        binf: IDC for t->inf
        m3t2: Third central moment
        t1: First time scale
        t2: Second time scale
        ai: Per-class arrival rates
        dvt3: Per-class variance differences at resolution t3
        t3: Third time scale

    Returns:
        M3PP as [D0, D1, D1_1, ..., D1_m]
    """
    from line_solver.lib.kpctoolbox import mmpp2_fitc_approx

    ai = np.asarray(ai, dtype=float).ravel()
    dvt3 = np.asarray(dvt3, dtype=float).ravel()
    m = len(ai)

    if abs(a - np.sum(ai)) > 1e-8:
        raise ValueError("Inconsistent per-class arrival rates")

    FIT = mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2)
    D0 = np.atleast_2d(np.asarray(FIT[0], dtype=float))
    D1 = np.atleast_2d(np.asarray(FIT[1], dtype=float))

    trivial = _m3pp2m_trivial_split(D0, D1, ai, a)
    if trivial is not None:
        return trivial

    coeffs = _m3pp_split_coeffs_dv(D1[0, 0], D1[1, 1], D0[0, 1], D0[1, 0], t3)
    x = _m3pp_split_solve(coeffs, ai, dvt3, a)
    q1_a, q1_x, q1_c, q2_a, q2_x, q2_c = coeffs

    result = [D0, D1]
    for i in range(m):
        q1 = q1_a * ai[i] + q1_x * x[i] + q1_c
        q2 = q2_a * ai[i] + q2_x * x[i] + q2_c
        result.append(D1 * np.array([[q1, 0.0], [0.0, q2]]))
    return result


def m3pp2m_fitc_approx_ag(a: float, bt1: float, bt2: float, binf: float,
                          m3t2: float, t1: float, t2: float,
                          ai: np.ndarray, gt3: np.ndarray, t3: float
                          ) -> List[np.ndarray]:
    """
    Fit the underlying MMPP(2) by optimization, then apply the 'ag' per-class
    split. Port of MATLAB ``m3a/m3pp/m3pp2m_fitc_approx_ag.m``.

    Args:
        a, bt1, bt2, binf, m3t2, t1, t2: aggregate counting characteristics
        ai: Per-class arrival rates, which must sum to a
        gt3: Per-class variance + covariance with all other classes, at t3
        t3: Third time scale

    Returns:
        M3PP as [D0, D1, D1_1, ..., D1_m]
    """
    from line_solver.lib.kpctoolbox import mmpp2_fitc_approx

    ai = np.asarray(ai, dtype=float).ravel()
    if abs(a - np.sum(ai)) > 1e-8:
        raise ValueError('Inconsistent per-class arrival rates.')

    FIT = mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2)
    return m3pp2m_fitc_approx_ag_multiclass(list(FIT), ai, gt3, t3)


def m3pp2m_fitc(a: float, bt1: float, bt2: float, binf: float, m3t2: float,
                t1: float, t2: float, ai: np.ndarray, dvt3: np.ndarray,
                t3: float) -> List[np.ndarray]:
    """
    Fit a second-order Marked MMPP (exact per-class delta-variance matching).

    Faithful port of MATLAB ``m3pp2m_fitc.m``. The underlying MMPP(2) that
    matches the joint-process moments (a, bt1, bt2, binf, m3t2) is obtained with
    ``mmpp2_fitc``; the per-class marking probabilities q(:,i) are then set with
    the closed-form solution of the exact per-class count-variance-difference
    constraint at resolution t3.

    Args:
        a: arrival rate
        bt1, bt2, binf: IDC at scales t1, t2 and for t->inf
        m3t2: third central moment at scale t2
        t1, t2, t3: time scales
        ai: (m,) per-class arrival rates
        dvt3: (m,) delta between the variance of class i and the variance of all
              other classes combined, at resolution t3

    Returns:
        M3PP as ``[D0, D1, D1_1, ..., D1_m]``, or ``[]`` if the joint-process
        fit is infeasible.
    """
    ai = np.asarray(ai, dtype=float).ravel()
    dvt3 = np.asarray(dvt3, dtype=float).ravel()
    m = len(ai)

    from line_solver.lib.kpctoolbox import mmpp2_fitc
    FIT = mmpp2_fitc(a, bt1, bt2, binf, m3t2, t1, t2)
    if FIT is None:
        return []
    D0, D1 = np.atleast_2d(FIT[0]), np.atleast_2d(FIT[1])
    if D0.size == 0:
        return []

    # Degenerate case: marked Poisson process (single-phase MMPP).
    if D0.shape[0] == 1:
        result = [D0, D1]
        for i in range(m):
            result.append(np.array([[ai[i]]]))
        return result

    l1 = D1[0, 0]
    l2 = D1[1, 1]
    r1 = D0[0, 1]
    r2 = D0[1, 0]

    t = t3
    x = (r1 * t) / 2.0 + (r2 * t) / 2.0
    se = np.sinh(x) * np.exp(-x)   # sinh(x)*exp(-x), matching MATLAB verbatim

    q = np.zeros((2, m))
    for i in range(m - 1):
        a_1 = ai[i]
        dv_1 = dvt3[i]
        q[0, i] = -(dv_1*r1**4 + dv_1*r2**4 - 2*a_1*r1**4*t - 2*a_1*r2**4*t
                    + 4*dv_1*r1*r2**3 + 4*dv_1*r1**3*r2 + l1*r2**4*t + l2*r1**4*t
                    + 6*dv_1*r1**2*r2**2 + 4*a_1*l1*r2**3*t - 4*a_1*l2*r2**3*t
                    - 8*a_1*r1*r2**3*t - 8*a_1*r1**3*r2*t + 3*l1*r1*r2**3*t
                    + l1*r1**3*r2*t + l2*r1*r2**3*t + 3*l2*r1**3*r2*t
                    - 12*a_1*r1**2*r2**2*t + 3*l1*r1**2*r2**2*t + 2*l1**2*r1*r2**2*t
                    + 2*l1**2*r1**2*r2*t + 3*l2*r1**2*r2**2*t + 2*l2**2*r1*r2**2*t
                    + 2*l2**2*r1**2*r2*t - 8*a_1*l1*r2**2*se + 8*a_1*l2*r2**2*se
                    - 4*l1**2*r1*r2*se - 4*l2**2*r1*r2*se + 8*a_1*l1*r1*r2**2*t
                    + 4*a_1*l1*r1**2*r2*t - 8*a_1*l2*r1*r2**2*t - 4*a_1*l2*r1**2*r2*t
                    - 4*l1*l2*r1*r2**2*t - 4*l1*l2*r1**2*r2*t - 8*a_1*l1*r1*r2*se
                    + 8*a_1*l2*r1*r2*se + 8*l1*l2*r1*r2*se) \
                   / (4*l1*r2*(r1 + r2)*(2*l1*se - 2*l2*se - l1*r1*t - l1*r2*t
                                         + l2*r1*t + l2*r2*t))
        q[1, i] = (dv_1*r1**4 + dv_1*r2**4 - 2*a_1*r1**4*t - 2*a_1*r2**4*t
                   + 4*dv_1*r1*r2**3 + 4*dv_1*r1**3*r2 + l1*r2**4*t + l2*r1**4*t
                   + 6*dv_1*r1**2*r2**2 - 4*a_1*l1*r1**3*t + 4*a_1*l2*r1**3*t
                   - 8*a_1*r1*r2**3*t - 8*a_1*r1**3*r2*t + 3*l1*r1*r2**3*t
                   + l1*r1**3*r2*t + l2*r1*r2**3*t + 3*l2*r1**3*r2*t
                   - 12*a_1*r1**2*r2**2*t + 3*l1*r1**2*r2**2*t + 2*l1**2*r1*r2**2*t
                   + 2*l1**2*r1**2*r2*t + 3*l2*r1**2*r2**2*t + 2*l2**2*r1*r2**2*t
                   + 2*l2**2*r1**2*r2*t + 8*a_1*l1*r1**2*se - 8*a_1*l2*r1**2*se
                   - 4*l1**2*r1*r2*se - 4*l2**2*r1*r2*se - 4*a_1*l1*r1*r2**2*t
                   - 8*a_1*l1*r1**2*r2*t + 4*a_1*l2*r1*r2**2*t + 8*a_1*l2*r1**2*r2*t
                   - 4*l1*l2*r1*r2**2*t - 4*l1*l2*r1**2*r2*t + 8*a_1*l1*r1*r2*se
                   - 8*a_1*l2*r1*r2*se + 8*l1*l2*r1*r2*se) \
                  / (4*(r1 + r2)*(l2**2*r1**2*t - 2*l2**2*r1*se - l1*l2*r1**2*t
                                  + l2**2*r1*r2*t + 2*l1*l2*r1*se - l1*l2*r1*r2*t))

    result = [D0, D1]
    for i in range(m - 1):
        result.append(np.diag([q[0, i], q[1, i]]) * D1)
    result.append(np.diag([1.0 - np.sum(q[0, :]), 1.0 - np.sum(q[1, :])]) * D1)
    return result


def m3pp2m_fitc_theoretical(mmap: List[np.ndarray], method: str = 'approx_delta',
                             t: Optional[float] = None,
                             tinf: Optional[float] = None) -> List[np.ndarray]:
    """
    Fit the theoretical characteristics of a MMAP(n,m) with a M3PP(2,m).

    Faithful port of MATLAB ``m3pp2m_fitc_theoretical.m``: computes the
    joint-process and per-class count moments of ``mmap`` at the chosen time
    scales and dispatches to the corresponding closed-form / approximate fitter.

    Args:
        mmap: the MMAP(n,m) to fit, as ``[D0, D1, D1_1, ..., D1_m]``.
        method: 'exact_delta', 'approx_delta' or 'approx_cov'.
        t: optional finite time scale (if given, t1=t2=t3=t).
        tinf: optional near-infinite time scale.

    Returns:
        Fitted M3PP as ``[D0, D1, D1_1, ..., D1_m]``.
    """
    from line_solver.api.mam.map_analysis import (
        map_count_mean, map_count_var, map_count_moment)
    from line_solver.api.mam.mmap_ops import mmap_count_mean, mmap_count_var

    D0 = np.asarray(mmap[0], dtype=float)
    D1 = np.asarray(mmap[1], dtype=float)
    m = len(mmap) - 2

    if method == 'approx_cov' and m > 2:
        raise ValueError('Approximate covariance fitting only supported for two classes.')

    if t is None:
        t1, t2, t3 = 1.0, 10.0, 10.0
        tinf = 1e4
    else:
        t1 = t2 = t3 = float(t)
        if tinf is None:
            tinf = 1e4

    # joint-process characteristics
    a = float(np.ravel(map_count_mean(D0, D1, t1))[0]) / t1
    bt1 = float(np.ravel(map_count_var(D0, D1, t1))[0]) / (a * t1)
    bt2 = float(np.ravel(map_count_var(D0, D1, t2))[0]) / (a * t2)
    binf = float(np.ravel(map_count_var(D0, D1, tinf))[0]) / (a * tinf)
    mt2 = np.ravel(map_count_moment(D0, D1, t2, np.array([1, 2, 3])))
    m3t2 = mt2[2] - 3 * mt2[1] * mt2[0] + 2 * mt2[0] ** 3

    # per-class rates
    ai = np.ravel(mmap_count_mean(mmap, 1.0))

    if method in ('exact_delta', 'approx_delta'):
        dvt3 = np.zeros(m)
        for i in range(m):
            mmap2 = [mmap[0], mmap[1], mmap[2 + i], D1 - np.asarray(mmap[2 + i], dtype=float)]
            Vt3 = np.ravel(mmap_count_var(mmap2, t3))
            dvt3[i] = Vt3[0] - Vt3[1]

    if method == 'exact_delta':
        return m3pp2m_fitc(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3)
    elif method == 'approx_delta':
        return m3pp2m_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3)
    elif method == 'approx_cov':
        from line_solver.lib.kpctoolbox import mmpp2_fitc_approx
        # count covariance between the two classes at scale t3:
        # Cov(N1,N2) = 1/2*(Var[N1+N2] - Var[N1] - Var[N2])
        var_tot = float(np.ravel(map_count_var(D0, D1, t3))[0])
        var_cls = np.ravel(mmap_count_var(mmap, t3))
        st3 = 0.5 * (var_tot - float(np.sum(var_cls)))
        # fit underlying MMPP(2), then the M3PP(2,2) covariance match
        fit_mmpp = mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2)
        return m3pp22_fitc_approx_cov_multiclass(list(fit_mmpp), ai, st3, t3)
    else:
        raise ValueError("Invalid method '%s'" % method)


def m3pp22_fitc_approx_cov_multiclass(mmpp: List[np.ndarray], ai, st3: float,
                                      t3: float) -> List[np.ndarray]:
    """
    Fit an M3PP(2,2) given the underlying MMPP(2), matching the per-class rates
    ``ai`` and the count covariance ``st3`` between the two classes at scale
    ``t3``. Port of MATLAB m3a m3pp22_fitc_approx_cov_multiclass.

    Args:
        mmpp: underlying MMPP(2) as [D0, D1].
        ai: rates of the two classes (length m, m <= 2).
        st3: count covariance between the two classes at scale t3.
        t3: third time scale.

    Returns:
        Fitted M3PP as [D0, D1, D1_1, ..., D1_m].
    """
    ai = np.ravel(np.asarray(ai, dtype=float))
    m = ai.shape[0]
    if m > 2:
        raise ValueError('No more than two classes supported')

    D0 = np.asarray(mmpp[0], dtype=float)
    D1 = np.asarray(mmpp[1], dtype=float)

    # degenerate case: Poisson process (single phase)
    if D0.shape[0] == 1:
        a = float(np.sum(ai))
        pi = ai / a
        FIT = [D0, D1]
        for i in range(m):
            FIT.append(pi[i] * D1)
        return FIT

    # a single class
    if m == 1:
        return [D0, D1, D1.copy()]

    l1 = D1[0, 0]
    l2 = D1[1, 1]
    r1 = D0[0, 1]
    r2 = D0[1, 0]
    t = t3
    a1 = ai[0]

    e = np.exp(-(r1 + r2) * t)
    w0 = (2 * r1 * (1 - e - (r1 + r2) * t) *
          (a1 ** 2 * (r1 + r2) - a1 * r2 * (l1 - l2))) / (r2 * (r1 + r2) ** 3)
    w1 = -(2 * r1 * (1 - e - (r1 + r2) * t) *
           (2 * a1 * l2 * (r1 + r2) - l2 * r2 * (l1 - l2))) / (r2 * (r1 + r2) ** 3)
    w2 = ((2 * l2 ** 2 * r2 * t) * (r1 + r2) + 2 * l2 ** 2 * r1 * (1 - e)) / \
         (r2 * (r1 + r2) ** 2) - (2 * l2 ** 2 * t) / r2
    w3 = (r1 + r2) / (l2 * r1)
    w4 = (l1 * r2) / (l2 * r1)

    # bounds for the first and the second root
    L1 = -np.inf
    L2 = -np.inf
    U1 = np.inf
    U2 = np.inf

    # if set, the first (second) root is never feasible
    infeasible1 = False
    infeasible2 = False

    # defined for convenience
    z = w0 - w1 ** 2 / (4 * w2)

    # impose square root argument is >= 0
    if w2 > 0:
        L1 = max(L1, z)
        L2 = max(L2, z)
    elif w2 < 0:
        U1 = min(U1, z)
        U2 = min(U2, z)

    # impose q2 >= 0
    if w1 >= 0:
        L1 = max(L1, w0)
    elif w2 < 0:
        infeasible1 = True
    if w1 <= 0:
        U2 = min(U2, w0)
    elif w2 > 0:
        infeasible2 = True

    # impose q1 >= 0
    tmp = 2 * a1 * w3 * w2 + w1
    if tmp >= 0:
        U1 = min(U1, z + tmp ** 2 / (4 * w2))
    elif w2 > 0:
        infeasible1 = True
    if tmp <= 0:
        L2 = max(L2, z + tmp ** 2 / (4 * w2))
    elif w2 < 0:
        infeasible2 = True

    # impose q2 <= 1
    tmp = 2 * w2 + w1
    if tmp >= 0:
        U1 = min(U1, z + tmp ** 2 / (4 * w2))
    elif w2 > 0:
        infeasible1 = True
    if tmp <= 0:
        L2 = max(L2, z + tmp ** 2 / (4 * w2))
    elif w2 < 0:
        infeasible2 = True

    # impose q1 <= 1
    tmp = 2 * a1 * w2 * w3 - 2 * w2 * w4 + w1
    if tmp >= 0:
        L1 = max(L1, z + tmp ** 2 / (4 * w2))
    elif w2 < 0:
        infeasible1 = True
    if tmp <= 0:
        U2 = min(U2, z + tmp ** 2 / (4 * w2))
    elif w2 > 0:
        infeasible2 = True

    if infeasible1 and infeasible2:
        raise ValueError('Empty feasibility region. This should not happen.')

    # compute feasible covariance
    if infeasible2:
        sigma = max(min(st3, U1), L1)
        root = 1
    elif infeasible1:
        sigma = max(min(st3, U2), L2)
        root = 2
    else:
        sigma1 = max(min(st3, U1), L1)
        sigma2 = max(min(st3, U2), L2)
        if abs(sigma1 - st3) < abs(sigma2 - st3):
            sigma = sigma1
            root = 1
        else:
            sigma = sigma2
            root = 2

    # compute parameters
    if root == 1:
        q2 = (-w1 + np.sqrt(w1 ** 2 - 4 * w2 * (w0 - sigma))) / (2 * w2)
    else:
        q2 = (-w1 - np.sqrt(w1 ** 2 - 4 * w2 * (w0 - sigma))) / (2 * w2)
    q1 = (a1 * (r1 + r2) - l2 * q2 * r1) / (l1 * r2)

    # check feasibility just to be sure
    tol = 1e-8
    if (q1 >= tol and q1 <= 1 + tol) and (q2 >= tol and q2 <= 1 + tol):
        q1 = min(max(q1, 0.0), 1.0)
        q2 = min(max(q2, 0.0), 1.0)
    else:
        raise ValueError('Parameters are infeasible. This should not happen.')

    # assemble M3PP[2] (elementwise split of D1 between the two classes)
    Q1 = np.array([[q1, 0.0], [0.0, q2]])
    Q2 = np.array([[1 - q1, 0.0], [0.0, 1 - q2]])
    return [D0, D1, D1 * Q1, D1 * Q2]


def m3pp22_fitc_approx_cov(a: float, bt1: float, bt2: float, binf: float,
                           m3t2: float, t1: float, t2: float,
                           ai: np.ndarray, st3: float, t3: float
                           ) -> List[np.ndarray]:
    """
    Fit a second-order Marked MMPP with two classes, matching the count
    covariance between them. Port of MATLAB
    ``m3a/m3pp/m3pp22_fitc_approx_cov.m``.

    Args:
        a: Total arrival rate
        bt1, bt2, binf: IDC at scales t1, t2 and for t->inf
        m3t2: Third central moment
        t1, t2, t3: Time scales
        ai: Rates of the two classes
        st3: Count covariance between the two classes at scale t3

    Returns:
        M3PP(2,2) as [D0, D1, D1_1, D1_2]
    """
    from line_solver.lib.kpctoolbox import mmpp2_fitc_approx

    ai = np.asarray(ai, dtype=float).ravel()
    if abs(a - np.sum(ai)) > 1e-8:
        raise ValueError('Inconsistent per-class arrival rates.')

    FIT = mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2)
    return m3pp22_fitc_approx_cov_multiclass(list(FIT), ai, st3, t3)


def m3pp2m_fitc_trace(T: np.ndarray, A: np.ndarray, method: str = 'approx_delta',
                      t1: Optional[float] = None, tinf: Optional[float] = None
                      ) -> List[np.ndarray]:
    """
    Fit a multi-class trace with an M3PP(2,m) on its counting-process
    characteristics. Port of MATLAB ``m3a/m3pp/m3pp2m_fitc_trace.m``.

    Args:
        T: Inter-arrival times
        A: Class labels
        method: 'exact_delta', 'approx_delta', 'approx_cov' or 'approx_ag'
        t1: Finite time scale (default 10*mean(T))
        tinf: Near-infinite time scale (default max(10*t1, (sum(T)-T[0])/100))

    Returns:
        M3PP as [D0, D1, D1_1, ..., D1_m]
    """
    from line_solver.api.trace import mtrace_iat2counts

    T = np.asarray(T, dtype=float).ravel()
    A = np.asarray(A).ravel()
    labels = np.unique(A)
    m = len(labels)

    if method == 'approx_cov' and m > 2:
        raise ValueError('Approximate covariance fitting only supported for two classes.')

    TC = np.cumsum(T)
    if t1 is None:
        t1 = 10.0 * float(np.mean(T))
        tinf = max(10.0 * t1, (TC[-1] - TC[0]) / 100.0)
    elif tinf is None:
        tinf = max(10.0 * t1, (TC[-1] - TC[0]) / 100.0)
    t2 = t1 + float(np.mean(T))
    t3 = tinf  # controls the approximation; the reference uses tinf, not t1

    mNt1 = mtrace_iat2counts(T, A, t1)
    mNt2 = mNt1
    mNtinf = mtrace_iat2counts(T, A, tinf)
    mNt3 = mNt1

    Nt1 = np.sum(mNt1, axis=1)
    Nt2 = np.sum(mNt2, axis=1)
    Ntinf = np.sum(mNtinf, axis=1)

    a = 1.0 / float(np.mean(T))
    ai = np.array([a * np.sum(A == labels[i]) / len(A) for i in range(m)])

    bt1 = float(np.var(Nt1)) / (a * t1)
    bt2 = bt1
    binf = float(np.var(Ntinf)) / (a * tinf)
    mt2 = np.array([np.mean(Nt2), np.mean(Nt2 ** 2), np.mean(Nt2 ** 3)])
    m3t2 = mt2[2] - 3 * mt2[1] * mt2[0] + 2 * mt2[0] ** 3

    if method in ('exact_delta', 'approx_delta'):
        dvt3 = np.zeros(m)
        for i in range(m):
            other = np.sum(mNt3, axis=1) - mNt3[:, i]
            dvt3[i] = float(np.var(mNt3[:, i])) - float(np.var(other))

    if method == 'exact_delta':
        return m3pp2m_fitc(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3)
    if method == 'approx_delta':
        return m3pp2m_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3)
    if method == 'approx_cov':
        V = float(np.var(np.sum(mNt3, axis=1)))
        vi = [float(np.var(mNt3[:, 0])), float(np.var(mNt3[:, 1]))]
        s = 0.5 * (V - sum(vi))
        return m3pp22_fitc_approx_cov(a, bt1, bt2, binf, m3t2, t1, t2, ai, s, t3)
    if method == 'approx_ag':
        st3 = np.zeros(m)
        for i in range(m):
            other = np.sum(mNt3, axis=1) - mNt3[:, i]
            st3[i] = float(np.cov(mNt3[:, i], other)[0, 1])
        gt3 = np.array([float(np.var(mNt3[:, i])) for i in range(m)]) + st3
        return m3pp2m_fitc_approx_ag(a, bt1, bt2, binf, m3t2, t1, t2, ai, gt3, t3)
    raise ValueError("Invalid method '%s'" % method)


def m3pp_superpos_fitc(av: np.ndarray, btv: np.ndarray, binfv: np.ndarray,
                       m3tv: np.ndarray, t: float, tinf: float
                       ) -> Tuple[List[np.ndarray], List[List[np.ndarray]]]:
    """
    Fit one second-order M3PP per class and superpose them. Port of MATLAB
    ``m3a/m3pp/m3pp_superpos_fitc.m``.

    Args:
        av: Per-process rates
        btv: Per-process IDC(t)
        binfv: Per-process IDC(inf)
        m3tv: Per-process third central moment of counts at t
        t: Finite time scale
        tinf: Near-infinite time scale

    Returns:
        (superposed M3PP, list of the fitted per-process M3PPs)
    """
    from line_solver.lib.kpctoolbox import mmpp2_fitc
    from line_solver.api.mam import mmap_super

    av = np.asarray(av, dtype=float).ravel()
    btv = np.asarray(btv, dtype=float).ravel()
    binfv = np.asarray(binfv, dtype=float).ravel()
    m3tv = np.asarray(m3tv, dtype=float).ravel()
    m = len(av)

    m3pps = []
    for i in range(m):
        mmpp = mmpp2_fitc(av[i], btv[i], btv[i], binfv[i], m3tv[i], t, tinf)
        D0 = np.atleast_2d(np.asarray(mmpp[0], dtype=float))
        D1 = np.atleast_2d(np.asarray(mmpp[1], dtype=float))
        m3pps.append([D0, D1, D1.copy()])

    fit = m3pps[0]
    for i in range(1, m):
        fit = mmap_super(fit, m3pps[i])
    return fit, m3pps


def m3pp_superpos_fitc_theoretical(MMAP: List[np.ndarray], t: float, tinf: float
                                   ) -> Tuple[List[np.ndarray], List[List[np.ndarray]]]:
    """
    Superpose one M3PP per class to fit the counting characteristics of a
    MMAP(n,m). Port of MATLAB ``m3a/m3pp/m3pp_superpos_fitc_theoretical.m``.

    Args:
        MMAP: Process to fit, as [D0, D1, D1_1, ..., D1_m]
        t: Finite time scale
        tinf: Near-infinite time scale

    Returns:
        (superposed M3PP, list of the fitted per-class M3PPs)
    """
    from line_solver.api.mam import (mmap_count_mean, mmap_count_idc,
                                     mmap_count_moment)

    m = len(MMAP) - 2
    av = np.ravel(mmap_count_mean(MMAP, 1.0))
    btv = np.ravel(mmap_count_idc(MMAP, t))
    binfv = np.ravel(mmap_count_idc(MMAP, tinf))
    mtv = np.asarray(mmap_count_moment(MMAP, t, np.array([1, 2, 3])), dtype=float)
    m3tv = np.array([mtv[2, i] - 3 * mtv[1, i] * mtv[0, i] + 2 * mtv[0, i] ** 3
                     for i in range(m)])
    return m3pp_superpos_fitc(av, btv, binfv, m3tv, t, tinf)


def m3pp_superpos_fitc_trace(T: np.ndarray, A: np.ndarray,
                             t: Optional[float] = None, tinf: Optional[float] = None
                             ) -> Tuple[List[np.ndarray], List[List[np.ndarray]]]:
    """
    Superpose one M3PP per class to fit a multi-class trace. Port of MATLAB
    ``m3a/m3pp/m3pp_superpos_fitc_trace.m``.

    Args:
        T: Inter-arrival times
        A: Class labels
        t: Finite time scale (default 10*mean(T))
        tinf: Near-infinite time scale

    Returns:
        (superposed M3PP, list of the fitted per-class M3PPs)
    """
    from line_solver.api.trace import mtrace_iat2counts

    T = np.asarray(T, dtype=float).ravel()
    A = np.asarray(A).ravel()
    if t is None:
        t = 10.0 * float(np.mean(T))
        tinf = max(10.0 * t, (float(np.sum(T)) - T[0]) / 100.0)
    elif tinf is None:
        tinf = max(10.0 * t, (float(np.sum(T)) - T[0]) / 100.0)

    a = 1.0 / float(np.mean(T))
    labels = np.unique(A)
    m = len(labels)
    pv = np.array([np.sum(A == labels[j]) / len(A) for j in range(m)])
    av = pv * a

    Nt = mtrace_iat2counts(T, A, t)
    Ninf = mtrace_iat2counts(T, A, tinf)

    btv = np.zeros(m)
    binfv = np.zeros(m)
    m3tv = np.zeros(m)
    for i in range(m):
        btv[i] = float(np.var(Nt[:, i])) / (av[i] * t)
        binfv[i] = float(np.var(Ninf[:, i])) / (av[i] * tinf)
        mt = np.array([np.mean(Nt[:, i]), np.mean(Nt[:, i] ** 2), np.mean(Nt[:, i] ** 3)])
        m3tv[i] = mt[2] - 3 * mt[1] * mt[0] + 2 * mt[0] ** 3

    return m3pp_superpos_fitc(av, btv, binfv, m3tv, t, tinf)


def _m3pp_compute_d(bt1: float, binf: float, t1: float) -> float:
    """
    Solve d = r1 + r2 for an MMPP(2) matching IDC(t1) = bt1 and IDC(inf) = binf.

    The reference solves z = w exp(w) with fsolve from w = 1; the closed form is
    the principal Lambert W branch, since z = -c exp(-c) with c > 1 has the
    other root at w = -c on the W_{-1} branch.
    """
    from scipy.special import lambertw

    if not (binf > bt1 and bt1 > 1):
        raise ValueError('No solution, infeasible IDC(t): IDC(%.2f) = %.3f, IDC(inf) = %.3f'
                         % (t1, bt1, binf))
    c = (binf - 1.0) / (binf - bt1)
    z = -c * np.exp(-c)
    w = float(np.real(lambertw(z, 0)))
    return (w + c) / t1


def m3pp22_interleave_fitc(av: np.ndarray, btv: np.ndarray, binfv: np.ndarray,
                           stv: np.ndarray, t: float
                           ) -> Tuple[List[np.ndarray], List[List[np.ndarray]]]:
    """
    Fit L pairs of classes into a single MMAP obtained by lumped superposition
    of L M3PP(2,2) processes. Port of MATLAB
    ``m3a/m3pp/m3pp22_interleave_fitc.m``.

    The off-diagonal rates come from a FEASIBILITY linear program: the reference
    passes a zero objective, so every feasible point is optimal and which vertex
    is returned is the LP solver's choice. Different solvers therefore hand
    different MMPP(2)s to the per-pair covariance split, and a covariance that
    one accepts another can report infeasible. That is a property of the
    reference, not of this port.

    Args:
        av: (L,2) per-class rates
        btv: (L,) IDC at resolution t of each pair
        binfv: (L,) asymptotic IDC of each pair
        stv: (L,) count covariance at resolution t within each pair
        t: Time scale

    Returns:
        (lumped superposition, list of the L M3PP(2,2) processes)
    """
    from scipy.optimize import linprog

    av = np.atleast_2d(np.asarray(av, dtype=float))
    btv = np.asarray(btv, dtype=float).ravel()
    binfv = np.asarray(binfv, dtype=float).ravel()
    stv = np.asarray(stv, dtype=float).ravel()
    L = av.shape[0]

    # bounds on the upper off-diagonal element of each MMPP(2)
    uv = np.zeros(L)
    dv = np.zeros(L)
    for i in range(L):
        a = float(np.sum(av[i, :]))
        d = _m3pp_compute_d(btv[i], binfv[i], t)
        z = (binfv[i] - 1.0) * d ** 3 * a
        uv[i] = d * z / (2 * a ** 2 * d ** 2 + z)
        dv[i] = d

    A_ub = np.zeros((2 * L, 2 * L))
    b_ub = np.zeros(2 * L)
    for i in range(L):
        base = 2 * i
        for j in range(L):
            if j >= i:
                A_ub[base, j] = 1.0
                A_ub[base + 1, j] = -1.0
        b_ub[base] = dv[i] - 1e-6
        b_ub[base + 1] = -uv[i] - 1e-6

    A_eq = np.zeros((L, 2 * L))
    b_eq = np.zeros(L)
    for i in range(L):
        for j in range(L):
            if j >= i:
                A_eq[i, j] = 1.0
            if j <= i:
                A_eq[i, j + L] = 1.0
        b_eq[i] = dv[i]

    sol = linprog(np.zeros(2 * L), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=[(0.0, None)] * (2 * L), method='highs')
    if not sol.success:
        raise ValueError('m3pp22_interleave_fitc: no feasible off-diagonal rates')
    x = np.asarray(sol.x, dtype=float).ravel()
    r = np.vstack([x[:L], x[L:(2 * L)]])

    m3pps = []
    for i in range(L):
        r1 = float(np.sum(r[0, i:]))
        r2 = float(np.sum(r[1, :(i + 1)]))
        a = float(np.sum(av[i, :]))
        d = r1 + r2
        z = (binfv[i] - 1.0) * d ** 3 * a
        delta = np.sqrt(z / (2 * r1 * r2))
        l2 = a - r2 / d * delta
        l1 = l2 + delta
        D0 = np.array([[0.0, r1], [r2, 0.0]])
        D1 = np.array([[l1, 0.0], [0.0, l2]])
        for j in range(2):
            D0[j, j] = -np.sum(D0[j, :] + D1[j, :])
        m3pps.append(m3pp22_fitc_approx_cov_multiclass([D0, D1], av[i, :], stv[i], t))

    return m3pp2m_interleave(m3pps), m3pps


def m3pp_superpos(m3pps: List[List[np.ndarray]]) -> List[np.ndarray]:
    """
    Compute the superposition of multiple M3PP processes.

    Args:
        m3pps: List of M3PP processes

    Returns:
        Superposed M3PP
    """
    from line_solver.api.mam import mmap_super

    if len(m3pps) == 0:
        raise ValueError("Empty M3PP list")
    if len(m3pps) == 1:
        return m3pps[0]

    result = m3pps[0]
    for i in range(1, len(m3pps)):
        result = mmap_super(result, m3pps[i])

    return result


__all__ = [
    'm3pp_rand',
    'm3pp2m_interleave',
    'm3pp2m_fitc_approx_ag_multiclass',
    'm3pp2m_fitc_approx_ag',
    'm3pp2m_fitc_approx',
    'm3pp2m_fitc',
    'm3pp2m_fitc_theoretical',
    'm3pp2m_fitc_trace',
    'm3pp22_fitc_approx_cov_multiclass',
    'm3pp22_fitc_approx_cov',
    'm3pp22_interleave_fitc',
    'm3pp_superpos',
    'm3pp_superpos_fitc',
    'm3pp_superpos_fitc_theoretical',
    'm3pp_superpos_fitc_trace',
]
