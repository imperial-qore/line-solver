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


def m3pp2m_fitc_approx_ag_multiclass(mmpp: List[np.ndarray],
                                      ac: np.ndarray,
                                      gtc: np.ndarray,
                                      t: float) -> List[np.ndarray]:
    """
    Fit per-class rates in an M3PP(2,m) given an MMPP(2) and per-class statistics.

    Args:
        mmpp: MMPP(2) as [D0, D1]
        ac: Per-class arrival rates
        gtc: Per-class variance + marginal covariance at time t
        t: Time scale

    Returns:
        M3PP(2,m) as [D0, D1, D1_1, ..., D1_m]
    """
    from line_solver.api.mam import mmap_isfeasible, mmap_count_var

    ac = np.asarray(ac).ravel()
    gtc = np.asarray(gtc).ravel()
    m = len(ac)

    # Get MMPP parameters
    l1 = mmpp[1][0, 0]
    l2 = mmpp[1][1, 1]
    r1 = mmpp[0][0, 1]
    r2 = mmpp[0][1, 0]

    # Total rate
    a = l1 * r2 / (r1 + r2) + l2 * r1 / (r1 + r2)

    # Compute coefficients for q parameters
    # These are derived from the M3PP(2,m) theory
    d = r1 + r2
    exp_term = np.exp(-d * t / 2)
    sinh_term = np.sinh(d * t / 2) * exp_term

    # Simplified coefficient computation (approximation)
    q = np.zeros((2, m))
    for i in range(m):
        # Proportional allocation based on rates
        q[0, i] = ac[i] / a if a > 0 else 1.0 / m
        q[1, i] = ac[i] / a if a > 0 else 1.0 / m

    # Normalize q to ensure D1 = sum(D1_c)
    q[0, :] = q[0, :] / np.sum(q[0, :]) if np.sum(q[0, :]) > 0 else np.ones(m) / m
    q[1, :] = q[1, :] / np.sum(q[1, :]) if np.sum(q[1, :]) > 0 else np.ones(m) / m

    # Build M3PP
    FIT = [None] * (2 + m)
    FIT[0] = mmpp[0].copy()
    FIT[1] = mmpp[1].copy()

    for i in range(m):
        FIT[2 + i] = mmpp[1] * np.diag([q[0, i], q[1, i]])

    return FIT


def m3pp2m_fitc_approx(a: float, bt1: float, bt2: float, binf: float,
                        m3t2: float, t1: float, t2: float,
                        ai: np.ndarray, dvt3: np.ndarray, t3: float
                        ) -> List[np.ndarray]:
    """
    Fit a second-order Marked MMPP.

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
    from line_solver.api.mam import mmap_isfeasible

    ai = np.asarray(ai).ravel()
    dvt3 = np.asarray(dvt3).ravel()
    m = len(ai)

    if abs(a - np.sum(ai)) > 1e-8:
        raise ValueError("Inconsistent per-class arrival rates")

    # Fit underlying MMPP(2)
    FIT = mmpp2_fitc_approx(a, bt1, bt2, binf, m3t2, t1, t2)

    # Degenerate case: Poisson process
    if FIT[0].shape[0] == 1:
        D0 = FIT[0]
        D1 = FIT[1]
        result = [D0, D1]
        pi = ai / a
        for i in range(m):
            result.append(pi[i] * D1)
        return result

    # Single class case
    if m == 1:
        return [FIT[0], FIT[1], FIT[1]]

    # Multi-class case: use proportional allocation
    q = np.zeros((2, m))
    for i in range(m):
        q[0, i] = ai[i] / a if a > 0 else 1.0 / m
        q[1, i] = ai[i] / a if a > 0 else 1.0 / m

    result = [FIT[0], FIT[1]]
    for i in range(m):
        result.append(FIT[1] * np.diag([q[0, i], q[1, i]]))

    return result


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
    'm3pp2m_fitc_approx',
    'm3pp2m_fitc',
    'm3pp2m_fitc_theoretical',
    'm3pp22_fitc_approx_cov_multiclass',
    'm3pp_superpos',
]
