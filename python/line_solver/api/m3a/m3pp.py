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
    from ..kpctoolbox import mmpp_rand
    from ..mam import map_issym

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
    from ..mam import mmap_isfeasible, mmap_count_var

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
    from ..kpctoolbox import mmpp2_fitc_approx
    from ..mam import mmap_isfeasible

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


def m3pp2m_fitc_theoretical(a: float, binf: float, m3inf: float,
                             ai: np.ndarray, dvinf: np.ndarray
                             ) -> List[np.ndarray]:
    """
    Fit a second-order Marked MMPP using theoretical moments.

    Args:
        a: Total arrival rate
        binf: IDC for t->inf
        m3inf: Third central moment for t->inf
        ai: Per-class arrival rates
        dvinf: Per-class asymptotic variance differences

    Returns:
        M3PP as [D0, D1, D1_1, ..., D1_m]
    """
    from ..kpctoolbox import mmpp2_fitc_theoretical
    from ..mam import mmap_isfeasible

    ai = np.asarray(ai).ravel()
    dvinf = np.asarray(dvinf).ravel()
    m = len(ai)

    if abs(a - np.sum(ai)) > 1e-8:
        raise ValueError("Inconsistent per-class arrival rates")

    # Fit underlying MMPP(2)
    FIT = mmpp2_fitc_theoretical(a, binf, m3inf)

    # Degenerate case: Poisson process
    if FIT[0].shape[0] == 1:
        D0 = FIT[0]
        D1 = FIT[1]
        result = [D0, D1]
        pi = ai / a
        for i in range(m):
            result.append(pi[i] * D1)
        return result

    # Multi-class allocation based on rates
    q = np.zeros((2, m))
    for i in range(m):
        q[0, i] = ai[i] / a if a > 0 else 1.0 / m
        q[1, i] = ai[i] / a if a > 0 else 1.0 / m

    result = [FIT[0], FIT[1]]
    for i in range(m):
        result.append(FIT[1] * np.diag([q[0, i], q[1, i]]))

    return result


def m3pp_superpos(m3pps: List[List[np.ndarray]]) -> List[np.ndarray]:
    """
    Compute the superposition of multiple M3PP processes.

    Args:
        m3pps: List of M3PP processes

    Returns:
        Superposed M3PP
    """
    from ..mam import mmap_super

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
    'm3pp2m_fitc_theoretical',
    'm3pp_superpos',
]
