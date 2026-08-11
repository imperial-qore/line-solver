"""
APH(2) moment matching.

Port of the M3A reference toolbox (matlab/lib/m3a/m3a/aph2/), used by the
order-1 MMAP mixture compression to build one two-phase component per class.

Functions
---------
- aph2_assemble: build the APH(2) [D0, D1] from its parameters
- aph2_fitall:   all feasible APH(2) fits of the first three moments
- aph2_adjust:   nearest APH(2)-feasible (M2, M3) for a given M1
- aph2_fit:      fit an APH(2), adjusting the moments if infeasible

Reference
---------
A. Bobbio, A. Horvath, M. Telek. Matching three moments with minimal acyclic
phase type distributions. Stochastic Models 21:303-326, 2005.
M. Telek, A. Heindl. Matching moments for acyclic discrete and continuous
phase-type distributions of second order. 2002.
"""

from typing import List, Tuple

import numpy as np

from ...lib.kpctoolbox.aph import aph_fit

__all__ = ['aph2_assemble', 'aph2_fitall', 'aph2_adjust', 'aph2_fit']


def aph2_assemble(l1: float, l2: float, p1: float) -> List[np.ndarray]:
    """
    Build the APH(2) with the given parameters. Mirrors aph2_assemble.m.

    Args:
        l1: mean of the first phase (the MATLAB source parameterizes by 1/rate)
        l2: mean of the second phase
        p1: probability of continuing to the second phase

    Returns:
        The APH(2) as [D0, D1].
    """
    D0 = np.array([[-1.0 / l1, (1.0 / l1) * p1],
                   [0.0, -1.0 / l2]])
    D1 = np.array([[(1.0 / l1) * (1.0 - p1), 0.0],
                   [1.0 / l2, 0.0]])
    return [D0, D1]


def aph2_fitall(M1: float, M2: float, M3: float) -> List[List[np.ndarray]]:
    """
    Fit every APH(2) matching the given first three moments. Mirrors
    aph2_fitall.m.

    Returns a list of 1 or 2 APH(2) fits. It is never empty: when no APH(2)
    solution is feasible the general APH fitter aph_fit(M1, M2, M3, 2) supplies
    a fallback, exactly as in the MATLAB and JAR sources.
    """
    degentol = 1e-8

    SCV = (M2 - M1 ** 2) / M1 ** 2

    # lower bound of M3, defined only for SCV <= 1 (for SCV > 1 the (1-SCV)^(3/2)
    # factor is complex; the MATLAB source relies on && short-circuiting past it).
    if SCV <= 1:
        M3lb = 3 * M1 ** 3 * (3 * SCV - 1 + np.sqrt(2.0) * (1 - SCV) ** 1.5)
        at_lower_bound = abs(M3 - M3lb) < degentol
    else:
        at_lower_bound = False

    if at_lower_bound:
        # at the lower bound of M3 the square-root argument vanishes
        tmp0 = 0.0
    else:
        tmp0 = (M3 ** 2 / 9.0 + ((8 * M1 ** 3) / 3.0 - 2 * M2 * M1) * M3
                - 3 * M1 ** 2 * M2 ** 2 + 2 * M2 ** 3)
        if tmp0 < 0:
            # infeasible: square root of a negative number
            fitted, _ = aph_fit(M1, M2, M3, 2)
            return [[fitted['D0'], fitted['D1']]]

    tmp1 = 3 * np.sqrt(tmp0)
    tmp2 = M3 - 3 * M1 * M2
    tmp3 = 6 * M2 - 12 * M1 ** 2

    # maximum number of solutions: when tmp0 == 0 the diagonal of D0 is uniform
    n = 1 if tmp0 == 0.0 else 2

    h1v = np.zeros(n)
    h2v = np.zeros(n)
    if n == 1:
        h2v[0] = tmp2 / tmp3
        h1v[0] = h2v[0]
    else:
        h2v[0] = (tmp2 + tmp1) / tmp3
        h2v[1] = (tmp2 - tmp1) / tmp3
        h1v[1] = h2v[0]
        h1v[0] = h2v[1]

    APHS = []
    for j in range(n):
        h1 = h1v[j]
        h2 = h2v[j]
        r1 = (M1 - h1) / h2
        if h1 > 0 and h2 > 0 and -degentol <= r1 <= 1 + degentol:
            # feasible (or almost feasible) solution
            r1 = max(min(r1, 1.0), 0.0)
            APHS.append(aph2_assemble(h1, h2, r1))

    if not APHS:
        fitted, _ = aph_fit(M1, M2, M3, 2)
        APHS = [[fitted['D0'], fitted['D1']]]

    return APHS


def aph2_adjust(M1: float, M2: float, M3: float,
                method: str = 'simple') -> Tuple[float, float]:
    """
    Find the APH(2)-feasible (M2, M3) closest to the given pair, holding M1
    fixed. Mirrors the 'simple' branch of aph2_adjust.m [Telek and Heindl,
    2002], which is the only branch aph2_fit uses.

    Returns:
        (M2a, M3a), the adjusted moments. M1 is never altered.
    """
    if method != 'simple':
        raise ValueError("Invalid method: %s" % method)

    tol = 1e-4  # tolerance for strict inequalities
    M1sq = M1 ** 2
    scv = (M2 - M1sq) / M1sq

    if scv < 0.5:
        # clamp to SCV = 1/2, the APH(2) lower bound
        M2a = 1.5 * M1 ** 2
        scva = (M2a - M1sq) / M1sq
    else:
        M2a = M2
        scva = scv

    if scva <= 1:
        lb = 3 * M1 ** 3 * (3 * scva - 1 + np.sqrt(2.0) * (1 - scva) ** 1.5)
        ub = 6 * M1 ** 3 * scva
        if M3 < lb:
            M3a = lb
        elif M3 > ub:
            M3a = ub
        else:
            M3a = M3
    else:
        lb = 1.5 * M1 ** 3 * (1 + scva) ** 2
        if M3 <= lb:
            M3a = lb * (1 + tol)
        else:
            M3a = M3

    return M2a, M3a


def aph2_fit(M1: float, M2: float, M3: float) -> List[np.ndarray]:
    """
    Fit an APH(2) to the first three moments, adjusting M2/M3 if they are
    infeasible. Mirrors aph2_fit.m.

    M1 is always matched exactly: aph2_adjust never alters it, and the aph_fit
    fallback rescales its result to mean M1.

    Args:
        M1, M2, M3: the first three raw moments.

    Returns:
        The fitted APH(2) as [D0, D1].
    """
    APHS = aph2_fitall(M1, M2, M3)

    if not APHS:
        # Unreachable with the current aph2_fitall, which always falls back to
        # aph_fit. Kept to mirror the MATLAB and JAR control flow.
        M2a, M3a = aph2_adjust(M1, M2, M3)
        APHS = aph2_fitall(M1, M2a, M3a)
        if not APHS:
            raise RuntimeError('Fitting APH(2): feasibility could not be restored')

    return APHS[0]
