"""MAMAP(2,2) fitting on a moment plus the class TRANSITION probabilities.

Port of the m3a ``mamap22`` family: ``mamap2m_can1_coefficients.m`` and
``mamap2m_can2_coefficients.m``, ``mamap22_fit_fs_multiclass.m`` and
``mamap22_fit_bs_multiclass.m``, and the six drivers
``mamap22_fit_gamma_{fs,bs}{,_mmap,_trace}.m``. Cross-checked against
``cpp/include/line/api/mam/mamap22_fit_{fs,bs}.h`` and
``jline/api/mam/Mamap2m_coefficients.java``. An MMAP is a list
``[D0, D1, D11, ..., D1m]``.

TWO CLASSES ONLY, by construction: sigma is a class-to-class transition
probability and the closed forms below invert its (1,1) entry alone. The
reference refuses otherwise and so does this.

Unlike the F+B fitter of :mod:`mamap2m`, which is affine in the marking and
therefore a quadratic program, the moment-plus-sigma relation is a RATIO: one of
the three marking probabilities is a quadratic in the moment over a linear
denominator, and the other two follow linearly from it. That inverse is closed
form and EXACT whenever its answer lands in the unit box.

WHAT IS NOT PORTED, and is refused by name rather than approximated. When the
closed form lands outside the box and the caller asked for ``adjust``, the
reference repairs it by solving a NONCONVEX program -- the marking equality is
bilinear in the unknowns -- with YALMIP's ``bmibnb``, a spatial branch-and-bound
that returns a GLOBAL optimum, over two sides (the moment below and above the
mean) and keeps the better. A local method would report a fit the reference
would not have chosen, so ``adjust=False`` (take the clamped closed form) is
offered instead. Every other branch IS ported, including all four degeneracies
and the one-variable quadratic repairs.
"""
import warnings
from typing import List, Optional, Sequence, Tuple

import numpy as np

from .amap2 import amap2_fit_gamma, amap2_fitall_gamma, _map_repair
from .maph2m import maph2m_fit_multiclass

DEGENTOL = 1e-6
DENUMTOL = 1e-12


class Mamap22Unsupported(NotImplementedError):
    """The branch the reference resolves with YALMIP's nonconvex bmibnb."""


# ---------------------------------------------------------------------------
# coefficient tables
# ---------------------------------------------------------------------------

def mamap2m_can1_coefficients(h1: float, h2: float, r1: float, r2: float
                              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Marking coefficients of the FIRST canonical AMAP(2) form (gamma > 0).

    G holds the per-flow contributions to the class probability, the sigma and
    the moments; U the quadratic terms; Y three determinants of G that decide
    whether the system is solvable. Indices here are 0-based, so ``G[0]`` is the
    reference's ``G(1)``.

    REFERENCE TYPO, corrected here as in the JAR and in C++.
    ``mamap2m_can1_coefficients.m`` assigns ``G(10)`` twice in consecutive
    lines, so ``(r1*r2^2)/(r1*r2-r2+1)`` is discarded and ``G(9)`` is left at
    zero. It belongs at ``G(9)``, which is the only reading under which the
    table has no hole. Y does not read it, so the three determinants are
    unaffected.
    """
    A = r2 * (r1 - 1.0) + 1.0
    B = r1 * r2 - r2 + 1.0
    if A == 0.0 or B == 0.0:
        raise ValueError('mamap2m_can1_coefficients: the canonical denominators vanish at '
                         'this (r1, r2); the marking coefficients are undefined there')

    G = np.zeros(15)
    U = np.zeros(12)
    Y = np.zeros(3)

    G[0] = 1.0 - r1 / A
    G[1] = -(r1 * (r2 - 1.0)) / B
    G[2] = (r1 * r2) / B
    G[3] = (r1 * (r1 - 1.0)) / A - r1 + 1.0
    G[4] = -(r1 * (r1 - 1.0) * (r2 - 1.0) * (r2 - 2.0)) / B
    G[5] = (r1 * r2 * (r1 - 1.0) * (r2 - 1.0)) / B
    G[6] = (r1 ** 2 * (r2 - 1.0) ** 2) / A
    G[7] = -(r1 * r2 * (r1 + 1.0) * (r2 - 1.0)) / B
    G[8] = (r1 * r2 ** 2) / B                       # the reference's discarded G(10)
    G[9] = h1 - (h1 * r1) / A
    G[10] = -(r1 * (r2 - 1.0) * (h1 + h2 - h1 * r2)) / B
    G[11] = (r1 * r2 * (h1 + h2 - h1 * r2)) / B
    G[12] = ((h1 + h2 * r1) * (r1 - 1.0) * (r2 - 1.0)) / B
    G[13] = -(r1 * (h1 + h2 * r1) * (r2 - 1.0)) / B
    G[14] = (h2 * r1 * r2) / B

    t1 = h1 - h2 + h2 * r1
    U[0] = B ** 2
    U[1] = -B * (2.0 * h1 - h1 * r1 - 2.0 * h1 * r2 + 3.0 * h2 * r1 - h2 * r1 ** 2
                 + h2 * r1 ** 2 * r2 + h1 * r1 * r2 - h2 * r1 * r2)
    U[2] = r1 * (r2 - 1.0) * t1 ** 2
    U[3] = B * (h2 ** 2 * r1 - h1 ** 2 * r2 + h1 ** 2 + h1 * h2 * r1 - h1 * h2 * r1 * r2)
    U[4] = -r1 * (r2 - 1.0) * B * t1
    U[5] = r1 * (r2 - 1.0) * (h1 - h1 * r2 + h2 * r1) * t1
    t2 = h2 - h1 * r2
    U[6] = B ** 2
    U[7] = -B * (2.0 * h1 - 2.0 * h1 * r2 + h2 * r1 - h1 * r1 * r2 ** 2
                 + h1 * r1 * r2 + h2 * r1 * r2)
    U[8] = r1 * t2 ** 2 * (r2 - 1.0)
    U[9] = B * (h2 ** 2 * r1 - h1 ** 2 * r2 + h1 ** 2 + h1 * h2 * r1 - h1 * h2 * r1 * r2)
    U[10] = -r1 * t2 * (r2 - 1.0) * B
    U[11] = r1 * t2 * (r2 - 1.0) * (h1 - h1 * r2 + h2 * r1)

    Y[0] = (G[0] * G[10] * G[14] - G[0] * G[11] * G[13] - G[1] * G[9] * G[14]
            + G[1] * G[11] * G[12] + G[2] * G[9] * G[13] - G[2] * G[10] * G[12])
    Y[1] = G[2] * G[12] - G[0] * G[14]
    Y[2] = G[9] * G[2] - G[11] * G[0]
    return G, U, Y


def mamap2m_can2_coefficients(h1: float, h2: float, r1: float, r2: float
                              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Marking coefficients of the SECOND canonical AMAP(2) form (gamma < 0).

    The reference names the three tables E, V and Z; the roles are those of G, U
    and Y in :func:`mamap2m_can1_coefficients`.
    """
    A = r2 * (r1 - 1.0) - r1 + 2.0
    B = r1 * (r2 - 1.0) - r2 + 2.0
    C = r1 + r2 - r1 * r2 - 2.0
    if A == 0.0 or B == 0.0 or C == 0.0:
        raise ValueError('mamap2m_can2_coefficients: the canonical denominators vanish at '
                         'this (r1, r2); the marking coefficients are undefined there')

    E = np.zeros(14)
    V = np.zeros(12)
    Z = np.zeros(3)

    E[0] = 1.0 - 1.0 / A
    E[1] = -(r2 - 1.0) / B
    E[2] = r2 / B
    E[3] = (r2 - 2.0) / B - r2 + 2.0
    E[4] = r2 - r2 / B
    E[5] = -(r1 * (r2 - 1.0) ** 2) / C
    E[6] = -r2 - (r2 * (2.0 * r2 - 3.0)) / B
    E[7] = r2 ** 2 / A
    E[8] = h1 - h1 / A
    E[9] = h1 * (r2 - 1.0) - ((r2 - 1.0) * (2.0 * h1 + h2 - h1 * r2)) / B
    E[10] = (r2 * (2.0 * h1 + h2 - h1 * r2)) / B - h1 * r2
    E[11] = h2 - h2 / A
    E[12] = ((h1 + h2 * r1) * (r2 - 1.0)) / C
    E[13] = (h2 * r2) / B

    t1 = h1 - h2 + h2 * r1
    t2 = h1 - h2 - h1 * r1 + h1 * r1 * r2
    V[0] = -C ** 2
    V[1] = -C * (2.0 * h1 + 2.0 * h2 - h1 * r2 - h2 * r2 + h2 * r1 * r2)
    V[2] = h2 * (2.0 * h1 - h1 * r2 + h2 * r1) * C
    V[3] = (r2 - 1.0) * t1 ** 2
    V[4] = t1 * (2.0 * r2 - r1 * r2 + r1 * r2 ** 2 - r2 ** 2)
    V[5] = -(h1 * r2 + h2 * r2 - h1 * r2 ** 2) * t1
    V[6] = -C ** 2
    V[7] = -C * (2.0 * h1 + 2.0 * h2 - h1 * r2 - h2 * r2 + h1 * r1 * r2 ** 2 - h1 * r1 * r2)
    V[8] = h1 * C * (2.0 * h2 + h1 * r1 - h2 * r2 + h1 * r1 * r2 ** 2 - 2.0 * h1 * r1 * r2)
    V[9] = (r2 - 1.0) * t2 ** 2
    V[10] = -r2 * t2 * C
    V[11] = -r2 * (h1 + h2 - h1 * r2) * t2

    Z[0] = (E[9] * E[11] * E[2] - E[9] * E[13] * E[0] - E[10] * E[11] * E[1]
            + E[10] * E[12] * E[0] - E[12] * E[2] * E[8] + E[13] * E[1] * E[8])
    Z[1] = E[11] * E[1] - E[12] * E[0]
    Z[2] = E[9] * E[0] - E[1] * E[8]
    return E, V, Z


# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------

def _canonical_form(D0: np.ndarray, D1: np.ndarray, who: str) -> int:
    """Which of the two canonical acyclic forms this AMAP(2) is in."""
    if D0.shape[0] != 2:
        raise ValueError('%s: Underlying MAP must be of second-order.' % who)
    if D0[1, 0] != 0:
        raise ValueError('%s: Underlying MAP must be acyclic' % who)
    if D1[0, 1] == 0:
        return 1
    if D1[0, 0] == 0:
        return 2
    raise ValueError('%s: Underlying MAP must be in canonical acyclic form' % who)


def _rebuild(h1: float, h2: float, r1: float, r2: float, form: int) -> List[np.ndarray]:
    """The AMAP(2) of the perturbed (h1, h2, r1, r2), in the given form."""
    from line_solver.api.mam.map_analysis import map_normalize
    D0 = np.array([[-1.0 / h1, r1 / h1], [0.0, -1.0 / h2]])
    if form == 1:
        D1 = np.array([[(1.0 - r1) / h1, 0.0], [r2 / h2, (1.0 - r2) / h2]])
    else:
        D1 = np.array([[0.0, (1.0 - r1) / h1], [r2 / h2, (1.0 - r2) / h2]])
    return list(map_normalize(D0, D1))


def _mark(D1: np.ndarray, form: int, q1: float, q2: float, q3: float) -> List[np.ndarray]:
    """The two marked matrices D11, D12 of a marking (q1, q2, q3).

    The marking is clamped into [0, 1] first, as the reference's `fix` does: it
    is admitted within `feastol` of the box and must be a probability here.
    """
    q1 = min(max(q1, 0.0), 1.0)
    q2 = min(max(q2, 0.0), 1.0)
    q3 = min(max(q3, 0.0), 1.0)
    if form == 1:
        m1 = np.array([[q1, 0.0], [q2, q3]])
        m2 = np.array([[1.0 - q1, 0.0], [1.0 - q2, 1.0 - q3]])
    else:
        m1 = np.array([[0.0, q1], [q2, q3]])
        m2 = np.array([[0.0, 1.0 - q1], [1.0 - q2, 1.0 - q3]])
    return [D1 * m1, D1 * m2]


def _project_interval(coefs: Sequence[float], offs: Sequence[float], target: float,
                      who: str) -> float:
    """Minimize (x/target - 1)^2 over {x in [1e-6, 1e6]: 0 <= c_i x + o_i <= 1}.

    The reference solves this with an interior-point QP over one variable, whose
    objective ``0.5 x' H x + h' x`` with ``H = 2/target^2`` and ``h = -2/target``
    is ``(x/target - 1)^2 - 1``. On an interval the minimizer is the projection
    of ``target``, so this is the same optimum, not a substitute for it.
    """
    lo, hi = 1e-6, 1e6
    for c, o in zip(coefs, offs):
        if abs(c) < DENUMTOL:
            continue
        a, b = -o / c, (1.0 - o) / c
        lo = max(lo, min(a, b))
        hi = min(hi, max(a, b))
    if not lo <= hi:
        raise ValueError('%s: the one-variable repair has an empty feasible interval '
                         'for these targets' % who)
    return min(max(float(target), lo), hi)


def _fwd(mmap: Sequence[np.ndarray]) -> np.ndarray:
    from line_solver.api.mam.mmap_ops import mmap_forward_moment
    return np.asarray(mmap_forward_moment(list(mmap), [1]), dtype=float).ravel()


def _bwd(mmap: Sequence[np.ndarray]) -> np.ndarray:
    from line_solver.api.mam.mmap_ops import mmap_backward_moment
    return np.asarray(mmap_backward_moment(list(mmap), [1]), dtype=float).ravel()


def _sig(mmap: Sequence[np.ndarray]) -> np.ndarray:
    from line_solver.api.mam.mmap_ops import mmap_sigma
    return np.atleast_2d(np.asarray(mmap_sigma(list(mmap)), dtype=float))


# ---------------------------------------------------------------------------
# forward moment + sigma
# ---------------------------------------------------------------------------

def mamap22_fit_fs_multiclass(map_: Sequence[np.ndarray], p: Sequence[float],
                              F: Sequence[float], S,
                              classWeights: Optional[Sequence[float]] = None,
                              fsWeights: Optional[Sequence[float]] = None,
                              adjust: bool = True
                              ) -> Tuple[List[np.ndarray], np.ndarray, np.ndarray, bool]:
    """Mark a canonical AMAP(2) on its FORWARD moment and its sigma.

    The class probabilities are matched exactly; the forward moments and the
    one-step class transition probabilities as closely as the form allows.

    Args:
        map_: the underlying second-order AMAP as [D0, D1]
        p: the two class probabilities
        F: the two target forward moments
        S: the target class transition matrix; only S[0, 0] enters
        classWeights: per-class weights in the objective (default: unit)
        fsWeights: (forward, sigma) weights (default: unit)
        adjust: repair an infeasible closed form; the repair is unported

    Returns:
        (mmap, fF, fS, exact)
    """
    who = 'Fitting MAMAP(2,2) F+S'
    D0 = np.asarray(map_[0], dtype=float)
    D1 = np.asarray(map_[1], dtype=float)
    form = _canonical_form(D0, D1, who)

    p = np.asarray(p, dtype=float).ravel()
    F = np.asarray(F, dtype=float).ravel()
    S = np.atleast_2d(np.asarray(S, dtype=float))
    if p.size != 2:
        raise ValueError('%s: only two classes supported' % who)
    if F.size != 2:
        raise ValueError('%s: one forward moment per class is required' % who)
    cw = np.ones(2) if classWeights is None else np.asarray(classWeights, float).ravel()
    fw = np.ones(2) if fsWeights is None else np.asarray(fsWeights, float).ravel()

    feastol = 1e-3
    h1 = -1.0 / D0[0, 0]
    h2 = -1.0 / D0[1, 1]
    r1 = D0[0, 1] * h1
    r2 = D1[1, 1] * h2

    def feasible(*qs) -> bool:
        return all(-feastol <= q <= 1.0 + feastol for q in qs)

    # ---- the Poisson perturbation, which keeps the two-state structure ----
    if ((form == 1 and (r1 < DEGENTOL or r2 > 1 - DEGENTOL
                        or abs(h1 - h2 + h2 * r1) < DEGENTOL))
            or (form == 2 and (r2 > 1 - DEGENTOL or abs(h1 - h2 + h2 * r1) < DEGENTOL))):
        if r1 < DEGENTOL:
            r1 = DEGENTOL
        if r2 > 1 - DEGENTOL:
            r2 = 1 - DEGENTOL
        if abs(h1 - h2 + h2 * r1) < DEGENTOL:
            h1 = h2 * (1 - r1) + DEGENTOL
        D0, D1 = _rebuild(h1, h2, r1, r2, form)

    exact = False
    q1 = q2 = q3 = float('nan')

    if form == 2 and r2 < DEGENTOL and abs(1 - r1) < DEGENTOL:
        # DEGENERATE PHASE-TYPE: one degree of freedom, match p only
        q1 = q2 = q3 = p[0]
    elif form == 1 and r2 < DEGENTOL:
        # CANONICAL PHASE-TYPE. The reference delegates to the BACKWARD fitter
        # and has no forward targets for it, so it uses the ordinary mean for
        # both and warns; the warning is raised rather than dropped.
        from line_solver.api.mam.aph2_fitting import aph2_fit_map
        from line_solver.api.mam.map_analysis import map_mean, map_normalize, map_scv
        aph1 = D1.copy()
        aph1[1, 1] = 0.0
        aph = list(map_normalize(D0, aph1))
        if map_scv(aph[0], aph[1]) < 1.0 + DEGENTOL:
            # hypoexponential: the canonical form stops being informative
            aph = aph2_fit_map([D0, D1])
        mean = map_mean(aph[0], aph[1])
        warnings.warn('%s: setting backward moments to ordinary moments, you should try '
                      'to fit B+S' % who, RuntimeWarning, stacklevel=2)
        mmap, _ = maph2m_fit_multiclass(aph, p, np.array([mean, mean]), cw)
        return mmap, _fwd(mmap), _sig(mmap), False
    elif abs(1 - r1) < DEGENTOL:
        # NON-CANONICAL PHASE-TYPE: two degrees of freedom, fit p and F, not S
        def degen_forward(vF1):
            return (p[0],                                          # meaningless here
                    p[0] * (h2 - vF1) / (h1 * (r2 - 1.0)),
                    p[0] * (h1 + h2 - vF1) / (h1 * r2))
        q1, q2, q3 = degen_forward(F[0])
        if not feasible(q1, q2, q3):
            q2_F = -p[0] / (h1 * (r2 - 1.0))
            q2_0 = p[0] * h2 / (h1 * (r2 - 1.0))
            q3_F = -p[0] / (h1 * r2)
            q3_0 = p[0] * (h1 + h2) / (h1 * r2)
            q1, q2, q3 = degen_forward(
                _project_interval([q2_F, q3_F], [q2_0, q3_0], F[0], who))
    elif form == 2 and r2 < DEGENTOL:
        # DEGENERATE CASE FOR gamma < 0: fit F or sigma, whichever weighs more
        if fw[0] > fw[1]:
            def degen_forward2(vF1):
                return (p[0] * (r1 - 2.0) * (h1 + h2 * r1 - vF1)
                        / ((r1 - 1.0) * (h1 - h2 + h2 * r1)),
                        -p[0] * (vF1 - h2) * (r1 - 2.0) / (h1 - h2 + h2 * r1),
                        p[0])                                      # meaningless here
            q1, q2, q3 = degen_forward2(F[0])
            if not feasible(q1, q2, q3):
                q1_F = -p[0] * (r1 - 2.0) / ((r1 - 1.0) * (h1 - h2 + h2 * r1))
                q1_0 = p[0] * (r1 - 2.0) * (h1 + h2 * r1) / ((r1 - 1.0) * (h1 - h2 + h2 * r1))
                q2_F = -p[0] * (r1 - 2.0) / (h1 - h2 + h2 * r1)
                q2_0 = p[0] * (r1 - 2.0) * h2 / (h1 - h2 + h2 * r1)
                q1, q2, q3 = degen_forward2(
                    _project_interval([q1_F, q2_F], [q1_0, q2_0], F[0], who))
        else:
            def degen_transition(vS11):
                # sqrt(p(1)^2 - S11) is COMPLEX above p(1)^2, where the reference's
                # own feasibility test then fails on the real part and sends it to
                # the repair. Report that as "not real" rather than raising, so the
                # repair below is reachable from the same inputs.
                root = p[0] ** 2 - vS11
                s = np.sqrt(max(root, 0.0))
                return (p[0] + s / (r1 - 1.0), p[0] + s, p[0]), root >= 0.0
            (q1, q2, q3), real = degen_transition(S[0, 0])
            if not real or not feasible(q1, q2, q3):
                # The reference states this repair as a YALMIP program, but it is
                # a one-variable convex QP -- minimize (x - S11)^2 over the
                # interval [max(0, p1^2(1-(1-r1)^2), p1^2-(1-p1)^2), p1^2] -- so
                # the projection below IS its global optimum.
                lo = max(0.0, p[0] ** 2 * (1.0 - (1.0 - r1) ** 2),
                         p[0] ** 2 - (1.0 - p[0]) ** 2)
                hi = p[0] ** 2
                if lo > hi:
                    raise ValueError('%s: the sigma repair of the gamma < 0 degeneracy is '
                                     'infeasible for this (p, r1)' % who)
                (q1, q2, q3), _ = degen_transition(min(max(float(S[0, 0]), lo), hi))
    else:
        # FULL FORM: the closed-form inverse
        if form == 1:
            G, U, Y = mamap2m_can1_coefficients(h1, h2, r1, r2)
            den = p[0] * (U[4] * F[0] + U[5])
            if abs(den) < DENUMTOL:
                q1 = q2 = q3 = p[0]
            else:
                q2 = (U[0] * F[0] ** 2 * p[0] ** 2 + U[1] * F[0] * p[0] ** 2
                      + U[2] * S[0, 0] + U[3] * p[0] ** 2) / den
                q1 = -(G[14] * p[0] - F[0] * G[2] * p[0]
                       + (G[2] * G[13] - G[1] * G[14]) * q2) / Y[1]
                q3 = +(G[12] * p[0] - F[0] * G[0] * p[0]
                       + (G[0] * G[13] - G[1] * G[12]) * q2) / Y[1]
        else:
            E, V, Z = mamap2m_can2_coefficients(h1, h2, r1, r2)
            den = V[4] * F[0] * p[0] + V[5] * p[0]
            if abs(den) < DENUMTOL:
                q1 = q2 = q3 = p[0]
            else:
                q3 = (V[0] * F[0] ** 2 * p[0] ** 2 + V[1] * F[0] * p[0] ** 2
                      + V[2] * p[0] ** 2 + V[3] * S[0, 0]) / den
                q1 = -(E[12] * p[0] - F[0] * E[1] * p[0]
                       + (E[1] * E[13] - E[2] * E[12]) * q3) / Z[1]
                q2 = +(E[11] * p[0] - F[0] * E[0] * p[0]
                       + (E[0] * E[13] - E[2] * E[11]) * q3) / Z[1]

        if feasible(q1, q2, q3):
            exact = True
        elif adjust:
            raise Mamap22Unsupported(
                '%s: the closed-form forward-plus-sigma inverse is infeasible for these '
                'targets, and the reference repairs it with a NONCONVEX bilinear program '
                'solved by YALMIP bmibnb, a spatial branch-and-bound returning a GLOBAL '
                'optimum. That solver is not ported; a local method would report a fit the '
                'reference would not have chosen. Pass adjust=False to take the clamped '
                'closed form, or use mamap2m_fit_fb_multiclass.' % who)

    if adjust and not feasible(q1, q2, q3) and not exact:
        raise ValueError('%s: Feasibility could not be restored: q1 = %e, q2 = %e, q3 = %e'
                         % (who, q1, q2, q3))

    mmap = [D0, D1] + _mark(D1, form, q1, q2, q3)
    return mmap, _fwd(mmap), _sig(mmap), exact


# ---------------------------------------------------------------------------
# backward moment + sigma
# ---------------------------------------------------------------------------

def mamap22_fit_bs_multiclass(map_: Sequence[np.ndarray], p: Sequence[float],
                              B: Sequence[float], S,
                              classWeights: Optional[Sequence[float]] = None,
                              bsWeights: Optional[Sequence[float]] = None,
                              adjust: bool = True
                              ) -> Tuple[List[np.ndarray], np.ndarray, np.ndarray, bool]:
    """Mark a canonical AMAP(2) on its BACKWARD moment and its sigma.

    The twin of :func:`mamap22_fit_fs_multiclass`; only the inverse differs,
    reading the backward block of the coefficient tables (U(7..12), G(10..12),
    Y(3)) where the forward fitter reads U(1..6), G(13..15) and Y(2).

    Args:
        map_: the underlying second-order AMAP as [D0, D1]
        p: the two class probabilities
        B: the two target backward moments
        S: the target class transition matrix; only S[0, 0] enters
        classWeights: per-class weights in the objective (default: unit)
        bsWeights: (backward, sigma) weights (default: unit)
        adjust: repair an infeasible closed form; the repair is unported

    Returns:
        (mmap, fB, fS, exact)
    """
    who = 'Fitting MAMAP(2,2) B+S'
    D0 = np.asarray(map_[0], dtype=float)
    D1 = np.asarray(map_[1], dtype=float)
    form = _canonical_form(D0, D1, who)

    p = np.asarray(p, dtype=float).ravel()
    B = np.asarray(B, dtype=float).ravel()
    S = np.atleast_2d(np.asarray(S, dtype=float))
    if p.size != 2:
        raise ValueError('%s: fitting backward and transition probabilities only supports '
                         'two classes' % who)
    if B.size != 2:
        raise ValueError('%s: one backward moment per class is required' % who)
    cw = np.ones(2) if classWeights is None else np.asarray(classWeights, float).ravel()
    bw = np.ones(2) if bsWeights is None else np.asarray(bsWeights, float).ravel()

    feastol = 1e-4
    h1 = -1.0 / D0[0, 0]
    h2 = -1.0 / D0[1, 1]
    r1 = D0[0, 1] * h1
    r2 = D1[1, 1] * h2

    def feasible(*qs) -> bool:
        return all(-feastol <= q <= 1.0 + feastol for q in qs)

    # ---- the Poisson perturbation ----------------------------------------
    if ((form == 1 and (r1 < DEGENTOL or r2 > 1 - DEGENTOL
                        or abs(h2 - h1 * r2) < DEGENTOL))
            or (form == 2 and (r2 > 1 - DEGENTOL
                               or abs(h1 - h2 - h1 * r1 + h1 * r1 * r2) < DEGENTOL))):
        if r1 < DEGENTOL:
            r1 = DEGENTOL
        if r2 > 1 - DEGENTOL:
            r2 = 1 - DEGENTOL
        if form == 1 and abs(h2 - h1 * r2) < DEGENTOL:
            h2 = h1 * r2 + DEGENTOL
        if form == 2 and abs(h1 - h2 - h1 * r1 + h1 * r1 * r2) < DEGENTOL:
            h1 = (h2 + DEGENTOL) / (1.0 - r1 + r1 * r2)
        D0, D1 = _rebuild(h1, h2, r1, r2, form)

    exact = False
    q1 = q2 = q3 = float('nan')

    if form == 2 and r2 < DEGENTOL and abs(1 - r1) < DEGENTOL:
        # DEGENERATE PHASE-TYPE: only the class probabilities are identifiable
        q1 = q2 = q3 = p[0]
    elif form == 1 and r2 < DEGENTOL:
        # CANONICAL PHASE-TYPE: the problem IS the MAPH fit
        from line_solver.api.mam.map_analysis import map_normalize
        aph1 = D1.copy()
        aph1[1, 1] = 0.0
        aph = list(map_normalize(D0, aph1))
        mmap, _ = maph2m_fit_multiclass(aph, p, B, cw)
        return mmap, _bwd(mmap), _sig(mmap), False
    elif abs(1 - r1) < DEGENTOL:
        # NON-CANONICAL PHASE-TYPE: refit the timing as an APH(2), then mark it
        from line_solver.api.mam.aph2_fitting import aph2_fit_map
        aph = aph2_fit_map([D0, D1])
        mmap, _ = maph2m_fit_multiclass(aph, p, B, cw)
        return mmap, _bwd(mmap), _sig(mmap), False
    elif form == 2 and r2 < DEGENTOL:
        # DEGENERATE CASE FOR gamma < 0: fit B or sigma, whichever weighs more
        if bw[0] > bw[1]:
            def degen_backward(vB1):
                return ((p[0] * (r1 - 2.0) * (h2 - vB1 + h1 * r1))
                        / ((r1 - 1.0) * (h2 - h1 + h1 * r1)),
                        -(p[0] * (vB1 - h1) * (r1 - 2.0)) / (h2 - h1 + h1 * r1),
                        p[0])                                      # meaningless here
            q1, q2, q3 = degen_backward(B[0])
            if not feasible(q1, q2, q3):
                q1_B = -p[0] * (r1 - 2.0) / ((r1 - 1.0) * (h2 - h1 + h1 * r1))
                q1_0 = p[0] * (r1 - 2.0) * (h2 + h1 * r1) / ((r1 - 1.0) * (h2 - h1 + h1 * r1))
                q2_B = -p[0] * (r1 - 2.0) / (h2 - h1 + h1 * r1)
                q2_0 = p[0] * (r1 - 2.0) * h1 / (h2 - h1 + h1 * r1)
                q1, q2, q3 = degen_backward(
                    _project_interval([q1_B, q2_B], [q1_0, q2_0], B[0], who))
        else:
            def degen_transition(vS11):
                # sqrt(p(1)^2 - S11) is COMPLEX above p(1)^2, where the reference's
                # own feasibility test then fails on the real part and sends it to
                # the repair. Report that as "not real" rather than raising, so the
                # repair below is reachable from the same inputs.
                root = p[0] ** 2 - vS11
                s = np.sqrt(max(root, 0.0))
                return (p[0] + s / (r1 - 1.0), p[0] + s, p[0]), root >= 0.0
            (q1, q2, q3), real = degen_transition(S[0, 0])
            if not real or not feasible(q1, q2, q3):
                safetytol = 1e-10
                q1lb = p[0] ** 2 * (1.0 - (1.0 - r1) ** 2)
                q2ub = p[0] ** 2 - (1.0 - p[0]) ** 2
                s11 = float(S[0, 0])
                if s11 <= q1lb:
                    s11 = q1lb + safetytol
                elif s11 >= q2ub:
                    s11 = q2ub - safetytol
                (q1, q2, q3), _ = degen_transition(s11)
    else:
        # FULL FORM: the closed-form inverse
        if form == 1:
            G, U, Y = mamap2m_can1_coefficients(h1, h2, r1, r2)
            den = U[10] * B[0] * p[0] + U[11] * p[0]
            if abs(den) < DENUMTOL:
                q1 = q2 = q3 = p[0]
            else:
                q2 = (U[6] * B[0] ** 2 * p[0] ** 2 + U[7] * B[0] * p[0] ** 2
                      + U[8] * S[0, 0] + U[9] * p[0] ** 2) / den
                q1 = -(G[11] * p[0] - B[0] * G[2] * p[0]
                       + (G[2] * G[10] - G[1] * G[11]) * q2) / Y[2]
                q3 = +(G[9] * p[0] - B[0] * G[0] * p[0]
                       + (G[0] * G[10] - G[1] * G[9]) * q2) / Y[2]
        else:
            E, V, Z = mamap2m_can2_coefficients(h1, h2, r1, r2)
            den = V[10] * B[0] * p[0] + V[11] * p[0]
            if abs(den) < DENUMTOL:
                q1 = q2 = q3 = p[0]
            else:
                q3 = (V[6] * B[0] ** 2 * p[0] ** 2 + V[7] * B[0] * p[0] ** 2
                      + V[8] * p[0] ** 2 + V[9] * S[0, 0]) / den
                q1 = +(E[9] * p[0] - B[0] * E[1] * p[0]
                       + (E[1] * E[10] - E[2] * E[9]) * q3) / Z[2]
                q2 = -(E[8] * p[0] - B[0] * E[0] * p[0]
                       + (E[0] * E[10] - E[2] * E[8]) * q3) / Z[2]

        if feasible(q1, q2, q3):
            exact = True
        elif adjust:
            raise Mamap22Unsupported(
                '%s: the closed-form backward-plus-sigma inverse is infeasible for these '
                'targets, and the reference repairs it with a NONCONVEX bilinear program '
                'solved by YALMIP bmibnb, a spatial branch-and-bound returning a GLOBAL '
                'optimum. That solver is not ported; a local method would report a fit the '
                'reference would not have chosen. Pass adjust=False to take the clamped '
                'closed form, weight the forward moment above sigma so '
                'mamap2m_fit_fb_multiclass applies, or relax the targets.' % who)

    if adjust and not feasible(q1, q2, q3) and not exact:
        raise ValueError('%s: Feasibility could not be restored' % who)

    mmap = [D0, D1] + _mark(D1, form, q1, q2, q3)
    return mmap, _bwd(mmap), _sig(mmap), exact


# ---------------------------------------------------------------------------
# drivers over every AMAP(2) form
# ---------------------------------------------------------------------------

def _amap2_forms(M1: float, M2: float, M3: float, GAMMA: float, P: np.ndarray):
    """The candidate AMAP(2) forms, or a marked Poisson when none is second order.

    Returns (maps, poisson): exactly one of the two is not None. When the moment
    set admits only a one-state process the reference perturbs M2 and M3 just
    above the exponential to recover a two-state form, and falls back to a
    marked Poisson only if that also fails.
    """
    _, maps = amap2_fit_gamma(M1, M2, M3, GAMMA)
    if len(maps) == 1 and np.asarray(maps[0][0]).shape[0] == 1:
        M2a = M2 * (1 + 1e-4)
        M3a = M3 * (M2a / M2) ** 1.5
        maps2 = amap2_fitall_gamma(M1, M2a, M3a, GAMMA)
        if maps2:
            return [_map_repair(x[0], x[1]) for x in maps2], None
        D0 = np.asarray(maps[0][0], float)
        D1 = np.asarray(maps[0][1], float)
        return None, [D0.copy(), D1.copy()] + [D1 * P[c] for c in range(P.size)]
    return maps, None


def mamap22_fit_gamma_fs(M1: float, M2: float, M3: float, GAMMA: float,
                         P: Sequence[float], F: Sequence[float], S
                         ) -> List[np.ndarray]:
    """MAMAP(2,2) matching the moments and gamma, p exactly, and F and sigma
    approximately. Fits over every AMAP(2) form and keeps the closest."""
    P = np.asarray(P, float).ravel()
    F = np.asarray(F, float).ravel()
    S = np.atleast_2d(np.asarray(S, float))

    maps, poisson = _amap2_forms(M1, M2, M3, GAMMA, P)
    if poisson is not None:
        return poisson

    best, best_err = None, np.inf
    for m in maps:
        try:
            mmap, fF, fS, _ = mamap22_fit_fs_multiclass(m, P, F, S)
        except (Mamap22Unsupported, ValueError):
            continue
        # the reference scores on the FIRST class alone: fF(1) and fS(1,1)
        err = (fF[0] / F[0] - 1) ** 2 + (fS[0, 0] / S[0, 0] - 1) ** 2
        if err < best_err:
            best, best_err = mmap, err
    if best is None:
        raise ValueError('mamap22_fit_gamma_fs: no AMAP(2) form admits a feasible '
                         'forward-plus-sigma marking for these targets')
    return best


def mamap22_fit_gamma_bs(M1: float, M2: float, M3: float, GAMMA: float,
                         P: Sequence[float], B: Sequence[float], S
                         ) -> List[np.ndarray]:
    """MAMAP(2,2) matching the moments and gamma, p exactly, and B and sigma
    approximately. Fits over every AMAP(2) form and keeps the closest."""
    P = np.asarray(P, float).ravel()
    B = np.asarray(B, float).ravel()
    S = np.atleast_2d(np.asarray(S, float))

    maps, poisson = _amap2_forms(M1, M2, M3, GAMMA, P)
    if poisson is not None:
        return poisson

    best, best_err = None, np.inf
    for m in maps:
        try:
            mmap, fB, fS, _ = mamap22_fit_bs_multiclass(m, P, B, S)
        except (Mamap22Unsupported, ValueError):
            continue
        err = (fB[0] / B[0] - 1) ** 2 + (fS[0, 0] / S[0, 0] - 1) ** 2
        if err < best_err:
            best, best_err = mmap, err
    if best is None:
        raise ValueError('mamap22_fit_gamma_bs: no AMAP(2) form admits a feasible '
                         'backward-plus-sigma marking for these targets')
    return best


def _mmap_characteristics(mmap: Sequence[np.ndarray]):
    from line_solver.api.mam.map_analysis import map_gamma, map_moment
    from line_solver.api.mam.mmap_ops import mmap_pc
    D0 = np.asarray(mmap[0], float)
    D1 = np.asarray(mmap[1], float)
    return (map_moment(D0, D1, 1), map_moment(D0, D1, 2), map_moment(D0, D1, 3),
            map_gamma(D0, D1), np.asarray(mmap_pc(list(mmap)), float).ravel())


def mamap22_fit_gamma_fs_mmap(mmap: Sequence[np.ndarray]) -> List[np.ndarray]:
    """:func:`mamap22_fit_gamma_fs` driven from an MMAP[2] of arbitrary order."""
    M1, M2, M3, GAMMA, P = _mmap_characteristics(mmap)
    return mamap22_fit_gamma_fs(M1, M2, M3, GAMMA, P, _fwd(mmap), _sig(mmap))


def mamap22_fit_gamma_bs_mmap(mmap: Sequence[np.ndarray]) -> List[np.ndarray]:
    """:func:`mamap22_fit_gamma_bs` driven from an MMAP[2] of arbitrary order."""
    M1, M2, M3, GAMMA, P = _mmap_characteristics(mmap)
    return mamap22_fit_gamma_bs(M1, M2, M3, GAMMA, P, _bwd(mmap), _sig(mmap))


def _trace_characteristics(T: Sequence[float], A: Sequence[int]):
    from line_solver.api.trace.trace_analysis import mtrace_pc, trace_gamma
    T = np.asarray(T, dtype=float).ravel()
    A = np.asarray(A).ravel()
    if T.size == 0 or T.size != A.size:
        raise ValueError('the trace and its labels must agree in length')
    gamma = np.asarray(trace_gamma(T), dtype=float).ravel()
    return (float(np.mean(T)), float(np.mean(T ** 2)), float(np.mean(T ** 3)),
            float(gamma[0]) if gamma.size else 0.0,
            np.asarray(mtrace_pc(T, A), dtype=float).ravel(), T, A)


def mamap22_fit_gamma_fs_trace(T: Sequence[float], C: Sequence[int]) -> List[np.ndarray]:
    """:func:`mamap22_fit_gamma_fs` driven from a marked trace."""
    from line_solver.api.trace.trace_analysis import mtrace_forward_moment, mtrace_sigma
    M1, M2, M3, GAMMA, P, T, C = _trace_characteristics(T, C)
    F = np.asarray(mtrace_forward_moment(T, C, [1]), dtype=float).ravel()
    S = np.atleast_2d(np.asarray(mtrace_sigma(T, C), dtype=float))
    return mamap22_fit_gamma_fs(M1, M2, M3, GAMMA, P, F, S)


def mamap22_fit_gamma_bs_trace(T: Sequence[float], A: Sequence[int]) -> List[np.ndarray]:
    """:func:`mamap22_fit_gamma_bs` driven from a marked trace."""
    from line_solver.api.trace.trace_analysis import mtrace_backward_moment, mtrace_sigma
    M1, M2, M3, GAMMA, P, T, A = _trace_characteristics(T, A)
    B = np.asarray(mtrace_backward_moment(T, A, [1]), dtype=float).ravel()
    S = np.atleast_2d(np.asarray(mtrace_sigma(T, A), dtype=float))
    return mamap22_fit_gamma_bs(M1, M2, M3, GAMMA, P, B, S)
