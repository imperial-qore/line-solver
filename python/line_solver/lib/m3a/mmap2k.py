"""MMAP(2,K): closed-form fit of a marked MAP of second order with K classes.

The fit factorizes into two independent inverse problems:

1. the underlying MAP(2), obtained from (M1, M2, M3, GAMMA) by the acyclic
   canonical inverse of ``amap2_fit_gamma``. Order two loses nothing here:
   every MAP(2) is equivalent to one of the two acyclic canonical forms;
2. the marking, which is LINEAR in the class characteristics. Writing the
   marking as fractions of the aggregate D1,

       form 1:  D1^(c) = D1 .* [[q1c, 0  ], [q2c, q3c]]
       form 2:  D1^(c) = D1 .* [[0,   q1c], [q2c, q3c]]

   the map (q1c, q2c, q3c) -> (p_c, p_c F_c, p_c B_c) is linear with a 3x3
   matrix that depends only on (h1, h2, r1, r2) and is the SAME for every
   class, so the inverse is one pre-computed 3x3 block applied per class and
   costs nothing as K grows. The block was derived in SageMath, see
   sage/proofs/mmap2k_marking_inverse.py; the offset of the affine map is zero,
   so the formulas below are purely linear in (1, F_c, B_c) scaled by p_c.

Degrees of freedom: the canonical form has 4 + 3(K-1) parameters and realizes
3K+1 independent characteristics, the same count as a general MMAP(2,K)
(verified by Jacobian rank for K up to 6), so nothing is given up by working in
the canonical form.

The inverse is singular exactly on the six degenerate loci r1 in {0,1},
r2 in {0,1}, h1 - h2 + h2 r1 = 0 and (form 1) h1 r2 - h2 = 0 or (form 2)
h1 r1 r2 - h1 r1 + h1 - h2 = 0, which are the branches
``mamap2m_fit_fb_multiclass`` special-cases. There, and whenever the closed
form leaves the unit box, this module falls back to that quadratic program.
"""
from typing import List, Optional, Sequence, Tuple

import numpy as np

DEGENTOL = 1e-8
FEASTOL = 1e-8


def _canonical_parameters(D0: np.ndarray, D1: np.ndarray) -> Tuple[int, float, float, float, float]:
    """(form, h1, h2, r1, r2) of a canonical acyclic MAP(2)."""
    if D0.shape[0] != 2:
        raise ValueError('mmap2k_fit: the underlying MAP must be of second order')
    if D0[1, 0] != 0:
        raise ValueError('mmap2k_fit: the underlying MAP must be acyclic')
    if D1[0, 1] == 0:
        form = 1
    elif D1[0, 0] == 0:
        form = 2
    else:
        raise ValueError('mmap2k_fit: the underlying MAP must be in canonical acyclic form')
    h1 = -1.0 / D0[0, 0]
    h2 = -1.0 / D0[1, 1]
    r1 = D0[0, 1] * h1
    r2 = D1[1, 1] * h2
    return form, h1, h2, r1, r2


def marking_inverse(form: int, h1: float, h2: float, r1: float, r2: float,
                    p: float, Fc: float, Bc: float) -> Tuple[float, float, float]:
    """Closed-form marking fractions of one class.

    Pre-computed inverse of the linear map (q1, q2, q3) -> (p, p F, p B); see
    the module docstring. Raises ZeroDivisionError-free ValueError on the
    singular loci so the caller can fall back.
    """
    if form == 1:
        d1 = (r2 - 1) * (r1 - 1) * (h1 * r2 - h2)
        d2 = (r2 - 1) * r1
        d3 = r1 * r2 * (h1 + h2 * r1 - h2)
        if min(abs(d1), abs(d2), abs(d3), abs(h1 + h2 * r1 - h2), abs(h1 * r2 - h2)) < DEGENTOL:
            raise ValueError('mmap2k_fit: degenerate underlying form')
        W = r1 * r2 - r2 + 1
        q1 = p * W * ((h1 * r2 - h1 - h2) + Bc) / d1
        q2 = p * W * ((h1 ** 2 * (r2 - 1) + h1 * h2 * r1 * (r2 - 1) - h2 ** 2 * r1)
                      / ((h1 + h2 * r1 - h2) * (h1 * r2 - h2))
                      - Fc / (h1 + h2 * r1 - h2)
                      + Bc / (h1 * r2 - h2)) / d2
        q3 = p * W * ((h1 + h2 * r1) - Fc) / d3
    else:
        U = h1 * r1 * r2 - h1 * r1 + h1 - h2
        d1 = (r2 - 1) * (r1 - 1) * U
        d2 = (r2 - 1) * (h1 + h2 * r1 - h2)
        if min(abs(d1), abs(d2), abs(r2), abs(U), abs(h1 + h2 * r1 - h2)) < DEGENTOL:
            raise ValueError('mmap2k_fit: degenerate underlying form')
        V = r1 * r2 - r1 - r2 + 2
        q1 = p * V * ((h1 * r1 * r2 - h1 * r1 - h2) + Bc) / d1
        q2 = p * V * (h2 - Fc) / d2
        q3 = p * V * ((h1 ** 2 + h1 * h2 * r1 * r2 - h2 ** 2) / ((h1 + h2 * r1 - h2) * U)
                      - Fc / (h1 + h2 * r1 - h2)
                      - Bc / U) / r2
    return q1, q2, q3


def _assemble(D0: np.ndarray, D1: np.ndarray, form: int,
              q: np.ndarray) -> List[np.ndarray]:
    out = [np.array(D0, dtype=float), np.array(D1, dtype=float)]
    K = q.shape[1]
    for c in range(K):
        if form == 1:
            mask = np.array([[q[0, c], 0.0], [q[1, c], q[2, c]]])
        else:
            mask = np.array([[0.0, q[0, c]], [q[1, c], q[2, c]]])
        out.append(D1 * mask)
    return out


def mmap2k_fit(M1: float, M2: float, M3: float, GAMMA: float,
               P: Sequence[float], F: Sequence[float], B: Sequence[float],
               exact_only: bool = False) -> List[np.ndarray]:
    """Fit an MMAP(2,K) to the inter-arrival moments, the decay rate and the
    per-class probabilities, forward moments and backward moments.

    The three inter-arrival moments, the decay rate and the class
    probabilities are matched exactly; the forward and backward moments are
    matched exactly whenever the closed-form marking lands in the unit box,
    which is the generic case. Otherwise the fit falls back to the quadratic
    program of ``mamap2m_fit_fb_multiclass`` (unless exact_only is set).

    Args:
        M1, M2, M3: raw moments of the inter-arrival times
        GAMMA: autocorrelation decay rate
        P: class probabilities, summing to one
        F: first-order forward moments, with sum_c P_c F_c = M1
        B: first-order backward moments, with sum_c P_c B_c = M1
        exact_only: raise instead of falling back to the quadratic program

    Returns:
        The MMAP as [D0, D1, D1^(1), ..., D1^(K)]
    """
    from .amap2 import amap2_fit_gamma

    P = np.asarray(P, dtype=float).ravel()
    F = np.asarray(F, dtype=float).ravel()
    B = np.asarray(B, dtype=float).ravel()
    K = P.size
    if F.size != K or B.size != K:
        raise ValueError('mmap2k_fit: P, F and B must have the same length')

    _, candidates = amap2_fit_gamma(M1, M2, M3, GAMMA)

    best = None
    best_err = np.inf
    for cand in candidates:
        D0 = np.asarray(cand[0], dtype=float)
        D1 = np.asarray(cand[1], dtype=float)
        if D0.shape[0] != 2:
            continue
        try:
            form, h1, h2, r1, r2 = _canonical_parameters(D0, D1)
            q = np.empty((3, K))
            for c in range(K):
                q[:, c] = marking_inverse(form, h1, h2, r1, r2, P[c], F[c], B[c])
        except ValueError:
            continue
        # feasibility: fractions in the unit box and summing to one per phase
        box = float(min(q.min(), 1.0 - q.max()))
        closure = float(np.max(np.abs(q.sum(axis=1) - 1.0)))
        err = max(-box, 0.0) + closure
        if err < best_err:
            best_err = err
            best = (D0, D1, form, q)

    if best is not None and best_err <= FEASTOL:
        D0, D1, form, q = best
        return _assemble(D0, D1, form, np.clip(q, 0.0, 1.0))

    if exact_only:
        raise ValueError('mmap2k_fit: no exact closed-form marking for these '
                         'characteristics (worst violation %.3g)' % best_err)

    from .mamap2m import mamap2m_fit_gamma_fb
    return mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B)
