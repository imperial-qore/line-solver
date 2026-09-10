"""Drift guard for the native-Python m3a multiclass fitters.

``amap2_*``, ``maph2m_*`` and ``mamap2m_*`` are ports of the MATLAB m3a files of
the same names. Assertions are of two kinds:

- INVARIANTS: the fit reproduces the moments, decay rate and class
  probabilities it is meant to match, and is a valid MMAP.
- PINNED VALUES: matrices taken from MATLAB for the same input, so a change
  that stays self-consistent but is no longer the same algorithm is caught.

See _kb/03-api-layer.md, "Closed-form fitter audit".
"""
import math

import numpy as np
import pytest

from line_solver.lib.m3a.amap2 import amap2_fit_gamma, amap2_fitall_gamma
from line_solver.lib.m3a.maph2m import maph2m_fit, maph2m_fit_multiclass
from line_solver.lib.m3a.mamap2m import (mamap2m_fit, mamap2m_fit_fb_multiclass,
                                         mamap2m_fit_gamma_fb, mamap2m_fit_trace)
from line_solver.lib.m3a.utils import validate_mmap
from line_solver.api.mam.mmap_ops import mmap_pc


def _stats(D0, D1, k=3):
    D0 = np.asarray(D0, float)
    D1 = np.asarray(D1, float)
    n = D0.shape[0]
    A = np.linalg.inv(-D0)
    P = A @ D1
    one = np.ones(n)
    w, V = np.linalg.eig(P.T)
    pi = np.real(V[:, np.argmin(abs(w - 1))])
    pi = pi / pi.sum()
    mom = [float(math.factorial(i) * (pi @ np.linalg.matrix_power(A, i) @ one))
           for i in range(1, k + 1)]
    ev = np.sort(np.real(np.linalg.eigvals(P)))[::-1]
    g2 = float(ev[1]) if n > 1 else 0.0
    return mom, g2


CASES = [
    (1, 3, 15, 0.3, [0.6, 0.4], [0.55, 0.45], [1.0769230769230769, 0.88461538461538458]),
    (1, 5, 45, 0.5, [0.7, 0.3], [0.6, 0.4], [1.05, 0.88333333333333333]),
    (1, 2.5, 10, -0.2, [0.6, 0.4], [0.52, 0.48], [1.02, 0.97]),
]
S22 = np.array([[0.4, 0.2], [0.2, 0.2]])


@pytest.mark.parametrize("M1,M2,M3,GAMMA", [(1, 3, 15, 0.3), (1, 5, 45, 0.5), (1, 2.5, 10, -0.2)])
def test_amap2_fit_gamma_matches_moments_and_decay(M1, M2, M3, GAMMA):
    amap, amaps = amap2_fit_gamma(M1, M2, M3, GAMMA)
    assert len(amaps) >= 1
    mom, g2 = _stats(amap[0], amap[1])
    assert mom[0] == pytest.approx(M1, rel=1e-9)
    assert mom[1] == pytest.approx(M2, rel=1e-9)
    assert mom[2] == pytest.approx(M3, rel=1e-9)
    assert g2 == pytest.approx(GAMMA, rel=1e-8, abs=1e-12)
    # every candidate form matches the same characteristics
    for cand in amaps:
        m, g = _stats(cand[0], cand[1])
        assert m[0] == pytest.approx(M1, rel=1e-8)
        assert g == pytest.approx(GAMMA, rel=1e-7, abs=1e-12)


def test_amap2_fit_gamma_matches_matlab():
    # Pinned against MATLAB amap2_fit_gamma.m
    amap, _ = amap2_fit_gamma(1, 3, 15, 0.3)
    D0 = np.asarray(amap[0], float)
    D1 = np.asarray(amap[1], float)
    assert D0[0, 0] == pytest.approx(-3.4142135623730963, abs=1e-12)
    assert D0[0, 1] == pytest.approx(1.0701829114820407, abs=1e-12)
    assert D1[0, 0] == pytest.approx(2.3440306508910553, abs=1e-12)
    assert D1[1, 1] == pytest.approx(0.25596934910894475, abs=1e-12)


@pytest.mark.parametrize("M1,M2,M3,GAMMA,P,F,B", CASES)
def test_maph2m_fit_matches_class_probabilities(M1, M2, M3, GAMMA, P, F, B):
    maph = maph2m_fit(M1, M2, M3, P, B)
    assert validate_mmap(maph)
    pc = np.asarray(mmap_pc(maph), float).ravel()
    assert pc == pytest.approx(np.asarray(P), rel=1e-8)
    mom, _ = _stats(maph[0], maph[1])
    assert mom[0] == pytest.approx(M1, rel=1e-8)
    assert mom[1] == pytest.approx(M2, rel=1e-8)


def test_maph2m_fit_matches_matlab():
    # Pinned against MATLAB maph2m_fit.m
    maph = maph2m_fit(1, 3, 15, [0.6, 0.4], [1.0769230769230769, 0.88461538461538458])
    assert np.asarray(maph[2])[0, 0] == pytest.approx(1.107692307692308, abs=1e-9)
    assert np.asarray(maph[2])[1, 0] == pytest.approx(0.38970696064135169, abs=1e-9)
    assert np.asarray(maph[3])[0, 0] == pytest.approx(0.89230769230769247, abs=1e-9)
    assert np.asarray(maph[3])[1, 0] == pytest.approx(0.1960794769855532, abs=1e-9)


def test_maph2m_fit_degenerate_form_returns_backward_moments():
    # (2, 9, 60) drives the degenerate branch, where the MATLAB reference used
    # to leave fB unassigned; the fit must still come back complete.
    maph = maph2m_fit(2, 9, 60, [0.5, 0.5], [2.1, 1.9])
    assert validate_mmap(maph)
    assert np.asarray(maph[2])[1, 0] == pytest.approx(0.5, abs=1e-9)
    assert np.asarray(maph[3])[1, 0] == pytest.approx(0.5, abs=1e-9)


@pytest.mark.parametrize("M1,M2,M3,GAMMA,P,F,B", CASES)
def test_mamap2m_fit_matches_moments_decay_and_class_probabilities(M1, M2, M3, GAMMA, P, F, B):
    mmap = mamap2m_fit(M1, M2, M3, GAMMA, P, F, B, S22)
    assert validate_mmap(mmap)
    mom, g2 = _stats(mmap[0], mmap[1])
    assert mom[0] == pytest.approx(M1, rel=1e-8)
    assert mom[1] == pytest.approx(M2, rel=1e-8)
    assert mom[2] == pytest.approx(M3, rel=1e-8)
    assert g2 == pytest.approx(GAMMA, rel=1e-7, abs=1e-12)
    pc = np.asarray(mmap_pc(mmap), float).ravel()
    assert pc == pytest.approx(np.asarray(P), rel=1e-6)


def test_mamap2m_fit_matches_matlab():
    # Pinned against MATLAB mamap2m_fit.m
    mmap = mamap2m_fit(1, 3, 15, 0.3, [0.6, 0.4], [0.55, 0.45],
                       [1.0769230769230769, 0.88461538461538458], S22)
    assert np.asarray(mmap[2])[0, 0] == pytest.approx(1.2823540227600763, abs=1e-9)
    assert np.asarray(mmap[2])[1, 0] == pytest.approx(0.1337376115324069, abs=1e-9)
    assert np.asarray(mmap[2])[1, 1] == pytest.approx(0.2559693491089442, abs=1e-9)
    assert np.asarray(mmap[3])[0, 0] == pytest.approx(1.061676628130982, abs=1e-9)
    assert np.asarray(mmap[3])[1, 0] == pytest.approx(0.19607947698555314, abs=1e-9)


def test_mamap2m_fit_gamma_fb_matches_matlab():
    mmap = mamap2m_fit_gamma_fb(1, 3, 15, 0.3, [0.6, 0.4], [0.55, 0.45],
                                [1.0769230769230769, 0.88461538461538458])
    assert np.asarray(mmap[2])[1, 1] == pytest.approx(1.1170665688778327, abs=1e-9)
    assert np.asarray(mmap[3])[1, 1] == pytest.approx(1.2269640820132226, abs=1e-9)


def test_mamap2m_fit_fb_multiclass_reports_realized_moments():
    amap, _ = amap2_fit_gamma(1, 3, 15, 0.3)
    mmap, fF, fB = mamap2m_fit_fb_multiclass(
        amap, [0.6, 0.4], [0.55, 0.45], [1.0769230769230769, 0.88461538461538458])
    assert validate_mmap(mmap)
    assert fF.size == 2 and fB.size == 2
    # the realized backward moments must satisfy sum_c p_c b_c = M1
    assert float(np.asarray([0.6, 0.4]) @ fB) == pytest.approx(1.0, rel=1e-6)


def test_trace_fitting_matches_the_trace_statistics():
    # the fabricated placeholder this replaced ignored the variance entirely
    rng = np.random.default_rng(7)
    n = 20000
    u = rng.random(n)
    S = np.where(u < 0.2, rng.exponential(2.6, n), rng.exponential(0.6, n))
    S = S / S.mean()
    C = (rng.random(n) > 0.6).astype(int)
    mmap = mamap2m_fit_trace(S, C)
    assert validate_mmap(mmap)
    mom, _ = _stats(mmap[0], mmap[1], 2)
    scv = (mom[1] - mom[0] ** 2) / mom[0] ** 2
    assert mom[0] == pytest.approx(float(S.mean()), rel=1e-6)
    assert scv == pytest.approx(float(S.var() / S.mean() ** 2), rel=1e-6)
    pc = np.asarray(mmap_pc(mmap), float).ravel()
    assert pc[0] == pytest.approx(float(np.mean(C == 0)), rel=1e-6)


@pytest.mark.parametrize("K", [2, 3, 4, 6, 10])
@pytest.mark.parametrize("form", [1, 2])
def test_mmap2k_fit_recovers_a_canonical_mmap_exactly(K, form):
    """The closed-form MMAP(2,K) fit must return the process it was given.

    A random canonical MMAP(2,K) is built, its characteristics measured, and the
    fit must reproduce all of them: three moments, decay rate, class
    probabilities, forward and backward moments. The marking inverse is linear
    and K-independent (see io/sage/proofs/mmap2k_marking_inverse.py), so this is an
    exact-recovery test, not a tolerance one.
    """
    from line_solver.lib.m3a.amap2 import amap2_assemble
    from line_solver.lib.m3a.mmap2k import mmap2k_fit
    from line_solver.api.mam.map_analysis import map_gamma2, map_moment
    from line_solver.api.mam.mmap_ops import (mmap_backward_moment,
                                              mmap_forward_moment, mmap_pc)

    rng = np.random.default_rng(1000 + 10 * K + form)
    h1 = 0.4 + 0.6 * rng.random()
    h2 = 1.2 + 0.8 * rng.random()
    r1 = 0.25 + 0.5 * rng.random()
    r2 = 0.2 + 0.5 * rng.random()
    D0, D1 = (np.asarray(x, float) for x in amap2_assemble(h1, h2, r1, r2, form))
    q = rng.random((3, K))
    q /= q.sum(axis=1, keepdims=True)
    src = [D0, D1]
    for c in range(K):
        mask = (np.array([[q[0, c], 0.0], [q[1, c], q[2, c]]]) if form == 1
                else np.array([[0.0, q[0, c]], [q[1, c], q[2, c]]]))
        src.append(D1 * mask)

    m1 = map_moment(D0, D1, 1)
    m2 = map_moment(D0, D1, 2)
    m3 = map_moment(D0, D1, 3)
    g = map_gamma2(D0, D1)
    P = np.ravel(mmap_pc(src))
    Fm = np.ravel(mmap_forward_moment(src, [1]))
    Bm = np.ravel(mmap_backward_moment(src, [1]))

    fit = mmap2k_fit(m1, m2, m3, g, P, Fm, Bm, exact_only=True)
    fD0 = np.asarray(fit[0], float)
    fD1 = np.asarray(fit[1], float)
    assert map_moment(fD0, fD1, 1) == pytest.approx(m1, rel=1e-9)
    assert map_moment(fD0, fD1, 2) == pytest.approx(m2, rel=1e-9)
    assert map_moment(fD0, fD1, 3) == pytest.approx(m3, rel=1e-9)
    assert map_gamma2(fD0, fD1) == pytest.approx(g, rel=1e-8)
    assert np.ravel(mmap_pc(fit)) == pytest.approx(P, rel=1e-8)
    assert np.ravel(mmap_forward_moment(fit, [1])) == pytest.approx(Fm, rel=1e-8)
    assert np.ravel(mmap_backward_moment(fit, [1])) == pytest.approx(Bm, rel=1e-8)
    assert validate_mmap(fit)


@pytest.mark.parametrize("n", [2, 3])
@pytest.mark.parametrize("K", [2, 3, 4, 6])
def test_mmap3k_fit_recovers_the_marking_exactly(n, K):
    """Order-3 marking fit: exact recovery of a random marked MAP(n).

    Order two needs (p, F, B); order three needs the second-order backward
    moment as well, because the canonical D1 has one more nonzero. The
    independent characteristic set was derived in
    io/sage/proofs/mmap3k_marking_inverse.py.
    """
    from line_solver.lib.m3a.mmap3k import mmap3k_fit

    rng = np.random.default_rng(2000 + 10 * n + K)
    h = 0.4 + rng.random(n)
    r = 0.2 + 0.5 * rng.random(max(n - 1, 1))
    s = 0.2 + 0.5 * rng.random()
    D0 = np.zeros((n, n))
    D1 = np.zeros((n, n))
    for i in range(n):
        D0[i, i] = -1.0 / h[i]
        if i + 1 < n:
            D0[i, i + 1] = r[i] / h[i]
    for i in range(n):
        D1[i, 0] = (1 - (r[i] if i + 1 < n else s)) / h[i]
    D1[n - 1, n - 1] = s / h[n - 1]

    nz = [(i, j) for i in range(n) for j in range(n) if D1[i, j] != 0.0]
    z = len(nz)
    q = rng.random((z, K))
    q /= q.sum(axis=1, keepdims=True)
    src = [D0, D1]
    for c in range(K):
        Dc = np.zeros_like(D1)
        for jj, (i, j) in enumerate(nz):
            Dc[i, j] = D1[i, j] * q[jj, c]
        src.append(Dc)

    def chars(mmap):
        A = np.linalg.inv(-np.asarray(mmap[0], float))
        P = A @ np.asarray(mmap[1], float)
        M = (P.T - np.eye(n)).copy()
        M[n - 1, :] = 1.0
        rhs = np.zeros(n)
        rhs[n - 1] = 1.0
        pie = np.linalg.solve(M, rhs)
        one = np.ones(n)
        out = []
        for c in range(K):
            Dc = np.asarray(mmap[2 + c], float)
            p = pie @ A @ Dc @ one
            out.append((p, (pie @ A @ Dc @ A @ one) / p,
                        (pie @ A @ A @ Dc @ one) / p,
                        (pie @ np.linalg.matrix_power(A, 3) @ Dc @ one) / p))
        return np.array(out)

    tgt = chars(src)
    fit = mmap3k_fit(D0, D1, tgt[:, 0], tgt[:, 1], tgt[:, 2],
                     tgt[:, 3] if z > 3 else None, exact_only=True)
    got = chars(fit)
    for k in range(4):
        assert got[:, k] == pytest.approx(tgt[:, k], rel=1e-9)
    assert validate_mmap(fit)
