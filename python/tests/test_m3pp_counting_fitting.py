"""Drift guard for the native-Python M3PP counting-process fitters.

``m3pp2m_fitc_approx``, ``m3pp2m_fitc_approx_ag_multiclass``, the trace and
superposition entry points and ``m3pp22_interleave_fitc`` are ports of the
MATLAB m3a files of the same names. What is asserted is the algebra the fits
are built on, not a generator taken from another implementation:

- THE COVARIANCE SPLIT IS AN EXACT INVERSION. Handed the count covariance of a
  genuine M3PP(2,2), ``m3pp22_fitc_approx_cov_multiclass`` must return that
  process's own per-phase marking probabilities, not an approximation.
- THE PER-CLASS TARGETS ARE RECOVERED when they are those of an actual
  M3PP(2,m): the least-squares split then reaches its own targets exactly.
- A MARKING IS A PARTITION: sum_c Dc = D1 entry by entry.
- THE RATES ARE MATCHED EXACTLY in every fitter of the family, by construction.
- THE ORDERS ARE STRUCTURAL: superposing L two-phase processes gives the
  product chain of order 2^L, interleaving lumps them onto L + 1 levels.

The Poisson regime is out of the family's domain: IDC is then 1 at every scale,
the fitters return a reducible MMPP(2), and MATLAB's ``ctmc_solve`` refuses it
with the same message the native one does.

See _kb/03-api-layer.md, "M3PP counting-process fitters".
"""
import numpy as np
import pytest

from line_solver.api.mam import (map_count_mean, mmap_count_mcov, mmap_count_mean,
                                 mmap_count_moment, mmap_count_var, mmap_isfeasible)
from line_solver.api.trace import mtrace_iat2counts
from line_solver.lib.kpctoolbox import mmpp2_fitc_approx
from line_solver.lib.m3a.m3pp import (_m3pp_split_coeffs_dv, _m3pp_split_solve,
                                      m3pp22_fitc_approx_cov_multiclass,
                                      m3pp22_interleave_fitc, m3pp2m_fitc_approx_ag_multiclass,
                                      m3pp2m_fitc_trace, m3pp2m_interleave,
                                      m3pp_superpos_fitc, m3pp_superpos_fitc_theoretical,
                                      m3pp_superpos_fitc_trace)
from line_solver.lib.m3a.fit import M3aFitOptions, m3a_fit, m3afit_init


def _genuine_m3pp22(l1=2.0, l2=0.5, r1=0.3, r2=0.7, q1=0.7, q2=0.3):
    """An MMPP(2) with a per-phase Bernoulli marking, i.e. an exact M3PP(2,2)."""
    D0 = np.array([[-(l1 + r1), r1], [r2, -(l2 + r2)]])
    D1 = np.array([[l1, 0.0], [0.0, l2]])
    return [D0, D1,
            D1 * np.array([[q1, 0.0], [0.0, q2]]),
            D1 * np.array([[1 - q1, 0.0], [0.0, 1 - q2]])]


def _marking_defect(m):
    s = np.zeros_like(np.asarray(m[1], float))
    for c in range(len(m) - 2):
        s = s + np.asarray(m[2 + c], float)
    return float(np.max(np.abs(s - np.asarray(m[1], float))))


def _bursty_trace(n=20000, seed=3):
    """MMPP(2)-modulated two-class trace: the fast phase favours class 1."""
    rng = np.random.default_rng(seed)
    lam = [4.0, 0.4]
    sw = [0.05, 0.05]
    T = np.zeros(n)
    C = np.zeros(n, dtype=int)
    ph = 0
    for i in range(n):
        te = rng.exponential(1.0 / lam[ph])
        ts = rng.exponential(1.0 / sw[ph])
        if ts < te:
            ph = 1 - ph
            T[i] = ts
        else:
            T[i] = te
        C[i] = 1 if rng.random() < (0.3 if ph == 0 else 0.7) else 2
    return T, C


def test_mmap_count_moment_matches_the_poisson_closed_form():
    a, t = 3.0, 0.7
    p = np.array([0.25, 0.75])
    MP = [np.array([[-a]]), np.array([[a]]),
          np.array([[p[0] * a]]), np.array([[p[1] * a]])]
    M = mmap_count_moment(MP, t, np.array([1, 2, 3]))
    lam = p * a * t
    assert np.allclose(M[0], lam, rtol=1e-6)
    assert np.allclose(M[1], lam + lam ** 2, rtol=1e-6)
    assert np.allclose(M[2], lam + 3 * lam ** 2 + lam ** 3, rtol=1e-6)


def test_mtrace_iat2counts_anchors_every_window_at_an_arrival():
    # unit inter-arrivals, alternating labels; a window of 2.0 covers the next
    # two arrivals, one of each class, until it runs off the end of the trace
    T = np.ones(6)
    A = np.array([1, 2, 1, 2, 1, 2])
    C = mtrace_iat2counts(T, A, 2.0)
    assert C.shape == (4, 2)
    assert np.array_equal(C, np.ones((4, 2), dtype=C.dtype))


def test_covariance_split_inverts_its_own_covariance():
    truth = _genuine_m3pp22()
    t3 = 1.0
    ai = np.ravel(mmap_count_mean(truth, 1.0))
    st3 = mmap_count_mcov(truth, t3)[0, 1]

    F = m3pp22_fitc_approx_cov_multiclass([truth[0], truth[1]], ai, st3, t3)
    assert F[2][0, 0] / truth[1][0, 0] == pytest.approx(0.7, rel=1e-10)
    assert F[2][1, 1] / truth[1][1, 1] == pytest.approx(0.3, rel=1e-10)
    assert mmap_count_mcov(F, t3)[0, 1] == pytest.approx(st3, rel=1e-10)
    assert _marking_defect(F) < 1e-12
    assert np.allclose(np.ravel(mmap_count_mean(F, 1.0)), ai, rtol=1e-12)


def test_ag_split_recovers_the_targets_of_an_actual_m3pp():
    truth = _genuine_m3pp22()
    t3 = 1.0
    ai = np.ravel(mmap_count_mean(truth, 1.0))
    V = np.ravel(mmap_count_var(truth, t3))
    S = mmap_count_mcov(truth, t3)
    gt3 = np.array([V[i] + sum(S[i, j] for j in range(2) if j != i) for i in range(2)])

    F = m3pp2m_fitc_approx_ag_multiclass([truth[0], truth[1]], ai, gt3, t3)
    assert _marking_defect(F) < 1e-12
    assert np.allclose(np.ravel(mmap_count_mean(F, 1.0)), ai, rtol=1e-10)
    assert mmap_isfeasible(F)

    Vf = np.ravel(mmap_count_var(F, t3))
    Sf = mmap_count_mcov(F, t3)
    gf = np.array([Vf[i] + sum(Sf[i, j] for j in range(2) if j != i) for i in range(2)])
    assert np.allclose(gf, gt3, rtol=1e-8)


def test_dv_split_recovers_the_targets_of_an_actual_m3pp():
    truth = _genuine_m3pp22()
    t3 = 1.0
    ai = np.ravel(mmap_count_mean(truth, 1.0))
    dv = np.zeros(2)
    for i in range(2):
        two = [truth[0], truth[1], truth[2 + i], truth[1] - truth[2 + i]]
        V = np.ravel(mmap_count_var(two, t3))
        dv[i] = V[0] - V[1]
    a = float(np.ravel(map_count_mean(truth[0], truth[1], 1.0))[0])

    cf = _m3pp_split_coeffs_dv(truth[1][0, 0], truth[1][1, 1],
                               truth[0][0, 1], truth[0][1, 0], t3)
    x = _m3pp_split_solve(cf, ai, dv, a)
    assert np.allclose(x, dv, rtol=1e-8)

    q1 = np.array([cf[0] * ai[i] + cf[1] * x[i] + cf[2] for i in range(2)])
    q2 = np.array([cf[3] * ai[i] + cf[4] * x[i] + cf[5] for i in range(2)])
    assert q1.sum() == pytest.approx(1.0, abs=1e-10)
    assert q2.sum() == pytest.approx(1.0, abs=1e-10)


def test_mmpp2_fitc_approx_matches_the_requested_rate_exactly():
    from line_solver.api.mam import map_lambda
    for a in (0.5, 2.0):
        D0, D1 = mmpp2_fitc_approx(a, 3.0, 3.2, 8.0, 5.0, 1.0, 2.0)
        assert map_lambda(D0, D1) == pytest.approx(a, rel=1e-10)


def test_superposition_gives_the_product_chain_and_matches_every_rate():
    av = np.array([1.0, 2.0])
    fit, parts = m3pp_superpos_fitc(av, np.array([2.0, 3.0]), np.array([5.0, 6.0]),
                                    np.array([0.5, 0.8]), 1.0, 100.0)
    assert len(parts) == 2
    assert fit[0].shape[0] == 4          # 2 x 2, the product chain
    assert len(fit) - 2 == 2
    assert np.allclose(np.ravel(mmap_count_mean(fit, 1.0)), av, rtol=1e-10)
    assert _marking_defect(fit) < 1e-12


def test_superposition_theoretical_reproduces_the_class_rates_of_its_target():
    truth = _genuine_m3pp22()
    want = np.ravel(mmap_count_mean(truth, 1.0))
    fit, _ = m3pp_superpos_fitc_theoretical(truth, 1.0, 100.0)
    assert fit[0].shape[0] == 4
    assert np.allclose(np.ravel(mmap_count_mean(fit, 1.0)), want, rtol=1e-10)


def test_interleaving_lumps_onto_one_level_per_component_plus_one():
    parts = [_genuine_m3pp22(2.0, 0.5, 0.9, 0.4, 0.7, 0.3),
             _genuine_m3pp22(3.0, 0.8, 1.4, 0.9, 0.4, 0.6)]
    s = m3pp2m_interleave(parts)
    assert s[0].shape[0] == 3
    assert len(s) - 2 == 4
    assert _marking_defect(s) < 1e-12
    rows = np.sum(np.asarray(s[0], float) + np.asarray(s[1], float), axis=1)
    assert np.allclose(rows, 0.0, atol=1e-12)


def test_interleave_fitc_matches_every_per_class_rate():
    av = np.array([[0.6, 0.4], [1.2, 0.8]])
    lump, parts = m3pp22_interleave_fitc(av, np.array([2.0, 2.5]), np.array([6.0, 7.0]),
                                         np.array([1.0, 1.0]), 1.0)
    assert len(parts) == 2
    assert lump[0].shape[0] == 3
    assert len(lump) - 2 == 4
    assert _marking_defect(lump) < 1e-12
    assert np.allclose(np.ravel(mmap_count_mean(lump, 1.0)), av.ravel(), rtol=1e-10)
    assert mmap_isfeasible(lump)


def test_interleave_fitc_refuses_an_infeasible_idc_pair():
    av = np.array([[0.5, 0.5]])
    with pytest.raises(ValueError):
        # IDC(t) below one is sub-Poisson, outside the family
        m3pp22_interleave_fitc(av, np.array([0.5]), np.array([3.0]), np.array([0.0]), 1.0)


def test_trace_fitters_read_the_rate_and_the_class_split_off_the_trace():
    T, C = _bursty_trace()
    a = 1.0 / float(np.mean(T))
    p1 = float(np.sum(C == 1)) / len(C)

    for method in ('approx_ag', 'approx_delta'):
        F = m3pp2m_fitc_trace(T, C, method)
        r = np.ravel(mmap_count_mean(F, 1.0))
        assert r[0] == pytest.approx(a * p1, rel=1e-6)
        assert r[1] == pytest.approx(a * (1.0 - p1), rel=1e-6)
        assert mmap_isfeasible(F)

    fit, _ = m3pp_superpos_fitc_trace(T, C)
    r = np.ravel(mmap_count_mean(fit, 1.0))
    assert r[0] == pytest.approx(a * p1, rel=1e-6)
    assert r[1] == pytest.approx(a * (1.0 - p1), rel=1e-6)

    with pytest.raises(ValueError):
        m3pp2m_fitc_trace(T, C, 'nonesuch')


def test_m3a_fit_counting_branches_reach_the_real_fitters():
    T, C = _bursty_trace()
    a = 1.0 / float(np.mean(T))
    mt = m3afit_init(T, C)

    # NumStates == 2 routes to m3pp2m_fitc_trace('approx_ag')
    M2 = m3a_fit(mt, M3aFitOptions(method=1, num_states=2))
    assert M2 is not None and M2[0].shape[0] == 2
    assert float(np.sum(mmap_count_mean(M2, 1.0))) == pytest.approx(a, rel=1e-6)

    # NumStates > 2 routes to m3pp_superpos_fitc_trace, whose order is set by
    # the superposition, not by the requested state count
    M3 = m3a_fit(mt, M3aFitOptions(method=1, num_states=3))
    assert M3 is not None and M3[0].shape[0] >= 2
    assert float(np.sum(mmap_count_mean(M3, 1.0))) == pytest.approx(a, rel=1e-6)
