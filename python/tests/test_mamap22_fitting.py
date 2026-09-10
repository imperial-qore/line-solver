"""Drift guard for the native-Python mamap22 fitters.

``mamap22_fit_fs_multiclass``, ``mamap22_fit_bs_multiclass`` and the six
``mamap22_fit_gamma_*`` drivers are ports of the MATLAB m3a files of the same
names. Assertions are of three kinds:

- ROUND TRIP: a marking is read back off the MMAP it produced, so the exact
  branch is inverted against its own forward map.
- PINNED VALUES: D11 taken from MATLAB for the same input, so a change that
  stays self-consistent but is no longer the same algorithm is caught. The
  reference was run with the YALMIP options moved into the repair branch, which
  is what lets the exact path run without a YALMIP installation.
- REFUSALS BY NAME: the nonconvex bmibnb repair is unported and must say so.

See _kb/03-api-layer.md, "Closed-form fitter audit".
"""
import numpy as np
import pytest

from line_solver.api.mam.map_analysis import map_normalize
from line_solver.api.mam.mmap_ops import (mmap_backward_moment, mmap_forward_moment,
                                          mmap_pc, mmap_sigma)
from line_solver.lib.m3a.mamap22 import (Mamap22Unsupported, mamap22_fit_bs_multiclass,
                                         mamap22_fit_fs_multiclass,
                                         mamap22_fit_gamma_bs_mmap,
                                         mamap22_fit_gamma_bs_trace,
                                         mamap22_fit_gamma_fs_mmap,
                                         mamap22_fit_gamma_fs_trace,
                                         mamap2m_can1_coefficients,
                                         mamap2m_can2_coefficients)


def _amap(h1, h2, r1, r2, form):
    """The AMAP(2) the reference's own Poisson-perturbation block reconstructs.

    NOTE the convention: that block writes `r2/h2` into D1(2,1) and `(1-r2)/h2`
    into D1(2,2), while the fitters READ `r2 = D1(2,2)*h2`. So the r2 passed here
    is the COMPLEMENT of the r2 the fitters see. The pinned cases are built this
    way because that is what the MATLAB reference was run on; a test that needs a
    particular READ r2 must pass 1 - r2 here.
    """
    D0 = np.array([[-1.0 / h1, r1 / h1], [0.0, -1.0 / h2]])
    if form == 1:
        D1 = np.array([[(1 - r1) / h1, 0.0], [r2 / h2, (1 - r2) / h2]])
    else:
        D1 = np.array([[0.0, (1 - r1) / h1], [r2 / h2, (1 - r2) / h2]])
    return list(map_normalize(D0, D1))


def _marked(h1, h2, r1, r2, form, q1, q2, q3):
    D0, D1 = _amap(h1, h2, r1, r2, form)
    if form == 1:
        m1 = np.array([[q1, 0.0], [q2, q3]])
        m2 = np.array([[1 - q1, 0.0], [1 - q2, 1 - q3]])
    else:
        m1 = np.array([[0.0, q1], [q2, q3]])
        m2 = np.array([[0.0, 1 - q1], [1 - q2, 1 - q3]])
    return [D0, D1, D1 * m1, D1 * m2]


def _targets(mmap):
    return (np.asarray(mmap_pc(mmap), float).ravel(),
            np.asarray(mmap_forward_moment(mmap, [1]), float).ravel(),
            np.asarray(mmap_backward_moment(mmap, [1]), float).ravel(),
            np.atleast_2d(np.asarray(mmap_sigma(mmap), float)))


def _is_mmap(mmap):
    D0, D1 = np.asarray(mmap[0], float), np.asarray(mmap[1], float)
    assert np.allclose(D0 + D1, sum(np.asarray(m, float) for m in mmap[2:]) + D0, atol=1e-12)
    assert np.allclose((D0 + D1).sum(axis=1), 0.0, atol=1e-10)
    for m in mmap[2:]:
        assert np.all(np.asarray(m, float) >= -1e-14)


# (case, h1, h2, r1, r2, form, q1, q2, q3) and the MATLAB D11, row-major
PINNED = [
    (1, 1.0, 3.0, 0.4, 0.3, 1, 0.7, 0.2, 0.5,
     [0.3, 0.0, 0.05, 0.116666666666667],
     [0.42, 0.0, 0.0199999999999995, 0.116666666666667]),
    (2, 1.0, 3.0, 0.4, 0.3, 2, 0.7, 0.2, 0.5,
     [0.0, 0.419999999999957, 0.02, 0.116666666666671],
     [0.0, 0.42, 0.0200000000000004, 0.116666666666666]),
    (3, 0.5, 2.0, 0.6, 0.15, 1, 0.3, 0.8, 0.45,
     [0.239999999999466, 0.0, 0.0600000000000336, 0.19125],
     [0.240000000000001, 0.0, 0.0600000000000014, 0.191249999999999]),
]


@pytest.mark.parametrize('case,h1,h2,r1,r2,form,q1,q2,q3,d11_fs,d11_bs', PINNED)
def test_exact_fit_matches_matlab(case, h1, h2, r1, r2, form, q1, q2, q3, d11_fs, d11_bs):
    src = _marked(h1, h2, r1, r2, form, q1, q2, q3)
    p, F, B, S = _targets(src)

    mfs, fF, fSf, exact_fs = mamap22_fit_fs_multiclass([src[0], src[1]], p, F, S)
    assert exact_fs
    _is_mmap(mfs)
    assert np.allclose(np.asarray(mfs[2], float).ravel(), d11_fs, atol=1e-9)
    assert np.allclose(fF, F, atol=1e-9)
    assert fSf[0, 0] == pytest.approx(S[0, 0], abs=1e-9)

    mbs, fB, fSb, exact_bs = mamap22_fit_bs_multiclass([src[0], src[1]], p, B, S)
    assert exact_bs
    _is_mmap(mbs)
    assert np.allclose(np.asarray(mbs[2], float).ravel(), d11_bs, atol=1e-9)
    assert np.allclose(fB, B, atol=1e-9)
    assert fSb[0, 0] == pytest.approx(S[0, 0], abs=1e-9)
    # the class probabilities are matched EXACTLY by both, by construction
    assert np.allclose(np.asarray(mmap_pc(mfs), float).ravel(), p, atol=1e-10)
    assert np.allclose(np.asarray(mmap_pc(mbs), float).ravel(), p, atol=1e-10)


def test_gamma_drivers_recover_the_source_marking():
    src = _marked(0.5, 2.0, 0.6, 0.15, 1, 0.3, 0.8, 0.45)
    for driver in (mamap22_fit_gamma_fs_mmap, mamap22_fit_gamma_bs_mmap):
        out = driver(src)
        _is_mmap(out)
        assert np.allclose(np.asarray(out[2], float), np.asarray(src[2], float), atol=1e-8)


def test_coefficient_tables_have_no_hole():
    # G(9) is the entry the reference discards by assigning G(10) twice; the
    # port fills it, and Y must not read it either way.
    G, U, Y = mamap2m_can1_coefficients(1.0, 3.0, 0.4, 0.3)
    assert G.shape == (15,) and U.shape == (12,) and Y.shape == (3,)
    assert G[8] == pytest.approx(0.4 * 0.3 ** 2 / (0.4 * 0.3 - 0.3 + 1.0))
    assert G[9] == pytest.approx(1.0 - 1.0 * 0.4 / (0.3 * (0.4 - 1.0) + 1.0))
    E, V, Z = mamap2m_can2_coefficients(1.0, 3.0, 0.4, 0.3)
    assert E.shape == (14,) and V.shape == (12,) and Z.shape == (3,)


def test_the_nonconvex_repair_is_refused_by_name():
    a = _amap(2.0, 0.6, 0.5, 0.4, 1)
    p = np.array([0.5, 0.5])
    F = np.array([50.0, 50.0])
    S = np.full((2, 2), 0.9)
    with pytest.raises(Mamap22Unsupported, match='bmibnb'):
        mamap22_fit_fs_multiclass(a, p, F, S)
    with pytest.raises(Mamap22Unsupported, match='bmibnb'):
        mamap22_fit_bs_multiclass(a, p, F, S)
    # declining the repair returns the clamped closed form as a valid MMAP
    out, _, _, exact = mamap22_fit_fs_multiclass(a, p, F, S, None, None, False)
    assert not exact
    _is_mmap(out)


def test_the_sigma_arm_of_the_gamma_negative_degeneracy_is_ported():
    # read r2 = 0 in form 2 (so 1.0 here, see _amap), sigma weighted above the
    # moment: the reference states this repair as a YALMIP program, but its
    # feasible set is an interval, so the projection is its global optimum.
    b = _amap(2.0, 0.6, 0.5, 1.0, 2)
    p = np.array([0.5, 0.5])
    F = np.array([50.0, 50.0])
    S = np.full((2, 2), 0.9)   # above the p1^2 = 0.25 ceiling, so it clamps
    out, _, fS, _ = mamap22_fit_fs_multiclass(b, p, F, S, None, [1.0, 2.0])
    _is_mmap(out)
    assert fS[0, 0] == pytest.approx(0.25, abs=1e-9)


def test_two_classes_only_and_the_canonical_form_are_enforced():
    a = _amap(2.0, 0.6, 0.5, 0.4, 1)
    S = np.full((2, 2), 0.25)
    with pytest.raises(ValueError):
        mamap22_fit_fs_multiclass(a, np.full(3, 1 / 3), np.ones(3), S)
    with pytest.raises(ValueError):
        mamap22_fit_bs_multiclass(a, np.full(3, 1 / 3), np.ones(3), S)
    non_canonical = [np.array([[-1.0, 0.5], [0.0, -2.0]]),
                     np.array([[0.3, 0.2], [1.0, 1.0]])]
    with pytest.raises(ValueError, match='canonical'):
        mamap22_fit_fs_multiclass(non_canonical, np.array([0.5, 0.5]),
                                  np.ones(2), S)
    cyclic = [np.array([[-1.0, 0.5], [0.1, -2.0]]),
              np.array([[0.5, 0.0], [1.0, 0.9]])]
    with pytest.raises(ValueError, match='acyclic'):
        mamap22_fit_bs_multiclass(cyclic, np.array([0.5, 0.5]), np.ones(2), S)


@pytest.mark.parametrize('seed', [11, 23, 101])
def test_the_trace_drivers_run_on_a_marked_map_path(seed):
    # The trace must be a MAMAP sample path, not a MAP path with independent
    # marks: the drivers fit F and sigma jointly, and an arbitrary (F, sigma)
    # pair need not be reachable by ANY AMAP(2) marking, in which case every
    # form refuses and only the unported bmibnb repair is left. Marking from the
    # DEPARTING PHASE keeps the targets inside the markable region.
    rng = np.random.default_rng(seed)
    D0 = np.array([[-2.0, 0.8], [0.0, -0.5]])
    D1 = np.array([[1.2, 0.0], [0.3, 0.2]])
    st, T, A = 0, [], []
    for _ in range(4000):
        tot = 0.0
        while True:
            w = np.concatenate([np.where(np.arange(2) == st, 0.0, D0[st, :]), D1[st, :]])
            w = np.clip(w, 0.0, None)
            tot += rng.exponential(1.0 / w.sum())
            j = rng.choice(4, p=w / w.sum())
            departing = st
            st = j if j < 2 else j - 2
            if j >= 2:
                break
        T.append(tot)
        A.append(1 if rng.random() < (0.75 if departing == 0 else 0.3) else 2)
    T, A = np.array(T), np.array(A)
    for driver in (mamap22_fit_gamma_fs_trace, mamap22_fit_gamma_bs_trace):
        try:
            out = driver(T, A)
        except ValueError as exc:
            # A sampled (F, sigma) pair need not be reachable by ANY AMAP(2)
            # marking; MATLAB reaches the same verdict on the same trace and
            # falls through to the unported bmibnb repair. Refusing by name is
            # the faithful outcome, so it is admitted -- but nothing else is.
            assert 'no AMAP(2) form admits' in str(exc)
            continue
        _is_mmap(out)
        assert np.asarray(out[0]).shape[0] == 2
        assert np.allclose(np.asarray(mmap_pc(out), float).ravel(),
                           [np.mean(A == 1), np.mean(A == 2)], atol=1e-8)
