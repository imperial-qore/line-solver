"""Regression tests for pfqn_sens (analytic sensitivities) and pfqn_momlin.

pfqn_sens returns exact derivatives of the mean measures {X,Q,U,R} of a closed
product-form network with respect to the service demands L(i,r) and the think
times Z(r). It dispatches between two internal kernels (sens.py:108):

  * CoMoM      -- selected for the repairman model (M == 1, all mi == 1, and
                  every populated class with Z > FineTol and L > FineTol);
  * diff. MVA  -- every other model (M > 1, zero demand/think time, mi != 1).

Because the dispatch is automatic, the two kernels have no overlapping public
entry point: the repairman regime is served only by CoMoM and everything else
only by MVA. The tests below pin their agreement on the overlap by calling the
private kernels directly, and independently referee both against central finite
differences of pfqn_mva, which is oblivious to which kernel produced a value.

Golden values are the MATLAB pfqn_sens output (matlab/src/api/pfqn/pfqn_sens.m),
which jline.api.pfqn.sens.Pfqn_sens reproduces; see test_pfqn_sens.m in
line-test.git for the MATLAB-side counterpart.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_sens, pfqn_momlin, pfqn_mva
from line_solver.api.pfqn.sens import _sens_comom, _sens_mva


# Repairman models inside the CoMoM/MVA overlap: M=1, Z>0, L>0.
REPAIRMAN_CASES = [
    (np.array([[1.0]]), np.array([5.0]), np.array([2.0])),
    (np.array([[1.0, 2.0]]), np.array([3.0, 4.0]), np.array([1.0, 3.0])),
    (np.array([[0.5, 1.5, 2.5]]), np.array([2.0, 3.0, 4.0]), np.array([1.0, 2.0, 0.5])),
    (np.array([[3.0, 0.7]]), np.array([8.0, 6.0]), np.array([0.4, 5.0])),
]


def _relerr(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    den = np.maximum(np.abs(a), np.abs(b))
    return np.max(np.abs(a - b) / np.where(den < 1e-12, 1.0, den))


def _mva_base(L, N, Z):
    """Base X and Q from pfqn_mva, which returns a 7-tuple (X at 0, Q at 2)."""
    r = pfqn_mva(np.asarray(L, float), np.asarray(N, float), np.asarray(Z, float))
    return np.asarray(r[0]).flatten(), np.asarray(r[2]).reshape(np.shape(L))


def _fd_jacobian(L, N, Z, h=1e-6):
    """Central differences of the pfqn_mva base measures wrt each L(i,r), Z(r).

    Parameter order mirrors pfqn_sens: every L(i,r) first, then every Z(r).
    """
    M, R = L.shape
    dX, dQ = [], []
    for i in range(M):
        for r in range(R):
            Lp, Lm = L.copy(), L.copy()
            Lp[i, r] += h
            Lm[i, r] -= h
            Xp, Qp = _mva_base(Lp, N, Z)
            Xm, Qm = _mva_base(Lm, N, Z)
            dX.append((Xp - Xm) / (2 * h))
            dQ.append((Qp - Qm) / (2 * h))
    for r in range(R):
        Zp, Zm = Z.copy(), Z.copy()
        Zp[r] += h
        Zm[r] -= h
        Xp, Qp = _mva_base(L, N, Zp)
        Xm, Qm = _mva_base(L, N, Zm)
        dX.append((Xp - Xm) / (2 * h))
        dQ.append((Qp - Qm) / (2 * h))
    return np.array(dX).T, np.stack(dQ, axis=-1)


# ---------------------------------------------------------------- dispatch

def test_dispatch_selects_comom_for_repairman():
    """M=1 with positive L and Z must be served by the CoMoM kernel."""
    L, N, Z = np.array([[1.0, 2.0]]), np.array([3.0, 4.0]), np.array([1.0, 3.0])
    assert _relerr(pfqn_sens(L, N, Z).Q, _sens_comom(L, N, Z).Q) < 1e-12


def test_dispatch_falls_back_to_mva_on_zero_think_time():
    """A populated class with Z=0 leaves the CoMoM regime, so MVA must serve it."""
    L, N, Z = np.array([[1.0, 2.0]]), np.array([3.0, 4.0]), np.array([0.0, 3.0])
    assert _relerr(pfqn_sens(L, N, Z).Q, _sens_mva(L, N, Z, np.ones(1)).Q) < 1e-12


# ------------------------------------------------- CoMoM vs MVA on overlap

@pytest.mark.parametrize("L,N,Z", REPAIRMAN_CASES)
def test_comom_matches_mva_base_measures(L, N, Z):
    sc, sm = _sens_comom(L, N, Z), _sens_mva(L, N, Z, np.ones(L.shape[0]))
    assert _relerr(sc.X, sm.X) < 1e-12
    assert _relerr(sc.Q, sm.Q) < 1e-12
    assert _relerr(sc.U, sm.U) < 1e-12
    assert _relerr(sc.R, sm.R) < 1e-12


@pytest.mark.parametrize("L,N,Z", REPAIRMAN_CASES)
def test_comom_matches_mva_jacobians(L, N, Z):
    sc, sm = _sens_comom(L, N, Z), _sens_mva(L, N, Z, np.ones(L.shape[0]))
    # Loosest observed residual is ~1e-11 on the badly conditioned Z=0.01 case;
    # 1e-8 leaves headroom without admitting a genuine kernel disagreement.
    assert _relerr(sc.dX, sm.dX) < 1e-8
    assert _relerr(sc.dQ, sm.dQ) < 1e-8
    assert _relerr(sc.dU, sm.dU) < 1e-8
    assert _relerr(sc.dR, sm.dR) < 1e-8


def test_comom_matches_mva_via_zero_demand_padding():
    """Padding a zero-demand station defeats the M==1 dispatch without changing
    the physics, so the public entry point must return the CoMoM answer at
    station 0 and an idle padded station."""
    L, N, Z = np.array([[1.0, 2.0]]), np.array([3.0, 4.0]), np.array([1.0, 3.0])
    comom = pfqn_sens(L, N, Z)
    padded = pfqn_sens(np.vstack([L, np.zeros((1, L.shape[1]))]), N, Z)
    assert _relerr(comom.X, padded.X) < 1e-8
    assert _relerr(np.asarray(comom.Q)[0, :], np.asarray(padded.Q)[0, :]) < 1e-8
    np.testing.assert_allclose(np.asarray(padded.Q)[1, :], 0.0, atol=1e-12)


# ------------------------------------------------ finite-difference referee

@pytest.mark.parametrize("L,N,Z", REPAIRMAN_CASES[:2])
def test_both_kernels_match_finite_differences(L, N, Z):
    """Neutral oracle: both kernels must track central differences of pfqn_mva.

    Tolerance is set by the O(h^2) truncation error of the difference scheme,
    not by the kernels, which is why both are held to the same bound.
    """
    fdX, fdQ = _fd_jacobian(L, N, Z)
    for sens in (_sens_comom(L, N, Z), _sens_mva(L, N, Z, np.ones(L.shape[0]))):
        assert _relerr(sens.dX, fdX) < 1e-6
        assert _relerr(sens.dQ, fdQ) < 1e-6


def test_mva_kernel_matches_finite_differences_multistation():
    """The MVA kernel owns every M>1 model, which CoMoM never sees."""
    L = np.array([[1.0, 2.0], [0.5, 0.3]])
    N, Z = np.array([3.0, 4.0]), np.array([1.0, 3.0])
    fdX, fdQ = _fd_jacobian(L, N, Z)
    sens = pfqn_sens(L, N, Z)
    assert _relerr(sens.dX, fdX) < 1e-6
    assert _relerr(sens.dQ, fdQ) < 1e-6


# --------------------------------------------------------- parameter list

def test_param_descriptors_order_and_content():
    L, N, Z = np.array([[1.0, 2.0]]), np.array([3.0, 4.0]), np.array([1.0, 3.0])
    params = pfqn_sens(L, N, Z).params
    assert len(params) == 1 * 2 + 2  # M*R demands, then R think times
    assert [p["type"] for p in params] == ["L", "L", "Z", "Z"]
    assert [p["jobclass"] for p in params] == [0, 1, 0, 1]
    # Native Python is 0-based and marks think-time parameters with station -1.
    # MATLAB uses field 'class' and station 0 for Z; see test_pfqn_sens.m.
    assert [p["station"] for p in params] == [0, 0, -1, -1]


# ------------------------------------------------- second moments (QCov/QVar)

def test_qvar_is_diagonal_of_qcov():
    L = np.array([[1.0, 2.0], [0.5, 0.3]])
    sens = pfqn_sens(L, np.array([3.0, 4.0]), np.array([1.0, 3.0]))
    QCov, QVar = np.asarray(sens.QCov), np.asarray(sens.QVar)
    M, R = L.shape
    diag = np.array([[QCov[i, r, i, r] for r in range(R)] for i in range(M)])
    np.testing.assert_allclose(diag, QVar, rtol=0, atol=0)


def test_qcov_is_symmetric_and_variance_nonnegative():
    L = np.array([[1.0, 2.0], [0.5, 0.3]])
    sens = pfqn_sens(L, np.array([3.0, 4.0]), np.array([1.0, 3.0]))
    QCov = np.asarray(sens.QCov)
    M, R = L.shape
    flat = QCov.reshape(M * R, M * R)
    np.testing.assert_allclose(flat, flat.T, rtol=0, atol=1e-14)
    assert np.all(np.asarray(sens.QVar) >= 0.0)


def test_qvar_is_bernoulli_at_unit_population():
    """With N=1 the per-station queue length is Bernoulli, so Var = q(1-q).

    An analytic check independent of both kernels and of finite differences.
    """
    L = np.array([[1.0], [0.5]])
    sens = pfqn_sens(L, np.array([1.0]), np.array([0.0]))
    q = np.asarray(sens.Q).flatten()
    np.testing.assert_allclose(np.asarray(sens.QVar).flatten(), q * (1.0 - q),
                               rtol=1e-12, atol=1e-14)


# ------------------------------------------------------- MATLAB parity golden

def test_matches_matlab_golden_repairman():
    """Locked to MATLAB pfqn_sens on L=[1 2], N=[3 4], Z=[1 3] (CoMoM path)."""
    sens = pfqn_sens(np.array([[1.0, 2.0]]), np.array([3.0, 4.0]),
                     np.array([1.0, 3.0]))
    np.testing.assert_allclose(np.asarray(sens.X).flatten(),
                               [0.449063352907, 0.275328756927], rtol=1e-11)
    np.testing.assert_allclose(np.asarray(sens.Q).flatten(),
                               [2.550936647093, 3.174013729220], rtol=1e-11)
    np.testing.assert_allclose(np.asarray(sens.dX)[0, 0], -0.407069671034, rtol=1e-11)
    np.testing.assert_allclose(np.asarray(sens.dQ)[0, 0, 0], 0.407069671034, rtol=1e-11)
    np.testing.assert_allclose(np.asarray(sens.QVar).flatten(),
                               [0.407069671034, 0.731773642686], rtol=1e-11)


# -------------------------------------------------------------- pfqn_momlin

def test_momlin_is_exact_at_unit_population():
    """Schweitzer-Bard is exact for N=1, so the linearizer must be too."""
    L, N, Z = np.array([[1.0], [0.5]]), np.array([1.0]), np.array([0.0])
    m = pfqn_momlin(L, N, Z)
    Xe, Qe = _mva_base(L, N, Z)
    np.testing.assert_allclose(np.asarray(m.Q).reshape(L.shape), Qe, rtol=1e-8)
    np.testing.assert_allclose(np.asarray(m.X).flatten(), Xe, rtol=1e-8)


def test_momlin_approximates_exact_mva():
    """momlin carries AMVA error and is NOT exact: assert it stays in band.

    Observed Schweitzer-Bard error on this model is ~11% on Q and ~4% on X;
    the bounds below are deliberately loose because the point is to catch a
    broken linearizer, not to pin an approximation to exact values.
    """
    L = np.array([[1.0, 2.0], [0.5, 0.3]])
    N, Z = np.array([3.0, 4.0]), np.array([1.0, 3.0])
    m = pfqn_momlin(L, N, Z)
    Xe, Qe = _mva_base(L, N, Z)
    assert _relerr(np.asarray(m.Q).reshape(L.shape), Qe) < 0.20
    assert _relerr(np.asarray(m.X).flatten(), Xe) < 0.20


def test_momlin_error_decreases_with_population():
    """Schweitzer-Bard is asymptotically exact, so error must shrink as N grows."""
    L = np.array([[1.0, 2.0], [0.5, 0.3]])
    Z = np.array([1.0, 3.0])
    small, big = np.array([3.0, 4.0]), np.array([20.0, 15.0])
    err = []
    for N in (small, big):
        _, Qe = _mva_base(L, N, Z)
        err.append(_relerr(np.asarray(pfqn_momlin(L, N, Z).Q).reshape(L.shape), Qe))
    assert err[1] < err[0]


def test_momlin_conserves_population():
    L = np.array([[1.0, 2.0], [0.5, 0.3]])
    N, Z = np.array([3.0, 4.0]), np.array([1.0, 3.0])
    m = pfqn_momlin(L, N, Z)
    Q = np.asarray(m.Q).reshape(L.shape)
    X = np.asarray(m.X).flatten()
    # Queued jobs plus those thinking must account for the whole population.
    np.testing.assert_allclose(Q.sum(axis=0) + X * Z, N, rtol=1e-6)


def test_momlin_second_moments_track_exact_sensitivities():
    """momlin QVar uses the same product-form identity as pfqn_sens, but on an
    AMVA fixed point, so it approximates rather than reproduces the exact QVar."""
    L = np.array([[1.0, 2.0], [0.5, 0.3]])
    N, Z = np.array([3.0, 4.0]), np.array([1.0, 3.0])
    m = pfqn_momlin(L, N, Z)
    exact = pfqn_sens(L, N, Z)
    assert np.all(np.asarray(m.QVar) >= 0.0)
    assert _relerr(np.asarray(m.QVar), np.asarray(exact.QVar)) < 0.30
