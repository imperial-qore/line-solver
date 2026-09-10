"""
Native-Python ports of the joint-dependent CLW inversions, the Harel-Namn-Sturm
bounds and Xia's asymptotic load-dependent normalizing constant.

Every check is an identity against an independently computed value, never a
transcribed constant:
  1. pfqn_clwoi / pfqn_clwjd invert the generating function whose coefficient
     the convolution pfqn_ncoi / pfqn_ncjd computes, so the two must agree to
     the inversion accuracy.
  2. pfqn_clwjd with lcut = 1 IS pfqn_clwoi.
  3. Harel G(n) is the normalizing constant of the load-independent network, so
     TH(n) = G(n-1)/G(n) must equal pfqn_ca and the exact pfqn_mva throughput,
     and LB <= X_exact <= UB(n).
  4. pfqn_xia with one server per station is the dominant-pole asymptotics of
     the load-independent constant, available in closed form; with several
     servers it must converge to the exact pfqn_gld constant as N grows.
"""

import math

import numpy as np
import pytest

from line_solver.api.pfqn import (
    pfqn_ca, pfqn_clwjd, pfqn_clwoi, pfqn_gld, pfqn_harel_bounds,
    pfqn_harel_lb, pfqn_harel_ub, pfqn_mva, pfqn_ncjd, pfqn_ncoi, pfqn_xia,
)

INV_TOL = 1e-5      # lattice-Poisson aliasing at the CLW default parameters:
                    # l=[1,2,2] leaves the R=3 case at ~8e-7, and which side of
                    # 1e-6 it lands on is libm noise (7.6e-7 on glibc 2.35,
                    # 1.8e-6 on 2.31), so 1e-6 was not a portable bound.


def _busy_classes(n):
    return float(np.sum(np.asarray(n) > 0))


def _saturating(n):
    return 1.0 + float(np.sum(np.minimum(np.asarray(n), 2)))


@pytest.mark.parametrize("Z,N,mus,visits", [
    ([1.0, 2.0], [4, 3], [_busy_classes], None),
    ([0.5, 0.5], [3, 2], [lambda n: 1.0], None),
    ([1.0, 0.0], [3, 3], [lambda n: 1.0 + _busy_classes(n), lambda n: 2.0],
     [np.array([0.8, 1.2]), np.array([1.0, 0.5])]),
    ([0.7, 0.3, 1.1], [2, 3, 2], [lambda n: 1.0 + 0.5 * _busy_classes(n)], None),
])
def test_clwoi_inverts_the_ncoi_convolution(Z, N, mus, visits):
    G_conv = pfqn_ncoi(Z, N, mus, visits)[0]
    G_inv = pfqn_clwoi(Z, N, mus, visits)[0]
    assert abs(G_inv - G_conv) <= INV_TOL * abs(G_conv)


@pytest.mark.parametrize("Z,N,lcut", [
    ([1.0, 2.0], [8, 6], [2, 2]),
    ([0.5, 0.5], [4, 3], [2, 2]),
])
def test_clwjd_inverts_the_ncjd_convolution(Z, N, lcut):
    G_conv = pfqn_ncjd(Z, N, [_saturating])[0]
    G_inv = pfqn_clwjd(Z, N, [_saturating], None, lcut)[0]
    assert abs(G_inv - G_conv) <= INV_TOL * abs(G_conv)


def test_clwjd_multiserver_station():
    # mu(n) = min(sum n, c) is a function of the vector clipped at c
    c = 3
    rate = lambda n: float(min(np.sum(n), c)) if np.sum(n) > 0 else 1.0
    G_conv = pfqn_ncjd([1.0, 1.0], [5, 4], [rate])[0]
    G_inv = pfqn_clwjd([1.0, 1.0], [5, 4], [rate], None, [c, c])[0]
    assert abs(G_inv - G_conv) <= INV_TOL * abs(G_conv)


def test_clwjd_at_unit_cutoff_is_clwoi():
    rate = lambda n: 1.0 + _busy_classes(n)
    assert (pfqn_clwjd([1.0, 2.0], [4, 3], [rate], None, 1)[0] ==
            pytest.approx(pfqn_clwoi([1.0, 2.0], [4, 3], [rate])[0], rel=1e-15))


def test_clw_pure_delay_network():
    G = pfqn_clwoi([1.5, 0.5], [3, 2], None)[0]
    assert G == pytest.approx(1.5 ** 3 / math.factorial(3) * 0.5 ** 2 / math.factorial(2),
                              rel=1e-8)


def test_clw_refuses_a_rate_that_varies_inside_a_support():
    # a rate that is not support-only would be inverted to a plausible but
    # wrong G, so it must be an error rather than a warning
    bad = lambda n: 1.0 + float(np.sum(n))
    with pytest.raises(ValueError):
        pfqn_clwoi([1.0], [3], [bad])
    with pytest.raises(ValueError):
        pfqn_clwjd([1.0], [3], [bad], None, 1)


@pytest.mark.parametrize("rho", [
    [0.9, 0.4, 0.2], [0.7, 0.7, 0.3, 0.1], [0.55, 0.5, 0.45], [1.2, 0.3],
])
def test_harel_throughputs_are_the_exact_ones(rho):
    L = np.array(rho).reshape(-1, 1)
    _, _, TH = pfqn_harel_bounds(np.array(rho), 7)
    for n in range(1, 8):
        Gn = pfqn_ca(L, np.array([n]))[0]
        Gnm1 = pfqn_ca(L, np.array([n - 1]))[0]
        assert TH[n - 1] == pytest.approx(Gnm1 / Gn, rel=1e-10)
        X = float(np.atleast_1d(pfqn_mva(L, np.array([n]))[0]).ravel()[0])
        assert TH[n - 1] == pytest.approx(X, rel=1e-9)


@pytest.mark.parametrize("rho", [
    [0.9, 0.4, 0.2], [0.7, 0.7, 0.3, 0.1], [0.55, 0.5, 0.45], [1.2, 0.3],
])
@pytest.mark.parametrize("N", [2, 3, 5, 8, 12, 20])
def test_harel_bounds_bracket_the_exact_throughput(rho, N):
    L = np.array(rho).reshape(-1, 1)
    LB, UB, TH = pfqn_harel_bounds(np.array(rho), N)
    X = float(np.atleast_1d(pfqn_mva(L, np.array([N]))[0]).ravel()[0])
    assert LB <= X * (1 + 1e-12)
    for n in range(2, min(N, 7) + 1):
        assert X <= UB[n - 1] * (1 + 1e-12)


def test_harel_single_entry_points_agree_with_the_bundle():
    rho = np.array([0.9, 0.4, 0.2])
    LB, UB, _ = pfqn_harel_bounds(rho, 10)
    assert pfqn_harel_lb(rho, 10) == pytest.approx(LB, rel=1e-15)
    assert pfqn_harel_ub(rho, 10, 5) == pytest.approx(UB[4], rel=1e-15)


def test_harel_refusals():
    rho = np.array([0.5, 0.5])
    with pytest.raises(ValueError):
        pfqn_harel_lb(rho, 4, 1.0)          # think time is refused, not folded in
    with pytest.raises(ValueError):
        pfqn_harel_ub(rho, 10, 8)           # the G(n) ceiling of the reference
    with pytest.raises(ValueError):
        pfqn_harel_bounds(np.array([0.5, -0.5]), 4)


@pytest.mark.parametrize("L", [[0.9, 0.4, 0.2], [1.2, 0.3, 0.1], [2.0, 1.0, 0.5, 0.25]])
@pytest.mark.parametrize("N", [5, 20, 60])
def test_xia_single_server_is_the_dominant_pole_asymptotics(L, N):
    L = np.array(L)
    Lmax = L.max()
    closed_form = ((N + L.size - 1) * math.log(Lmax) -
                   float(np.sum(np.log(Lmax - L[L < Lmax]))))
    assert pfqn_xia(L, N, np.ones(L.size)) == pytest.approx(closed_form, rel=1e-12)


@pytest.mark.parametrize("N,tol", [(10, 1e-4), (30, 1e-8), (80, 1e-10)])
def test_xia_multiserver_converges_to_the_exact_constant(N, tol):
    L = np.array([0.9, 0.4])
    s = np.array([2.0, 3.0])
    mu = np.zeros((2, N))
    for i in range(2):
        for n in range(1, N + 1):
            mu[i, n - 1] = min(n, s[i])
    exact = pfqn_gld(L.reshape(-1, 1), np.array([N]), mu)
    lG = exact.lG if hasattr(exact, 'lG') else exact[1]
    assert pfqn_xia(L, N, s) == pytest.approx(lG, abs=tol * max(1.0, abs(lG)))
