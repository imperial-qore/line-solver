"""
Birman-Kogan asymptotics (Stochastic Models 8(3):543-563, 1992).

The multiprogramming model of the paper's Table 3 is fully specified in print --
J = 2 device groups of 10 and 50 stations, service times .015 and .030, five jobs
per chain, branching .50/.50 and .6/.4 alternating, cpu rate mu_k = M*mu0 with
M = 50 -- and the paper tabulates four algorithms on it, so the expected values
below pin the algorithms rather than this implementation.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_bk, pfqn_bkue, pfqn_bklc, pfqn_ca, pfqn_kt

TH = [1 / 0.015, 1 / 0.030]
MG = [10, 50]
MPAR = 50


def branch(k, j):
    return 0.5 if k % 2 == 0 else (0.6 if j == 0 else 0.4)


def multiprogramming(K, mu0):
    """The Table 3 network: K dedicated cpus followed by the two device groups."""
    rows = []
    for k in range(K):
        r = np.zeros(K)
        r[k] = 1.0 / (MPAR * mu0)
        rows.append(r)
    for j in range(2):
        row = np.array([branch(k, j) / (TH[j] * MG[j]) for k in range(K)])
        for _ in range(MG[j]):
            rows.append(row.copy())
    return np.array(rows)


@pytest.mark.parametrize("K,want", [(2, (90.9, 95.0)), (3, (83.2, 86.0)),
                                    (4, (75.6, 76.9)), (5, (69.7, 70.2))])
def test_saddle_point_reproduces_table3_column(K, want):
    # U_k = x_k^0 / mu_0k by Corollary 1. x^0 does not depend on mu_0k, so the
    # mu_0 = 4 row fixes every other row of the column by a pure rescaling.
    L = multiprogramming(K, 4.0)
    N = np.full(K, 5.0)
    _, _, X, _, A, B = pfqn_bk(L, N, np.zeros(K))
    for k in range(2):
        assert min(X[k] * L[k, k], 1.0) * 100 == pytest.approx(want[k], abs=0.05)
    assert len(B) == 0  # no chain saturates at mu_0 = 4


def test_algorithm1_pins_every_saturated_chain():
    L = multiprogramming(2, 2.0)
    _, _, X, _, A, B = pfqn_bk(L, np.full(2, 5.0), np.zeros(2))
    assert len(A) == 0 and len(B) == 2
    for k in range(2):
        assert X[k] * L[k, k] >= 1.0


@pytest.mark.parametrize("K,mu0,mva,ue", [
    (2, 4.0, (70.2, 72.3), (69.8, 72.0)),
    (3, 4.0, (67.2, 69.1), (66.8, 68.7)),
    (2, 2.0, (93.1, 94.0), (93.4, 94.4)),
    (5, 1.5, (95.9, 96.4), (96.3, 96.8)),
])
def test_loadconcealment_reproduces_both_table3_columns(K, mu0, mva, ue):
    L = multiprogramming(K, mu0)
    N = np.full(K, 5.0)
    Z = np.zeros(K)
    Xm = pfqn_bklc(L, N, Z, 'mva', 1e-12, 2000)[0]
    Xu = pfqn_bklc(L, N, Z, 'ue', 1e-12, 2000)[0]
    for k in range(2):
        assert Xm[k] * L[k, k] * 100 == pytest.approx(mva[k], abs=0.2)
        assert Xu[k] * L[k, k] * 100 == pytest.approx(ue[k], abs=0.3)


def test_saddle_point_equals_kt_without_dedicated_stations():
    # With a think time in every chain there are no dedicated stations to keep out
    # of the exponent, so Proposition 1 IS the multidimensional saddle point that
    # pfqn_kt already computes.
    L = np.array([[1.0, 0.6], [0.8, 1.2], [0.5, 0.9]])
    N = np.array([10.0, 8.0])
    Z = np.array([2.0, 1.0])
    assert pfqn_bk(L, N, Z)[1] == pytest.approx(pfqn_kt(L, N, Z)[1], abs=1e-9)


def test_uniform_expansion_is_exact_on_one_slow_station_against_a_group():
    # The regime the uniform expansion is built for: a single dominant pole
    # against M >> 1 identical stations.
    L = np.vstack([[0.9], np.full((12, 1), 0.1)])
    N = np.array([20.0])
    Z = np.array([0.0])
    exact = pfqn_ca(L, N, Z)[1]
    assert pfqn_bkue(L.flatten(), 20.0, 0.0)[1] == pytest.approx(exact, abs=1e-4)
    assert pfqn_bk(L, N, Z)[1] == pytest.approx(exact, abs=1e-4)


def test_uniform_expansion_degenerates_without_a_group():
    L = np.array([[0.5], [0.4], [0.3]])
    N = np.array([60.0])
    Z = np.array([5.0])
    assert pfqn_bkue(L.flatten(), 60.0, 5.0)[1] == pytest.approx(pfqn_kt(L, N, Z)[1], abs=1e-9)


def test_nc_solver_methods_are_advertised_and_dispatch():
    from line_solver.api.pfqn.nc import pfqn_nc
    L = np.array([[1.0, 0.6], [0.8, 1.2], [0.5, 0.9]])
    N = np.array([10.0, 8.0])
    Z = np.array([2.0, 1.0])
    lg_ca = pfqn_nc(L, N, Z, method='ca')[1]
    for m in ('bk', 'bkue', 'lc', 'lc.ue'):
        lg = pfqn_nc(L, N, Z, method=m)[1]
        assert np.isfinite(lg)
        assert abs(lg - lg_ca) < 1.0
