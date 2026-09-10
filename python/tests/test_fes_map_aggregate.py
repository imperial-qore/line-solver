"""
Contract of the MAP flow-equivalent server of Casale, Mi, Cherkasova and Smirni,
IEEE Trans. Soft. Eng. 37(5), 2011, Section 5.2.

The mean inter-departure time of the aggregated subnetwork is the reciprocal of its
throughput, so for exponential service the recursion must reproduce exact MVA at
every population: the MAP flow-equivalent server degrades to the classic Norton
one, which is exact for product-form. Burstiness must survive the aggregation when
the subnetwork carries it. Expected values come from the MATLAB reference
(matlab/src/api/fes/fes_map_aggregate.m).
"""

import numpy as np
import pytest

from line_solver.api.fes import (fes_map_aggregate, fes_map_deaggregate, fes_map_grid,
                                 fes_map_interdeparture, fes_map_interp, fes_map_levels,
                                 fes_map_moments, fes_map_solve)
from line_solver.api.mam import map_exponential, map_idc, map_moment, map_scv, map2_fit_idc
from line_solver.api.pfqn import pfqn_mva

TOL = 1e-9


def _mva_throughput(demands, n, servers=None):
    L = np.array(demands).reshape(-1, 1)
    mi = np.ones(len(demands)) if servers is None else np.array(servers, dtype=float)
    out = pfqn_mva(L, np.array([n], dtype=float), np.array([0.0]), mi)
    return float(np.ravel(out[0])[0])


def test_interdeparture_mean_is_inverse_throughput():
    mu1, mu2 = 1.5, 1.0
    station = map_exponential(1 / mu1)
    fes = map_exponential(1 / mu2)
    for n in (1, 2, 5, 10):
        T0, T1 = fes_map_interdeparture(station, fes, n)
        e1, e2, e3, e11, idc = fes_map_moments(T0, T1)
        x = _mva_throughput([1 / mu1, 1 / mu2], n)
        assert abs(e1 - 1 / x) <= TOL / x


def test_euler_quadrature_converges_with_the_step():
    station = map_exponential(0.8)
    fes = map_exponential(1.3)
    T0, T1 = fes_map_interdeparture(station, fes, 6)
    exact = fes_map_moments(T0, T1, 'ssolve')[0]
    coarse = abs(fes_map_moments(T0, T1, 'euler', 0.1)[0] - exact)
    fine = abs(fes_map_moments(T0, T1, 'euler', 0.01)[0] - exact)
    assert coarse / exact < 5e-2
    assert fine < 0.2 * coarse


def test_exponential_subnetwork_reproduces_exact_mva():
    rates = [2.0, 1.5, 1.0]
    maps = [map_exponential(1 / r) for r in rates]
    nmax = 8
    fes, info = fes_map_aggregate(maps, [1, 1, 1], nmax)
    for k in range(1, nmax + 1):
        x = _mva_throughput([1 / r for r in rates], k)
        assert abs(info['throughput'][k - 1] - x) <= TOL * x


def test_delay_and_multiserver_subnetwork():
    maps = [map_exponential(0.5), map_exponential(1.0)]
    nmax = 6
    fes, info = fes_map_aggregate(maps, [np.inf, 2], nmax)
    for k in range(1, nmax + 1):
        p = np.ones(k + 1)
        for j in range(1, k + 1):
            p[j] = p[j - 1] * ((k - j + 1) * 2.0) / min(j, 2)
        p = p / p.sum()
        x = float(np.sum(p * np.minimum(np.arange(k + 1), 2)))
        assert abs(info['throughput'][k - 1] - x) <= 1e-8 * x


def test_burstiness_survives_aggregation():
    bursty = (np.array([[-1.9, 0.0], [0.0, -0.1]]),
              np.array([[1.71, 0.19], [0.01, 0.09]]))
    maps = [bursty, map_exponential(1.2)]
    fes, info = fes_map_aggregate(maps, [1, 1], 5)
    for k in range(5):
        assert info['status'][k] == 0
        assert abs(map_idc(fes[k][0], fes[k][1]) - info['moments'][3, k]) <= 1e-6 * info['moments'][3, k]
        assert abs(map_moment(fes[k][0], fes[k][1], 1) - info['moments'][0, k]) <= 1e-9
        assert map_scv(fes[k][0], fes[k][1]) > 1


def test_fit_falls_back_to_exponential_when_not_bursty():
    fit, status = map2_fit_idc(1.0, 1.5, 4.0, 0.7)
    assert status == 1
    assert abs(map_moment(fit[0], fit[1], 1) - 1.0) <= TOL


def test_grid_skips_levels_only_above_twenty():
    assert fes_map_grid(20).size == 20
    grid = fes_map_grid(100)
    assert grid.size < 100
    assert grid[0] == 1 and grid[-1] == 100


def test_interp_is_shape_preserving():
    x = np.array([1.0, 2.0, 5.0, 9.0])
    y = np.array([1.0, 2.0, 2.5, 2.6])
    xq = np.linspace(1, 9, 41)
    yq = fes_map_interp(x, y, xq).ravel()
    assert np.all(np.diff(yq) >= -1e-12)
    assert yq.max() <= y.max() + 1e-12


def test_reduced_model_solve_matches_exact_mva():
    # the closed model of a delay and an aggregated exponential subnetwork is
    # product-form, so the reduced solve must return exact MVA
    rates = [2.0, 1.5, 1.0]
    z = 0.8
    maps = [map_exponential(1 / r) for r in rates]
    demands = [1 / r for r in rates]
    for n in (1, 3, 6, 10):
        fes, _ = fes_map_aggregate(maps, [1, 1, 1], n)
        x, r, q, pk = fes_map_solve(fes, map_exponential(z), n)
        L = np.array(demands).reshape(-1, 1)
        out = pfqn_mva(L, np.array([float(n)]), np.array([z]), np.ones(3))
        assert abs(x - float(np.ravel(out[0])[0])) <= TOL * float(np.ravel(out[0])[0])
        assert abs(r - float(np.sum(out[4]))) <= 1e-8 * float(np.sum(out[4]))
        assert abs(pk.sum() - 1.0) <= 1e-12


def test_deaggregation_matches_exact_mva_per_station():
    rates = [2.0, 1.5, 1.0]
    z = 0.8
    n = 8
    maps = [map_exponential(1 / r) for r in rates]
    demands = np.array([1 / r for r in rates])
    fes, _ = fes_map_aggregate(maps, [1, 1, 1], n)
    _, _, _, pk = fes_map_solve(fes, map_exponential(z), n)
    qn, un, xn, rn = fes_map_deaggregate(pk, demands, [1, 1, 1], [False, False, False])

    out = pfqn_mva(demands.reshape(-1, 1), np.array([float(n)]), np.array([z]), np.ones(3))
    xref = float(np.ravel(out[0])[0])
    qref = np.ravel(np.asarray(out[2]))
    uref = np.ravel(np.asarray(out[3]))
    for i in range(3):
        assert abs(qn[i] - qref[i]) <= 1e-8 * qref[i]
        assert abs(un[i] - uref[i]) <= 1e-8 * uref[i]
        assert abs(xn[i] - xref) <= 1e-8 * xref
