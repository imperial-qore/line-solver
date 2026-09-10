"""
Regression tests: the ETAQA measures of BMAP/MAP/1 and MAP/BMAP/1.

Both entry points used to answer with SUBSTITUTES rather than MAMSolver's
ETAQA. `mg1_g_etaqa` ran a plain functional iteration where the reference runs
Bini-Meini cyclic reduction, `mg1_qlen_etaqa` reported
sum(pi1) + 2 sum(pi*) -- which assumes every level at or above 2 IS level 2 --
and the whole GI/M/1 side (R by naive iteration, pi matrix-geometric, moments
by a geometric guess over the tail) never touched GIM1_R_ETAQA at all. On the
four-phase cases below that read 1.0887 against MATLAB's 2.4961 and 0.3840
against MATLAB's -3.4342.

The oracles here are of two kinds and both matter:

 1. CLOSED FORMS. A unit-batch BMAP with Poisson arrivals and exponential
    service is an M/M/1, whose level is geometric, so E[N] and its higher
    moments are exact numbers. A genuine batch of size 1 or 2 is an M[X]/M/1,
    whose mean the unit-batch case cannot distinguish and a MAP surrogate of
    the same rate would miss entirely.
 2. MATLAB VALUES on the four-phase cases, including the NEGATIVE mean on the
    GI/M/1 side. That negative value is not a defect of this port: MATLAB's
    GIM1_qlen_ETAQA initializes its accumulator with the SCALAR A(3) where the
    third BLOCK is meant, which corrupts the last column of the moment system
    for more than one phase. It is reproduced deliberately, and pi and R -- the
    outputs the defect does not touch -- are pinned alongside it.
"""
import numpy as np
import pytest

from line_solver.api.mam import solver_mam_bmap_map_1, solver_mam_map_bmap_1

ARR_D0 = np.array([[-1.4, 0.2], [0.3, -0.8]])
ARR_D1 = np.array([[1.2, 0.0], [0.0, 0.5]])
SVC_D0 = np.array([[-3.0, 0.5], [0.1, -2.0]])
SVC_D1 = np.array([[2.0, 0.5], [1.4, 0.5]])


def test_bmap_map_1_unit_batch_is_mm1():
    D = [np.array([[-0.6]]), np.array([[0.6]])]
    r = solver_mam_bmap_map_1(D, np.array([[-1.0]]), np.array([[1.0]]))
    rho = 0.6
    assert r.mean_queue_length == pytest.approx(rho / (1 - rho), rel=1e-10)
    assert r.utilization == pytest.approx(rho, rel=1e-12)
    assert r.throughput == pytest.approx(0.6, rel=1e-12)
    assert r.mean_response_time == pytest.approx(2.5, rel=1e-10)
    # A rank-one A0 sends MG1_EG home with G = 1 before cyclic reduction runs.
    assert float(np.asarray(r.G).ravel()[0]) == pytest.approx(1.0, rel=1e-12)
    # The aggregates are (1-rho, rho(1-rho), rho^2).
    pi = np.asarray(r.pi).ravel()
    assert pi[0] == pytest.approx(0.4, rel=1e-10)
    assert pi[1] == pytest.approx(0.24, rel=1e-10)
    assert pi[2] == pytest.approx(0.36, rel=1e-10)


def test_bmap_map_1_batch_matches_mx_m_1():
    # Batch sizes 1 and 2 with equal probability, batch rate 0.3, lambda 0.45.
    D = [np.array([[-0.3]]), np.array([[0.15]]), np.array([[0.15]])]
    r = solver_mam_bmap_map_1(D, np.array([[-1.0]]), np.array([[1.0]]))
    rho, ex, exx1 = 0.45, 1.5, 1.0
    exact = rho / (1 - rho) + rho * (exx1 / ex) / (2 * (1 - rho))
    assert r.mean_queue_length == pytest.approx(exact, rel=1e-10)
    assert r.throughput == pytest.approx(0.45, rel=1e-12)


def test_bmap_map_1_matches_matlab_on_a_map_input():
    D = [ARR_D0, 0.5 * ARR_D1, 0.5 * ARR_D1]
    r = solver_mam_bmap_map_1(D, SVC_D0, SVC_D1)
    assert r.mean_queue_length == pytest.approx(2.496123595502, rel=1e-10)
    assert r.utilization == pytest.approx(0.610619469027, rel=1e-10)
    assert r.mean_response_time == pytest.approx(1.808785214132, rel=1e-10)
    assert r.throughput == pytest.approx(1.38, rel=1e-12)

    pi = np.asarray(r.pi).ravel()
    assert np.sum(pi) == pytest.approx(1.0, rel=1e-12)
    assert pi[0] == pytest.approx(0.118941068513, rel=1e-9)
    assert pi[4] == pytest.approx(0.052543965796, rel=1e-9)
    assert pi[8] == pytest.approx(0.188514965691, rel=1e-9)

    # G is stochastic: a positive recurrent chain leaves a level downwards with
    # probability one.
    G = np.asarray(r.G)
    assert np.max(np.abs(np.sum(G, axis=1) - 1.0)) < 1e-9
    assert G[0, 0] == pytest.approx(0.686039381261, rel=1e-9)
    assert G[3, 3] == pytest.approx(0.217998974669, rel=1e-9)


def test_bmap_map_1_higher_moments_match_matlab():
    from line_solver.lib.thirdparty.smc import (mg1_g_etaqa, mg1_pi_etaqa,
                                                mg1_qlen_etaqa)
    ma = ms = 2
    m = ma * ms
    K = 2
    D = [ARR_D0, 0.5 * ARR_D1, 0.5 * ARR_D1]
    A = np.zeros((m, m * (K + 2)))
    A[:, 0:m] = np.kron(np.eye(ma), SVC_D1)
    A[:, m:2*m] = np.kron(D[0], np.eye(ms)) + np.kron(np.eye(ma), SVC_D0)
    for k in range(1, K + 1):
        A[:, (k+1)*m:(k+2)*m] = np.kron(D[k], np.eye(ms))
    B = np.zeros((m, m * (K + 1)))
    B[:, 0:m] = np.kron(D[0], np.eye(ms)) + np.kron(np.eye(ma), SVC_D0 + SVC_D1)
    for k in range(1, K + 1):
        B[:, k*m:(k+1)*m] = np.kron(D[k], np.eye(ms))

    G = mg1_g_etaqa(A)
    pi = mg1_pi_etaqa(B, A, G)
    moments = [mg1_qlen_etaqa(B, A, pi, n) for n in (1, 2, 3)]
    assert moments[0] == pytest.approx(2.496123595502, rel=1e-10)
    assert moments[1] == pytest.approx(17.835504350850, rel=1e-9)
    assert moments[2] == pytest.approx(191.796380027147, rel=1e-9)


def test_map_bmap_1_unit_batch_is_mm1():
    D = [np.array([[-1.0]]), np.array([[1.0]])]
    r = solver_mam_map_bmap_1(np.array([[-0.6]]), np.array([[0.6]]), D)
    assert r.mean_queue_length == pytest.approx(1.5, rel=1e-10)
    assert r.mean_response_time == pytest.approx(2.5, rel=1e-10)
    assert float(np.asarray(r.R).ravel()[0]) == pytest.approx(0.6, rel=1e-10)
    pi = np.asarray(r.pi).ravel()
    assert pi[0] == pytest.approx(0.4, rel=1e-10)
    assert pi[2] == pytest.approx(0.36, rel=1e-10)


def test_map_bmap_1_single_phase_batch_service():
    D = [np.array([[-1.0]]), np.array([[0.5]]), np.array([[0.5]])]
    r = solver_mam_map_bmap_1(np.array([[-0.9]]), np.array([[0.9]]), D)
    assert r.mean_queue_length == pytest.approx(2.0611000442, rel=1e-9)
    assert r.utilization == pytest.approx(0.6, rel=1e-12)
    assert float(np.asarray(r.R).ravel()[0]) == pytest.approx(0.6733200531, rel=1e-9)


def test_map_bmap_1_matches_matlab_including_the_negative_mean():
    D = [SVC_D0, 0.5 * SVC_D1, 0.5 * SVC_D1]
    r = solver_mam_map_bmap_1(ARR_D0, ARR_D1, D)
    # Negative, and negative in MATLAB too, to twelve digits: see the module
    # docstring and lib/thirdparty/smc.gim1_qlen_etaqa.
    assert r.mean_queue_length == pytest.approx(-3.434197005219, rel=1e-9)
    assert r.utilization == pytest.approx(0.271386430678, rel=1e-10)
    assert r.throughput == pytest.approx(0.92, rel=1e-12)

    pi = np.asarray(r.pi).ravel()
    assert np.sum(pi) == pytest.approx(1.0, rel=1e-12)
    assert pi[0] == pytest.approx(0.227094812985, rel=1e-9)
    assert pi[4] == pytest.approx(0.082609049829, rel=1e-9)
    assert pi[8] == pytest.approx(0.050296137186, rel=1e-9)

    R = np.asarray(r.R)
    assert R[0, 0] == pytest.approx(0.319301430437, rel=1e-9)
    assert R[3, 3] == pytest.approx(0.187111846865, rel=1e-9)
