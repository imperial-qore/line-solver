"""
Known-answer tests for the formulas ported from A. Thomasian, "Analysis of
Fork/Join and Related Queueing Systems", ACM Computing Surveys 47(2),
Article 17, 2014.

Every expected value is either the survey's own worked number, an exact identity
the formula must reproduce, or the value the MATLAB fj_* twin returns
(matlab/src/api/fj/), which is the reference implementation.
"""

import math

import numpy as np

from line_solver.api.fjnative import (
    fj_amva, fj_char_max_blom, fj_char_max_discrete, fj_cox_fit, fj_dag_makespan,
    fj_delay_opt, fj_dispersion, fj_harmonic, fj_ism_green, fj_lst_max_het, fj_qgb,
    fj_respt_bulk, fj_respt_closed, fj_respt_nosplit, fj_serialization,
    fj_tsm_capacity, fj_xmax_coxian, fj_xmax_het, fj_xmax_hz, fj_xmax_hz_het,
    fj_xmax_moments_het,
)


def test_qgb_reduces_to_the_ordinary_geometric_bound():
    D = [1.0, 2.0, 3.0]
    Q, y = fj_qgb(D, [1, 1, 1], 7, 1.5)
    denom = 1.5 + 6.0 + 3.0 * 7.0
    for i in range(3):
        yi = D[i] * 7.0 / denom
        assert abs(Q[i] - (yi / (1 - yi) - yi ** 8 / (1 - yi))) < 1e-12


def test_amva_with_unit_fork_degrees_is_exact_mva():
    D = np.array([1.0, 2.0, 3.0])
    R, Q, X, U = fj_amva(D, [1, 1, 1], 7, 1.5)
    Qref = np.zeros(3)
    Xref = 0.0
    for m in range(1, 8):
        Rref = D * (1 + Qref)
        Xref = m / (1.5 + Rref.sum())
        Qref = Xref * Rref
    assert abs(X - Xref) < 1e-12
    assert np.allclose(Q, Qref, atol=1e-12)


def test_respt_closed_is_tight_at_two_branches():
    R, exact = fj_respt_closed(2, 0.4, 5)
    assert abs(R - 0.4 * (1.5 + 4)) < 1e-9
    assert exact


def test_xmax_het_matches_the_textbook_two_variable_answer():
    assert abs(fj_xmax_het([1.0, 2.0]) - (1 + 0.5 - 1 / 3)) < 1e-9
    assert abs(fj_xmax_het([2.0] * 4) - fj_harmonic(4) / 2) < 1e-9


def test_moment_recurrence_agrees_with_inclusion_exclusion():
    lam = [1.0, 2.0, 3.0, 5.0]
    m = fj_xmax_moments_het(lam, 3)
    for n in (1, 2, 3):
        assert abs(m[n - 1] - fj_xmax_het(lam, n)) < 1e-10


def test_transform_is_one_at_the_origin_and_its_slope_is_the_mean():
    lam = [1.0, 2.0, 3.0, 5.0]
    assert abs(fj_lst_max_het(lam, 0.0) - 1.0) < 1e-12
    h = 1e-5
    assert abs(-(fj_lst_max_het(lam, h) - 1.0) / h - fj_xmax_het(lam)) < 1e-4


def test_harrison_zertal_is_exact_for_the_exponential():
    m1 = 0.7
    x, _ = fj_xmax_hz(m1, 2 * m1 ** 2, 6)
    assert abs(x - fj_harmonic(6) * m1) < 1e-9


def test_harrison_zertal_recurrence_is_exact_for_iid_exponentials():
    K = 4
    lam = 1.3
    cdf = [lambda t, lam=lam: np.where(t > 0, 1 - np.exp(-lam * np.maximum(t, 0)), 0.0)] * K
    got = fj_xmax_hz_het([1 / lam] * K, [2 / lam ** 2] * K, cdf)
    assert abs(got - fj_harmonic(K) / lam) < 1e-6


def test_characteristic_maximum_bounds_the_lattice_maximum():
    MK, mK, exact = fj_char_max_discrete(8, 'geometric', 0.6)
    assert MK >= exact - 1e-9
    MK, mK, exact = fj_char_max_discrete(8, 'poisson', 4.0)
    assert MK >= exact - 1e-9


def test_blom_position_sits_inside_the_kruskal_weiss_bracket():
    mK, lo, hi = fj_char_max_blom(20)
    assert lo < mK < hi


def test_coxian_fit_reproduces_its_targets_and_degenerates_to_the_exponential():
    mu1, mu2, q, _, _ = fj_cox_fit(2.0, 1.5)
    x, m1, c2 = fj_xmax_coxian(1, mu1, mu2, q)
    assert abs(m1 - 2.0) < 1e-12
    assert abs(c2 - 1.5) < 1e-12
    assert abs(x - 2.0) < 1e-12
    x5, _, _ = fj_xmax_coxian(5, 1.0, 1.0, 0.0)
    assert abs(x5 - fj_harmonic(5)) < 1e-10


def test_dispersion_of_two_exponential_branches():
    mu = 1.7
    Edisp, Emax, Emin = fj_dispersion([1, 1], [mu, mu])
    assert abs(Emax - 1.5 / mu) < 1e-7
    assert abs(Emin - 0.5 / mu) < 1e-7
    assert abs(Edisp - 1.0 / mu) < 1e-7


def test_delaying_never_increases_the_dispersion():
    shape = [1, 3, 2]
    rate = [1.0, 2.0, 0.8]
    d0 = fj_dispersion(shape, rate)[0]
    d, Edisp, _ = fj_delay_opt(shape, rate)
    assert Edisp <= d0 + 1e-9
    assert abs(float(np.min(d))) < 1e-12


def test_no_splitting_collapses_to_mm1_at_one_task():
    R, _ = fj_respt_nosplit(1, 0.5, 1.4)
    assert abs(R - 1 / 0.9) < 1e-9


def test_bulk_arrivals_collapse_to_mm1_at_unit_batch_and_one_server():
    Rreq, Rtask, Q, _ = fj_respt_bulk(1, 0.5, 1.4, 1)
    rho = 0.5 / 1.4
    assert abs(Q - rho / (1 - rho)) < 1e-7
    assert abs(Rreq - 1 / 0.9) < 1e-7
    assert abs(Rtask - Rreq) < 1e-7


def test_independent_server_model_reproduces_the_surveys_utilization():
    W, R, out = fj_ism_green(1.0, 1.4, 4, [0.1, 0.2, 0.3, 0.4])
    # The survey quotes rho = 0.9286 for this instance
    assert abs(out['rho'] - 0.9285714285714286) < 1e-9
    assert W > 0
    assert 0 < out['pq'] < 1
    assert 0 < out['pd'] <= 1


def test_team_service_capacity_matches_the_surveys_worked_example():
    f = [0.25] * 4
    r = [1, 2, 3, 4]
    x = [1.0] * 4
    Lmax, Llp, Lfcfs, _ = fj_tsm_capacity(4, f, r, x)
    # Lambda_max = s / sum f r x = 4 / 2.5 = 1.6, attained because every state
    # carrying probability is full capacity
    assert abs(Lmax - 1.6) < 1e-12
    assert abs(Llp - 1.6) < 1e-7
    eps = 1e-6
    _, Llp2, _, _ = fj_tsm_capacity(4, [eps, 0.5 - eps, 0.5 - eps, eps], r, x)
    # The skewed frequency vector of the survey reaches only 4/3
    assert abs(Llp2 - 4 / 3) < 1e-4


def test_team_service_fcfs_capacity_of_the_two_server_case():
    _, _, Lfcfs, _ = fj_tsm_capacity(2, [0.5, 0.5], [1, 2], [0.5, 1 / 3])
    mu1, mu2, f1, f2 = 2.0, 3.0, 0.5, 0.5
    want = 2 * mu1 * mu2 / (f1 ** 2 * mu2 + 2 * f2 ** 2 * mu1 + 2 * f1 * f2 * (mu1 + mu2))
    assert abs(Lfcfs - want) < 1e-12


def test_serialization_blocking_probability():
    P, delay, Rtot = fj_serialization([0.2, 0.3], 1.0, 5)
    assert abs(P[0] - (1 - (1 - 0.2 / 1.5) ** 4)) < 1e-12


def test_dag_makespan_reproduces_both_worked_examples():
    pred = np.zeros((2, 2))
    rate = np.array([[1 / 10, 1 / 15], [1 / 20, 1 / 30]])
    C, I, Cend, E = fj_dag_makespan(pred, rate)
    # The survey's coupled two-task example gives 26.77 by hand, 80/3 exactly
    assert abs(C - 80 / 3) < 1e-9
    assert abs(Cend[0] - 40 / 3) < 1e-9
    assert abs(Cend[1] - 70 / 3) < 1e-9
    C2, _, _, _ = fj_dag_makespan(pred, np.array([1.0, 0.1]))
    # The state-truncation example gives 10.1
    assert abs(C2 - (1 / 1.1 + (1 / 1.1) * 10 + (0.1 / 1.1) * 1)) < 1e-9


def test_dag_makespan_of_a_chain_is_the_sum_of_the_means():
    pred = np.zeros((3, 3))
    pred[0, 1] = 1
    pred[1, 2] = 1
    C, I, Cend, E = fj_dag_makespan(pred, np.array([1.0, 2.0, 4.0]))
    assert abs(C - (1 + 0.5 + 0.25)) < 1e-12
    assert abs(I[2] - 1.5) < 1e-12
