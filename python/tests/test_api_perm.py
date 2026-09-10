"""Tests for the matrix permanent algorithms in line_solver.api.perm.

The exact solvers (inclusion-exclusion with multiplicities, Ryser in both
variants, permutation enumeration) must agree with each other to machine
precision and with the hand-computed permanents of small matrices; they are the
twins of MATLAB perm.m and of jline.lib.perm.{Permanent, RyzerPermanent,
NaivePermanent}.

The approximations are checked against the exact value with the tolerance each
method can guarantee: the Bethe permanent is a lower bound within a factor
sqrt(2)^n of the permanent (Gurvits, Vontobel), the Sinkhorn heuristic is a
mean-field estimate with no error bound and is only checked for order of
magnitude, and the two samplers are unbiased so their error is checked at a
loose multiple of the sampling standard error.
"""

import itertools
import math
import time

import pytest

import numpy as np

from line_solver.api.perm import (
    AdaPartSampler,
    BethePermanent,
    HeuristicPermanent,
    HuberLawSampler,
    NaivePermanent,
    NetworkNoThink,
    NetworkThink,
    Permanent,
    RyzerPermanent,
    SaddlePointPermanent,
    compute_permanent,
    perm,
    perm_bethe,
    perm_heur,
    perm_spm,
    permanent,
    preprocessing_ds,
)
from line_solver.api.perm.exact import _permanent_ryser

M2 = np.array([[1.0, 2.0], [3.0, 4.0]])
M3 = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
PERM_M2 = 10.0
PERM_M3 = 450.0


def _random_matrix(n, rng):
    return rng.random((n, n)) + 0.05


def test_exact_small_matrices():
    # perm([[1,2],[3,4]]) = 1*4 + 2*3, perm(M3) = 450 by direct expansion
    assert abs(compute_permanent(M2) - PERM_M2) < 1e-12
    assert abs(permanent(M2) - PERM_M2) < 1e-12
    assert abs(perm(M2) - PERM_M2) < 1e-12
    assert abs(compute_permanent(M3) - PERM_M3) < 1e-12
    assert abs(compute_permanent(M3, use_multiplicities=False) - PERM_M3) < 1e-12
    assert abs(_permanent_ryser(M3) - PERM_M3) < 1e-12


def _leibniz_permanent(a):
    """Brute-force permanent from its definition, independent of the module."""
    n = a.shape[0]
    total = 0.0
    for p in itertools.permutations(range(n)):
        term = 1.0
        for i in range(n):
            term *= a[i, p[i]]
        total += term
    return total


def test_ryser_sign_pinned_to_brute_force_for_odd_n():
    # Regression: _permanent_ryser applied the (-1)^n subset sign AND a further
    # (-1)^n, so it returned -per(A) for every odd n. Pin both parities against
    # the Leibniz definition, and require a positive value on a positive matrix
    rng = np.random.default_rng(2026)
    for n in range(1, 7):
        for _ in range(3):
            a = rng.random((n, n)) + 0.05
            expected = _leibniz_permanent(a)
            assert expected > 0.0
            assert abs(_permanent_ryser(a) - expected) < 1e-9 * expected
            assert abs(compute_permanent(a, use_multiplicities=False) - expected) < 1e-9 * expected
            assert abs(compute_permanent(a) - expected) < 1e-9 * expected
            assert abs(RyzerPermanent(a, "graycode", True).value - expected) < 1e-9 * expected
            assert abs(RyzerPermanent(a, "naive", True).value - expected) < 1e-9 * expected
            assert abs(NaivePermanent(a, True).value - expected) < 1e-9 * expected
            assert abs(Permanent(a, True).value - expected) < 1e-9 * expected


def test_exact_solvers_agree():
    rng = np.random.default_rng(42)
    for n in range(2, 9):
        a = _random_matrix(n, rng)
        reference = _permanent_ryser(a)
        assert abs(Permanent(a, True).value - reference) < 1e-9 * abs(reference)
        assert abs(RyzerPermanent(a, "graycode", True).value - reference) < 1e-9 * abs(reference)
        assert abs(RyzerPermanent(a, "naive", True).value - reference) < 1e-9 * abs(reference)
        if n <= 7:
            assert abs(NaivePermanent(a, True).value - reference) < 1e-9 * abs(reference)


def test_exact_repeated_columns_and_rows():
    # The multiplicity path collapses repeated columns and repeated rows
    repeated_cols = np.array([[1.0, 1.0, 2.0], [2.0, 2.0, 3.0], [3.0, 3.0, 4.0]])
    repeated_rows = np.array([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    for a in (repeated_cols, repeated_rows):
        reference = NaivePermanent(a, True).value
        assert abs(Permanent(a, True).value - reference) < 1e-9 * abs(reference)
        assert abs(RyzerPermanent(a, "graycode", True).value - reference) < 1e-9 * abs(reference)


def test_solver_instrumentation():
    solver = Permanent(M3)
    assert solver.value == 0.0
    solver.solve()
    assert abs(solver.value - PERM_M3) < 1e-12
    assert solver.time_ms >= 0.0
    assert solver.memory_bytes >= 0
    result = solver.get_result()
    assert result.n == 3
    assert abs(result.value - PERM_M3) < 1e-12


def test_bethe_is_a_lower_bound_within_the_gurvits_factor():
    # The message update excludes the diagonal entry rather than the target
    # entry, as in the JAR twin, so the sqrt(2)^n factor of Gurvits is only
    # asserted from n = 3 on; see _kb/03-api-layer.md
    rng = np.random.default_rng(7)
    for n in range(2, 8):
        a = _random_matrix(n, rng)
        reference = _permanent_ryser(a)
        estimate = BethePermanent(a, solve=True).value
        assert estimate > 0.0
        assert estimate <= reference * (1.0 + 1e-9)
        if n >= 3:
            assert estimate * np.sqrt(2.0) ** n >= reference
        assert abs(perm_bethe(a) - estimate) < 1e-9 * estimate


def test_heuristic_matches_the_matlab_formula():
    # perm_heur.m on a doubly stochastic matrix reduces to the van der Waerden
    # mean field n!/n^n, which the Gurvits capacity term reproduces exactly
    n = 4
    a = np.full((n, n), 1.0 / n)
    expected = math.factorial(n) / float(n) ** n
    assert abs(perm_heur(a) - expected) < 1e-9 * expected

    rng = np.random.default_rng(11)
    for n in range(2, 8):
        a = _random_matrix(n, rng)
        reference = _permanent_ryser(a)
        estimate = HeuristicPermanent(a, solve=True).value
        assert 0.1 * reference < estimate < 10.0 * reference


def test_heuristic_rejects_negative_entries():
    negative = np.array([[-1.0, 2.0], [3.0, 4.0]])
    try:
        HeuristicPermanent(negative, solve=True)
        raise AssertionError("negative entries must be rejected")
    except ValueError:
        pass


def test_huber_law_sampler_accuracy():
    # 4000 draws with acceptance p give a standard error sqrt((1-p)/(p*4000))
    # on the ratio; 6 standard errors plus a 2 percent floor is the tolerance
    rng = np.random.default_rng(3)
    draws = 4000
    for n in range(2, 9):
        a = _random_matrix(n, rng)
        reference = _permanent_ryser(a)
        sampler = HuberLawSampler(a, mode="sample", number_of_samples=draws,
                                  solve=True, seed=n)
        accept = sum(sampler.sample_accepted) / float(len(sampler.sample_accepted))
        tolerance = 0.02 + 6.0 * np.sqrt(max(1.0 - accept, 1e-6) / (accept * draws))
        assert abs(sampler.value - reference) <= tolerance * reference
        assert len(sampler.perm_step) == draws


def test_adapart_sampler_accuracy():
    # Same tolerance rule as the Huber-Law sampler, at the classic-mode budget
    rng = np.random.default_rng(5)
    accepted = 200
    for n in range(2, 7):
        a = _random_matrix(n, rng)
        reference = _permanent_ryser(a)
        sampler = AdaPartSampler(a, maximum_accepted_samples=accepted, solve=True, seed=n)
        total = len(sampler.sample_accepted)
        rate = accepted / float(total)
        tolerance = 0.02 + 6.0 * np.sqrt(max(1.0 - rate, 1e-6) / (rate * total))
        assert abs(sampler.value - reference) <= tolerance * reference
        assert len(sampler.perm_step) == total


def test_samplers_are_reproducible():
    a = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
    first = HuberLawSampler(a, mode="sample", number_of_samples=200, solve=True, seed=17).value
    second = HuberLawSampler(a, mode="sample", number_of_samples=200, solve=True, seed=17).value
    assert first == second
    first = AdaPartSampler(a, maximum_accepted_samples=20, solve=True, seed=17).value
    second = AdaPartSampler(a, maximum_accepted_samples=20, solve=True, seed=17).value
    assert first == second


def test_preprocessing_ds_is_doubly_stochastic_and_rescales():
    rng = np.random.default_rng(13)
    a = _random_matrix(5, rng)
    scaled, factor = preprocessing_ds(a)
    assert np.allclose(scaled.sum(axis=0), 1.0, atol=1e-3)
    assert np.allclose(scaled.sum(axis=1), 1.0, atol=1e-3)
    reference = _permanent_ryser(a)
    assert abs(_permanent_ryser(scaled) / factor - reference) < 1e-6 * reference


def test_network_marginal_matches_the_direct_permanent():
    # Two queues, one class, three jobs: the marginal of a per-queue state is
    # the permanent of the replicated demand matrix over the state factorials
    demands = np.array([[0.4, 0.6], [0.9, 0.2]])
    network = NetworkNoThink(2, 1, [3], demands, progress=False)
    state = [2, 1]
    result = network.marginal(NaivePermanent(np.eye(1)), state)
    replicated = np.array([
        [0.4, 0.4, 0.6],
        [0.4, 0.4, 0.6],
        [0.9, 0.9, 0.2],
    ])
    expected = _permanent_ryser(replicated) / (2.0 * 1.0)
    assert abs(result.probability - expected) < 1e-9 * expected

    marginals = network.generate_marginal(RyzerPermanent(np.eye(1)))
    assert len(marginals) == 4
    assert all(m.probability >= 0.0 for m in marginals.values())


def test_network_think_joint_carries_the_think_time():
    demands = np.array([[0.4, 0.6], [0.9, 0.2]])
    think = [1.5]
    network = NetworkThink(2, 1, [3], think, demands, progress=False)
    state = np.array([[2], [1]])
    # 2 jobs at queue 1 and 1 at queue 2 leave no job thinking, so the think
    # time enters with exponent zero and the joint matches the no-think model
    no_think = NetworkNoThink(2, 1, [3], demands, progress=False)
    assert abs(network.joint(state) - no_think.joint(state)) < 1e-12

    state = np.array([[1], [1]])
    ratio = network.joint(state) / no_think.joint(state)
    assert abs(ratio - think[0]) < 1e-12


# ---------------------------------------------------------------------------
# Zeros in the demand matrix: the four approximate engines require full support
# ---------------------------------------------------------------------------

ZERO_ROW = np.array([[1.0, 2.0, 3.0], [0.0, 0.0, 0.0], [4.0, 5.0, 6.0]])
ZERO_ENTRY = np.array([[1.0, 2.0], [0.0, 3.0]])


def test_saddle_point_is_exact_on_a_single_column():
    # h == 1 leaves no direction after the homogeneity is quotiented out, so the
    # Laplace factor is empty and the expansion returns the exact n! prod_k a_k1
    rng = np.random.default_rng(19)
    for n in range(1, 7):
        col = rng.random((n, 1)) + 0.1
        expected = math.factorial(n) * float(np.prod(col))
        assert abs(perm_spm(col, [n]) - expected) <= 1e-12 * expected


def test_saddle_point_closed_form_on_the_all_ones_matrix():
    # J_n scales to itself (xi = 1), so phi = n log n, the reduced Laplacian
    # I - J/n has determinant 1/n, and the estimate is the closed form below.
    # It is the analytic anchor the MATLAB, JAR and C++ twins are held to.
    for n in range(2, 11):
        expected = (2.0 * math.pi) ** (-0.5 * (n - 1)) * n ** (n + 0.5)
        assert abs(perm_spm(np.ones((n, n))) - expected) <= 1e-9 * expected


def test_saddle_point_beats_the_capacity_it_corrects():
    # exp(log_capacity) is the Gurvits capacity, an upper bound on the permanent.
    # The Gaussian factor is what turns that e^n-scale bound into a usable
    # estimate, so it must land strictly nearer the exact value.
    rng = np.random.default_rng(23)
    for n in range(3, 9):
        a = _random_matrix(n, rng)
        reference = _permanent_ryser(a)
        solver = SaddlePointPermanent(a, solve=True)
        capacity = math.exp(solver.log_capacity)
        assert capacity >= reference * (1.0 - 1e-9)
        assert abs(solver.value - reference) < abs(capacity - reference)
        assert abs(solver.log_value - math.log(solver.value)) < 1e-12
        # measured bias at unit multiplicities, near (e/sqrt(2 pi))^n
        assert reference <= solver.value <= reference * (math.e / math.sqrt(2.0 * math.pi)) ** n * 1.1


def test_saddle_point_sharpens_as_the_column_multiplicities_grow():
    # With h fixed this is a genuine asymptotic expansion in min(m). The matrix
    # whose 2k rows all equal (a, b) with m = (k, k) has permanent (2k)! (ab)^k,
    # and the expansion returns exactly 4^k / (C(2k,k) sqrt(pi k)) times it: an
    # analytic ratio, free of the matrix, that decreases to 1 like 1 + 1/(8k).
    a, b = 0.7, 1.3
    previous = np.inf
    for k in range(1, 9):
        rows = np.tile(np.array([[a, b]]), (2 * k, 1))
        exact = math.factorial(2 * k) * (a * b) ** k
        ratio = perm_spm(rows, [k, k]) / exact
        closed = 4.0 ** k / (math.comb(2 * k, k) * math.sqrt(math.pi * k))
        assert abs(ratio - closed) <= 1e-9 * closed
        assert abs(ratio - (1.0 + 1.0 / (8.0 * k))) < 0.01 / k
        assert ratio < previous
        previous = ratio


def test_saddle_point_agrees_with_the_exact_permanent_of_repeated_columns():
    # perm_spm(A, m) estimates perm(A with column l repeated m_l times), which
    # is exactly what perm(A, m) computes; the two must be within the bias.
    rng = np.random.default_rng(31)
    a = rng.random((6, 2)) + 0.05
    wide = np.hstack([np.repeat(a[:, [0]], 3, axis=1), np.repeat(a[:, [1]], 3, axis=1)])
    reference = _permanent_ryser(wide)
    estimate = perm_spm(a, [3, 3])
    assert reference < estimate < 1.2 * reference


def test_saddle_point_refuses_what_it_cannot_expand():
    with pytest.raises(ValueError):
        perm_spm(ZERO_ROW)                              # no full support
    with pytest.raises(ValueError):
        perm_spm(np.array([[-1.0, 2.0], [3.0, 4.0]]))   # negative entry
    with pytest.raises(ValueError):
        perm_spm(np.ones((4, 2)))                       # not square, no multiplicities
    with pytest.raises(ValueError):
        perm_spm(np.ones((4, 2)), [1, 1])               # sum(m) != number of rows
    assert perm_spm(np.zeros((0, 0))) == 1.0            # permanent of the empty matrix


def test_each_approximation_refuses_a_structural_zero():
    # An epsilon floor is not invertible: perm(max(A, eps)) = n! eps perm(rest)
    # against a true permanent of 0, and n! outruns eps by n = 18. Every
    # approximate engine must refuse rather than floor.
    assert _permanent_ryser(ZERO_ROW) == 0.0
    with pytest.raises(ValueError):
        perm_bethe(ZERO_ROW)
    with pytest.raises(ValueError):
        perm_heur(ZERO_ROW)
    with pytest.raises(ValueError):
        HuberLawSampler(ZERO_ROW, solve=True)
    with pytest.raises(ValueError):
        AdaPartSampler(ZERO_ROW, solve=True)
    with pytest.raises(ValueError):
        preprocessing_ds(ZERO_ROW)


def test_approximations_terminate_on_a_matrix_with_a_positive_permanent():
    # [[1, 2], [0, 3]] has permanent 3 > 0 but no total support. The samplers
    # used to spin here: a zero row sum makes the guarded Sinkhorn
    # normalization a no-op, the margin error never falls, and the loop has no
    # cap. Termination is what is asserted, not a value.
    assert _permanent_ryser(ZERO_ENTRY) == 3.0
    for factory in (lambda: perm_bethe(ZERO_ENTRY),
                    lambda: perm_heur(ZERO_ENTRY),
                    lambda: HuberLawSampler(ZERO_ENTRY, solve=True),
                    lambda: AdaPartSampler(ZERO_ENTRY, solve=True)):
        start = time.perf_counter()
        with pytest.raises(ValueError):
            factory()
        assert time.perf_counter() - start < 10.0


def test_sinkhorn_scaling_is_capped_and_raises():
    # Total support, not positivity, is what Sinkhorn needs. The block
    # triangular matrix [[J3, 0], [J3, J3]] has permanent 36 > 0 and no total
    # support, so it must be refused rather than iterated indefinitely.
    block = np.block([[np.ones((3, 3)), np.zeros((3, 3))],
                      [np.ones((3, 3)), np.ones((3, 3))]])
    assert _permanent_ryser(block) == 36.0
    with pytest.raises(ValueError):
        preprocessing_ds(block)


def test_maximum_weight_assignment_is_not_greedy():
    # The row-by-row greedy that used to stand in for the assignment takes the
    # larger entry of row 0 and leaves row 1 with the smaller one. alpha3 is
    # that assignment's weight and it sets the Huber-Law flooring level, so the
    # optimum is what the method's guarantee needs.
    a = np.array([[1.0, 2.0], [1e-9, 3.0]])
    sampler = HuberLawSampler(a, seed=0)
    assignment = sampler._hungarian_assignment(np.log(a))
    assert list(assignment) == [0, 1]
    greedy_weight = np.log(a[0, 1]) + np.log(a[1, 0])
    optimal_weight = np.log(a[0, 0]) + np.log(a[1, 1])
    assert optimal_weight > greedy_weight


def test_preprocessing_true_agrees_with_preprocessing_false():
    # preprocessing_ds returns f with DS = X A Y, so perm(A) = perm(DS) / f.
    # _marginal_from_matrix used to multiply, which is wrong by f^2: this state
    # returned 0.03029 against a truth of 0.464.
    demands = np.array([[0.4, 0.6], [0.9, 0.2]])
    network = NetworkNoThink(2, 1, [3], demands, progress=False)
    state = [2, 1]
    without = network.marginal(RyzerPermanent(np.eye(1)), state,
                               preprocessing=False).probability
    with_pre = network.marginal(RyzerPermanent(np.eye(1)), state,
                                preprocessing=True).probability
    assert abs(with_pre - without) < 1e-9 * without


def test_exact_engine_returns_zero_on_a_frobenius_koenig_state():
    # Class 2 never visits station 1, so a state placing more than N_1 jobs at
    # station 1 is structurally impossible: the replicated matrix carries a
    # zero block large enough to force the permanent to 0 by Frobenius-Koenig.
    # The exact engine must say so exactly, not approximately.
    replicated = np.array([
        [1.0, 1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0, 0.0],
        [0.5, 0.5, 2.0, 2.0],
    ])
    assert _permanent_ryser(replicated) == 0.0
    assert permanent(replicated) == 0.0
