"""
Tests for the Matrix Exponential (ME) and Rational Arrival Process (RAP)
distributions of the native Python implementation.

Mirrors the MATLAB coverage in
line-test.git/test/testsAdvFeatures/distributions/test_me_rap_distributions.m:
construction and validation, rejection of invalid representations, the
negative-entry cases that separate ME/RAP from PH/MAP, sampling accuracy and
CDF evaluation. The cross-language reference values are those asserted by
jar/src/test/java/jline/lang/processes/MECrossLanguageTest.java.

Note: the me_* API modules and tests/test_nc_mem.py are Maximum Entropy
(Kouvatsos), unrelated to the matrix-exponential distributions tested here.
"""
import math

import numpy as np
import pytest

from line_solver import MAP, ME, RAP
from line_solver.distributions import markovian as _markovian
from line_solver.lib.thirdparty.butools.ph import (CheckMEPositiveDensity,
                                                   CheckMERepresentation)

TOL = 1e-9
LOOSE_TOL = 1e-6

# Representations shared with the MATLAB test.
ALPHA_2 = [0.3, 0.7]
A_2 = [[-2.0, 1.5], [0.5, -3.0]]

# ME with negative off-diagonal entries (MATLAB test 6, line 264): a valid ME
# that is not a PH distribution. Eigenvalues of A are -2 and -3.
ALPHA_NEGOFF = [0.5, 0.5]
A_NEGOFF = [[-2.5, -0.5], [-0.5, -2.5]]

# ME with a negative entry in alpha (MATLAB test 6, line 319). The sum of alpha
# is 1, so the representation is a valid ME, but its density is negative near 0.
ALPHA_NEGENTRY = [1.5, -0.5]
A_NEGENTRY = [[-1.0, 1.0], [0.0, -2.0]]

# RAP with a negative off-diagonal entry in H0 (MATLAB test 6, line 300).
H0_NEGOFF = [[-3.0, -0.2], [0.1, -1.5]]
H1_NEGOFF = [[1.7, 1.5], [0.9, 0.5]]


def _raw_moments(alpha, A, k_max):
    """Raw moments m_k = k! * alpha * (-A)^-k * e of an ME distribution."""
    alpha = np.asarray(alpha, dtype=float)
    minv = np.linalg.inv(-np.asarray(A, dtype=float))
    e = np.ones(len(alpha))
    return [float(math.factorial(k) * alpha @ np.linalg.matrix_power(minv, k) @ e)
            for k in range(1, k_max + 1)]


class TestMEConstruction:
    def test_valid_me(self):
        me = ME(ALPHA_2, A_2)
        assert me.getNumberOfPhases() == 2
        assert me.getMean() > 0
        assert me.getVar() > 0
        assert me.getSCV() > 0
        # MECrossLanguageTest.testMEMomentsConsistency: SCV = var / mean^2
        assert me.getSCV() == pytest.approx(me.getVar() / me.getMean() ** 2, abs=TOL)

    def test_getters_return_inputs(self):
        me = ME(ALPHA_2, A_2)
        assert np.allclose(me.getAlpha(), ALPHA_2, atol=TOL)
        assert np.allclose(me.getA(), A_2, atol=TOL)

    def test_from_exp(self):
        rate = 2.0
        me = ME.fromExp(rate)
        assert me.getNumberOfPhases() == 1
        assert me.getMean() == pytest.approx(1.0 / rate, abs=LOOSE_TOL)
        assert me.getVar() == pytest.approx(1.0 / rate ** 2, abs=LOOSE_TOL)
        assert me.getSCV() == pytest.approx(1.0, abs=LOOSE_TOL)

    def test_from_erlang(self):
        k, rate = 3, 1.0
        me = ME.fromErlang(k, rate)
        assert me.getNumberOfPhases() == k
        assert me.getMean() == pytest.approx(k / rate, abs=LOOSE_TOL)
        assert me.getVar() == pytest.approx(k / rate ** 2, abs=LOOSE_TOL)
        assert me.getSCV() == pytest.approx(1.0 / k, abs=LOOSE_TOL)

    def test_from_hyperexp(self):
        p = [0.3, 0.7]
        rates = [1.0, 4.0]
        me = ME.fromHyperExp(p, rates)
        expected_mean = p[0] / rates[0] + p[1] / rates[1]
        m2 = 2.0 * (p[0] / rates[0] ** 2 + p[1] / rates[1] ** 2)
        expected_var = m2 - expected_mean ** 2
        assert me.getMean() == pytest.approx(expected_mean, abs=LOOSE_TOL)
        assert me.getVar() == pytest.approx(expected_var, abs=LOOSE_TOL)
        assert me.getSCV() == pytest.approx(expected_var / expected_mean ** 2,
                                            abs=LOOSE_TOL)

    def test_number_of_phases(self):
        assert ME.fromExp(1.0).getNumberOfPhases() == 1
        assert ME.fromErlang(2, 1.0).getNumberOfPhases() == 2
        assert ME.fromErlang(5, 1.0).getNumberOfPhases() == 5

    def test_process_representation(self):
        # MECrossLanguageTest.testMEProcessRepresentation: D0 = A, D1 = -A*e*alpha
        alpha = np.array([0.4, 0.6])
        A = np.array([[-2.0, 1.0], [0.5, -1.5]])
        me = ME(alpha, A)
        D0, D1 = me.getProcess()
        assert np.allclose(D0, A, atol=TOL)
        expected_D1 = np.outer(-(A @ np.ones(2)), alpha)
        assert np.allclose(D1, expected_D1, atol=LOOSE_TOL)
        assert np.allclose(me.getD0(), D0, atol=TOL)
        assert np.allclose(me.getD1(), D1, atol=TOL)

    def test_negative_off_diagonal_entries_accepted(self):
        # MATLAB test 6 (line 264): valid ME, not a PH distribution.
        eig = np.linalg.eigvals(np.array(A_NEGOFF))
        assert np.all(np.real(eig) < 0) and np.all(np.abs(np.imag(eig)) < TOL)
        me = ME(ALPHA_NEGOFF, A_NEGOFF)
        expected = _raw_moments(ALPHA_NEGOFF, A_NEGOFF, 2)
        assert me.getMean() == pytest.approx(expected[0], abs=LOOSE_TOL)
        assert me.getVar() == pytest.approx(expected[1] - expected[0] ** 2,
                                            abs=LOOSE_TOL)

    def test_negative_alpha_entry_accepted(self):
        # MATLAB test 6 (line 319): alpha has a negative entry but sums to 1.
        assert sum(ALPHA_NEGENTRY) == pytest.approx(1.0, abs=TOL)
        me = ME(ALPHA_NEGENTRY, A_NEGENTRY)
        expected = _raw_moments(ALPHA_NEGENTRY, A_NEGENTRY, 1)
        assert me.getMean() == pytest.approx(expected[0], abs=LOOSE_TOL)


class TestMERejection:
    def test_positive_eigenvalue(self):
        # MATLAB test 1: A = 2 has a positive eigenvalue.
        with pytest.raises(ValueError):
            ME([1.0], [[2.0]])

    def test_alpha_sum_above_one(self):
        with pytest.raises(ValueError):
            ME([0.8, 0.5], A_2)

    def test_alpha_sum_below_zero(self):
        with pytest.raises(ValueError):
            ME([1.5, -1.6], A_NEGENTRY)

    def test_complex_dominant_eigenvalue(self):
        # Eigenvalues -1 +- 5i: negative real parts, but not real.
        A = [[-1.0, 5.0], [-5.0, -1.0]]
        assert np.all(np.real(np.linalg.eigvals(np.array(A))) < 0)
        with pytest.raises(ValueError):
            ME([0.5, 0.5], A)

    def test_dimension_mismatch(self):
        with pytest.raises(ValueError):
            ME([0.3, 0.7], [[-1.0]])

    def test_non_square_matrix(self):
        with pytest.raises(ValueError):
            ME([0.3, 0.7], [[-1.0, 0.5, 0.5], [0.5, -1.0, 0.5]])


class TestMEPositiveDensity:
    """CheckMEPositiveDensity, ported from BuTools (MATLAB and JAR parity)."""

    @staticmethod
    def _check(alpha, A, max_size=1000):
        return bool(CheckMEPositiveDensity(np.matrix(np.atleast_2d(alpha)),
                                           np.matrix(np.array(A, dtype=float)),
                                           max_size, 1e-14))

    def test_phase_type_has_positive_density(self):
        assert self._check([1.0, 0.0, 0.0],
                           [[-1.0, 1.0, 0.0], [0.0, -1.0, 1.0], [0.0, 0.0, -1.0]])
        assert self._check([0.3, 0.7], [[-1.0, 0.0], [0.0, -4.0]])
        assert self._check(ALPHA_2, A_2)

    def test_negative_off_diagonal_me_has_positive_density(self):
        assert self._check(ALPHA_NEGOFF, A_NEGOFF)

    def test_negative_alpha_entry_has_negative_density(self):
        # f(0) = alpha * (-A*e) = -1 < 0, so no Markovian monocyclic equivalent
        # exists, yet the representation is a valid ME.
        alpha = np.array(ALPHA_NEGENTRY)
        A = np.array(A_NEGENTRY)
        assert float(alpha @ (-(A @ np.ones(2)))) < 0
        assert CheckMERepresentation(np.matrix(alpha.reshape(1, 2)),
                                     np.matrix(A), 1e-14)
        assert not self._check(ALPHA_NEGENTRY, A_NEGENTRY)

    def test_invalid_representation_reports_false(self):
        # The check calls CheckMERepresentation first, as MATLAB and the JAR do.
        assert not self._check([1.0], [[2.0]])

    def test_construction_warns_but_succeeds(self, monkeypatch):
        seen = []
        monkeypatch.setattr(_markovian, 'line_warning',
                            lambda caller, msg: seen.append((caller, msg)))
        me = ME(ALPHA_NEGENTRY, A_NEGENTRY)
        assert me.getNumberOfPhases() == 2
        assert len(seen) == 1 and seen[0][0] == 'ME'

    def test_valid_density_does_not_warn(self, monkeypatch):
        seen = []
        monkeypatch.setattr(_markovian, 'line_warning',
                            lambda caller, msg: seen.append((caller, msg)))
        ME(ALPHA_2, A_2)
        ME(ALPHA_NEGOFF, A_NEGOFF)
        assert seen == []


class TestMEFitMoments:
    def test_round_trip_three_moments(self):
        moms = [1.0, 2.5, 10.0]
        me = ME.fitMoments(moms)
        fitted = _raw_moments(me.getAlpha(), me.getA(), 3)
        assert fitted == pytest.approx(moms, rel=1e-6)

    def test_round_trip_hyperexp_moments(self):
        # First three moments of HyperExp(p=[0.3,0.7], rates=[1,4]), which is
        # exactly of order 2, the order implied by three moments.
        p = np.array([0.3, 0.7])
        rates = np.array([1.0, 4.0])
        moms = [float(math.factorial(k) * np.sum(p / rates ** k))
                for k in (1, 2, 3)]
        me = ME.fitMoments(moms)
        assert me.getNumberOfPhases() == 2
        fitted = _raw_moments(me.getAlpha(), me.getA(), 3)
        assert fitted == pytest.approx(moms, rel=1e-6)
        assert me.getMean() == pytest.approx(moms[0], rel=1e-6)

    def test_degenerate_moment_sequence_rejected(self):
        # The moments of Exp(1) are of order 1, so no order-2 ME (which three
        # moments call for) reproduces them. BuTools' van de Liefvoort
        # recursion detects the rank deficiency and refuses to fit.
        with pytest.raises(ValueError):
            ME.fitMoments([1.0, 2.0, 6.0])


class TestMETransforms:
    def test_cdf_matches_exponential(self):
        # MECrossLanguageTest.testMECDFConsistency
        me = ME.fromExp(1.0)
        for t in [0.0, 0.5, 1.0, 2.0, 5.0]:
            assert me.evalCDF(t) == pytest.approx(1.0 - math.exp(-t), abs=LOOSE_TOL)

    def test_lst_at_zero_is_one(self):
        for me in [ME(ALPHA_2, A_2), ME.fromErlang(3, 1.0),
                   ME(ALPHA_NEGOFF, A_NEGOFF)]:
            assert me.evalLST(0.0) == pytest.approx(1.0, abs=1e-12)

    def test_lst_derivative_at_zero_is_minus_mean(self):
        h = 1e-4
        for me in [ME(ALPHA_2, A_2), ME.fromHyperExp([0.3, 0.7], [1.0, 4.0]),
                   ME(ALPHA_NEGOFF, A_NEGOFF)]:
            deriv = (me.evalLST(h) - me.evalLST(-h)) / (2 * h)
            assert -deriv == pytest.approx(me.getMean(), rel=1e-6)

    def test_lst_second_derivative_at_zero_is_second_moment(self):
        h = 1e-3
        me = ME(ALPHA_2, A_2)
        second = (me.evalLST(h) - 2 * me.evalLST(0.0) + me.evalLST(-h)) / h ** 2
        m2 = me.getVar() + me.getMean() ** 2
        assert second == pytest.approx(m2, rel=1e-5)

    def test_pdf_integrates_to_cdf(self):
        me = ME(ALPHA_NEGOFF, A_NEGOFF)
        for upper in [0.2, 0.5, 1.0, 2.0]:
            grid = np.linspace(0.0, upper, 20001)
            pdf = np.array([me.evalPDF(float(x)) for x in grid])
            integral = np.trapezoid(pdf, grid) if hasattr(np, 'trapezoid') \
                else np.trapz(pdf, grid)
            assert integral == pytest.approx(me.evalCDF(upper), abs=1e-6)

    def test_cdf_is_monotone_and_bounded(self):
        me = ME(ALPHA_NEGOFF, A_NEGOFF)
        grid = np.linspace(0.0, 5.0, 200)
        vals = np.array([me.evalCDF(float(x)) for x in grid])
        assert vals[0] == pytest.approx(0.0, abs=TOL)
        assert np.all(np.diff(vals) >= -TOL)
        assert vals[-1] == pytest.approx(1.0, abs=1e-4)


class TestMESampling:
    def test_exponential_moments_at_fixed_seed(self):
        me = ME.fromExp(1.0)
        rng = np.random.default_rng(20260719)
        s = me.sample(20000, rng=rng)
        assert s.shape == (20000,)
        assert np.all(s >= 0)
        assert s.mean() == pytest.approx(me.getMean(), rel=0.03)
        scv = s.var() / s.mean() ** 2
        assert scv == pytest.approx(me.getSCV(), rel=0.05)

    def test_hyperexp_moments_and_ecdf_at_fixed_seed(self):
        me = ME.fromHyperExp([0.3, 0.7], [1.0, 4.0])
        rng = np.random.default_rng(7)
        s = me.sample(20000, rng=rng)
        assert s.mean() == pytest.approx(me.getMean(), rel=0.03)
        assert s.var() / s.mean() ** 2 == pytest.approx(me.getSCV(), rel=0.08)
        # Empirical CDF against the analytic one (Kolmogorov-Smirnov statistic).
        srt = np.sort(s)
        ecdf = np.arange(1, len(srt) + 1) / len(srt)
        analytic = np.array([me.evalCDF(float(x)) for x in srt[::20]])
        ks = np.max(np.abs(ecdf[::20] - analytic))
        assert ks < 0.02

    def test_negative_off_diagonal_me_sampling(self):
        me = ME(ALPHA_NEGOFF, A_NEGOFF)
        rng = np.random.default_rng(11)
        s = me.sample(20000, rng=rng)
        assert s.mean() == pytest.approx(me.getMean(), rel=0.03)
        assert s.var() / s.mean() ** 2 == pytest.approx(me.getSCV(), rel=0.08)

    def test_single_sample_default_shape(self):
        me = ME.fromErlang(2, 1.0)
        s = me.sample(1, rng=np.random.default_rng(3))
        assert s.shape == (1,)
        assert s[0] > 0


class TestRAPConstruction:
    def test_valid_rap(self):
        H0 = [[-2.0, 1.0], [0.5, -1.5]]
        H1 = [[0.5, 0.5], [0.5, 0.5]]
        rap = RAP(H0, H1)
        assert rap.getNumberOfPhases() == 2
        assert rap.getMean() > 0
        assert rap.getVar() > 0
        assert rap.getSCV() == pytest.approx(rap.getVar() / rap.getMean() ** 2,
                                             abs=TOL)
        # MECrossLanguageTest.testRAPProcessRepresentation: D0 = H0, D1 = H1
        assert np.allclose(rap.getD0(), H0, atol=TOL)
        assert np.allclose(rap.getD1(), H1, atol=TOL)
        assert np.allclose(rap.getH0(), H0, atol=TOL)
        assert np.allclose(rap.getH1(), H1, atol=TOL)

    def test_from_poisson(self):
        rate = 2.0
        rap = RAP.fromPoisson(rate)
        assert rap.getNumberOfPhases() == 1
        assert rap.getMean() == pytest.approx(1.0 / rate, abs=LOOSE_TOL)
        assert rap.getVar() == pytest.approx(1.0 / rate ** 2, abs=LOOSE_TOL)
        assert rap.getSCV() == pytest.approx(1.0, abs=LOOSE_TOL)

    def test_from_erlang(self):
        k, rate = 2, 1.0
        rap = RAP.fromErlang(k, rate)
        assert rap.getNumberOfPhases() == k
        assert rap.getMean() == pytest.approx(k / rate, abs=LOOSE_TOL)
        assert rap.getVar() == pytest.approx(k / rate ** 2, abs=LOOSE_TOL)
        assert rap.getSCV() == pytest.approx(1.0 / k, abs=LOOSE_TOL)
        assert RAP.fromErlang(4, 1.0).getNumberOfPhases() == 4

    def test_from_map(self):
        D0 = np.array([[-2.0, 1.0], [0.5, -1.5]])
        D1 = np.array([[0.5, 0.5], [0.5, 0.5]])
        m = MAP(D0, D1)
        rap = RAP.fromMAP(m)
        assert rap.getMean() == pytest.approx(m.getMean(), abs=LOOSE_TOL)
        assert rap.getVar() == pytest.approx(m.getVar(), abs=LOOSE_TOL)
        assert rap.getSCV() == pytest.approx(m.getSCV(), abs=LOOSE_TOL)

    def test_negative_off_diagonal_entries_accepted(self):
        # MATLAB test 6 (line 300): valid RAP, not a MAP.
        H0 = np.array(H0_NEGOFF)
        H1 = np.array(H1_NEGOFF)
        assert np.allclose((H0 + H1).sum(axis=1), 0.0, atol=TOL)
        rap = RAP(H0, H1)
        assert rap.getMean() > 0
        assert rap.getSCV() > 0
        assert rap.getD0()[0, 1] < 0


class TestRAPRejection:
    def test_nonzero_row_sums(self):
        # MATLAB test 2: row 1 of H0+H1 sums to -0.5.
        H0 = [[-2.0, 0.5], [0.5, -1.5]]
        H1 = [[1.0, 0.0], [0.5, 0.5]]
        with pytest.raises(ValueError):
            RAP(H0, H1)

    def test_positive_eigenvalue(self):
        with pytest.raises(ValueError):
            RAP([[1.0]], [[-1.0]])

    def test_size_mismatch(self):
        with pytest.raises(ValueError):
            RAP([[-1.0, 1.0], [0.5, -1.5]], [[1.0]])


class TestRAPSampling:
    def test_renewal_process_is_uncorrelated(self):
        # An Erlang renewal RAP has zero autocorrelation; the sampler must not
        # manufacture any.
        rap = RAP.fromErlang(2, 1.0)
        s = rap.sample(20000, rng=np.random.default_rng(101))
        assert s.mean() == pytest.approx(rap.getMean(), rel=0.03)
        x = s - s.mean()
        acf1 = float((x[:-1] * x[1:]).mean() / x.var())
        assert abs(acf1) < 0.03

    def test_lag1_autocorrelation_is_reproduced(self):
        # Strongly correlated MMPP2 read as a RAP: the conditional phase-vector
        # update is what makes the sampled lag-1 autocorrelation match the
        # analytic one. Drawing from the marginal alone would give acf1 = 0.
        q = 0.1
        H0 = np.array([[-(10.0 + q), q], [q, -(1.0 + q)]])
        H1 = np.array([[10.0, 0.0], [0.0, 1.0]])
        rap = RAP(H0, H1)
        acf = rap.getACF([1, 2])
        assert acf[0] == pytest.approx(0.35355123, abs=1e-6)

        s = rap.sample(50000, rng=np.random.default_rng(2024))
        assert s.mean() == pytest.approx(rap.getMean(), rel=0.05)
        x = s - s.mean()
        acf1 = float((x[:-1] * x[1:]).mean() / x.var())
        acf2 = float((x[:-2] * x[2:]).mean() / x.var())
        assert acf1 == pytest.approx(acf[0], abs=0.03)
        assert acf2 == pytest.approx(acf[1], abs=0.03)

    def test_negative_off_diagonal_rap_sampling(self):
        rap = RAP(H0_NEGOFF, H1_NEGOFF)
        s = rap.sample(20000, rng=np.random.default_rng(5))
        assert np.all(s > 0)
        assert s.mean() == pytest.approx(rap.getMean(), rel=0.05)


class TestCrossLanguageReferenceValues:
    """Values produced by MATLAB ME.m / RAP.m and by the JAR, to 15 digits.

    The MATLAB figures were obtained by running ME/RAP on the representations
    of MECrossLanguageTest under matlab -batch; Python reproduces them exactly.
    """

    def test_me_moments(self):
        me = ME(ALPHA_2, A_2)
        assert me.getMean() == pytest.approx(0.59047619047619, abs=1e-13)
        assert me.getVar() == pytest.approx(0.395102040816326, abs=1e-13)
        assert me.getSCV() == pytest.approx(1.13319458896982, abs=1e-13)

    def test_me_hyperexp_moments(self):
        me = ME.fromHyperExp([0.3, 0.7], [1.0, 4.0])
        assert me.getMean() == pytest.approx(0.475, abs=1e-13)
        assert me.getVar() == pytest.approx(0.461875, abs=1e-13)
        assert me.getSCV() == pytest.approx(2.04709141274238, abs=1e-13)

    def test_me_cdf_values(self):
        me = ME.fromExp(1.0)
        expected = {0.0: 0.0, 0.5: 0.393469340287367, 1.0: 0.632120558828558,
                    2.0: 0.864664716763387, 5.0: 0.993262053000915}
        for t, ref in expected.items():
            assert me.evalCDF(t) == pytest.approx(ref, abs=1e-13)

    def test_rap_moments(self):
        rap = RAP([[-2.0, 1.0], [0.5, -1.5]], [[0.5, 0.5], [0.5, 0.5]])
        assert rap.getMean() == pytest.approx(1.0, abs=1e-13)
        assert rap.getVar() == pytest.approx(1.0, abs=1e-13)
        assert rap.getSCV() == pytest.approx(1.0, abs=1e-13)

    def test_rap_negative_off_diagonal_moments_and_acf(self):
        rap = RAP(H0_NEGOFF, H1_NEGOFF)
        assert rap.getMean() == pytest.approx(0.458167330677291, abs=1e-13)
        assert rap.getSCV() == pytest.approx(1.23803469561872, abs=1e-13)
        assert rap.getACF(1)[0] == pytest.approx(-0.0106343027915219, abs=1e-13)

    def test_positive_density_flags_match_matlab(self):
        # MATLAB CheckMEPositiveDensity(.,.,1000) returns 1, 0, 1, 0 here.
        check = TestMEPositiveDensity._check
        assert check(ALPHA_NEGOFF, A_NEGOFF) is True
        assert check(ALPHA_NEGENTRY, A_NEGENTRY) is False
        assert check([1.0, 0.0, 0.0],
                     [[-1.0, 1.0, 0.0], [0.0, -1.0, 1.0], [0.0, 0.0, -1.0]]) is True
        assert check([1.0], [[2.0]]) is False


class TestRAPTransforms:
    def test_lst_at_zero_and_mean(self):
        rap = RAP(H0_NEGOFF, H1_NEGOFF)
        assert rap.evalLST(0.0) == pytest.approx(1.0, abs=1e-12)
        h = 1e-4
        deriv = (rap.evalLST(h) - rap.evalLST(-h)) / (2 * h)
        assert -deriv == pytest.approx(rap.getMean(), rel=1e-6)

    def test_cdf_matches_exponential(self):
        rap = RAP.fromPoisson(1.0)
        for t in [0.0, 0.5, 1.0, 2.0, 5.0]:
            assert rap.evalCDF(t) == pytest.approx(1.0 - math.exp(-t),
                                                   abs=LOOSE_TOL)
