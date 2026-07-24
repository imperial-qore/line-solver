"""
Tests for the unified iltcme inverse Laplace transform.
Validates CME, Euler, and Gaver methods against known analytical inversions.
"""

import math
import numpy as np
import pytest

from line_solver.api.lti import ilt, iltcme_load_params, laplace_invert_cme


# F(s) = 1/(1+s) -> f(t) = exp(-t)
def exp_lt(s):
    return 1.0 / (1.0 + s)


# F(s) = 1/(1+s^2) -> f(t) = sin(t)
def sin_lt(s):
    return 1.0 / (1.0 + s * s)


T_POINTS = np.array([0.5, 1.0, 2.0, 5.0])
MAX_FN_EVALS = 25


class TestCme:
    def test_exponential(self):
        result = ilt(exp_lt, T_POINTS, MAX_FN_EVALS, method='cme')
        assert len(result) == len(T_POINTS)
        for i, t in enumerate(T_POINTS):
            expected = math.exp(-t)
            assert abs(result[i] - expected) < 1e-3, \
                f"CME exp(-t) at t={t}: expected={expected:.10f}, got={result[i]:.10f}"

    def test_sine(self):
        result = ilt(sin_lt, T_POINTS, MAX_FN_EVALS, method='cme')
        assert len(result) == len(T_POINTS)
        for i, t in enumerate(T_POINTS):
            expected = math.sin(t)
            assert abs(result[i] - expected) < 0.01, \
                f"CME sin(t) at t={t}: expected={expected:.10f}, got={result[i]:.10f}"

    def test_parameter_loading(self):
        params = iltcme_load_params()
        assert len(params) > 0
        assert 'n' in params[0]
        assert 'cv2' in params[0]
        assert 'mu1' in params[0]

    def test_small_max_fn_evals(self):
        result = ilt(exp_lt, [1.0], 3, method='cme')
        assert len(result) == 1
        assert abs(result[0] - math.exp(-1.0)) < 0.05, \
            "CME with maxFnEvals=3 should still give reasonable result"

    def test_max_fn_evals_constraint(self):
        result_low = ilt(exp_lt, [1.0], 5, method='cme')
        result_high = ilt(exp_lt, [1.0], 100, method='cme')
        err_low = abs(result_low[0] - math.exp(-1.0))
        err_high = abs(result_high[0] - math.exp(-1.0))
        assert err_high <= err_low + 1e-10, \
            f"Higher maxFnEvals should give better accuracy: err5={err_low:.2e}, err100={err_high:.2e}"

    def test_laplace_invert_cme_single(self):
        result = laplace_invert_cme(exp_lt, 1.0, maxFnEvals=25)
        expected = math.exp(-1.0)
        assert abs(result - expected) < 1e-3


class TestEuler:
    def test_exponential(self):
        result = ilt(exp_lt, T_POINTS, MAX_FN_EVALS, method='euler')
        assert len(result) == len(T_POINTS)
        for i, t in enumerate(T_POINTS):
            expected = math.exp(-t)
            assert abs(result[i] - expected) < 1e-6, \
                f"Euler exp(-t) at t={t}: expected={expected:.10f}, got={result[i]:.10f}"

    def test_sine(self):
        result = ilt(sin_lt, T_POINTS, MAX_FN_EVALS, method='euler')
        assert len(result) == len(T_POINTS)
        for i, t in enumerate(T_POINTS):
            expected = math.sin(t)
            assert abs(result[i] - expected) < 1e-6, \
                f"Euler sin(t) at t={t}: expected={expected:.10f}, got={result[i]:.10f}"


class TestGaver:
    def test_exponential(self):
        result = ilt(exp_lt, T_POINTS, MAX_FN_EVALS, method='gaver')
        assert len(result) == len(T_POINTS)
        for i, t in enumerate(T_POINTS):
            expected = math.exp(-t)
            assert abs(result[i] - expected) < 0.15, \
                f"Gaver exp(-t) at t={t}: expected={expected:.10f}, got={result[i]:.10f}"

    def test_sine(self):
        result = ilt(sin_lt, T_POINTS, MAX_FN_EVALS, method='gaver')
        assert len(result) == len(T_POINTS)
        for i, t in enumerate(T_POINTS):
            expected = math.sin(t)
            assert abs(result[i] - expected) < 0.05, \
                f"Gaver sin(t) at t={t}: expected={expected:.10f}, got={result[i]:.10f}"


class TestGeneral:
    def test_default_method_is_cme(self):
        result_default = ilt(exp_lt, T_POINTS, MAX_FN_EVALS)
        result_cme = ilt(exp_lt, T_POINTS, MAX_FN_EVALS, method='cme')
        np.testing.assert_allclose(result_default, result_cme, atol=1e-15)

    def test_invalid_method(self):
        with pytest.raises(ValueError):
            ilt(exp_lt, T_POINTS, MAX_FN_EVALS, method='invalid')

    def test_scalar_input(self):
        result = ilt(exp_lt, 1.0, MAX_FN_EVALS)
        assert len(result) == 1
        assert abs(result[0] - math.exp(-1.0)) < 1e-3


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
