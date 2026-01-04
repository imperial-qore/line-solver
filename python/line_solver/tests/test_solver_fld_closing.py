"""
Unit tests for closing method and FCFS approximation utilities.

Tests the iterative closing method implementation, Coxian fitting,
and moment matching utilities.
"""

import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver.solvers.solver_fld.ode.closing_ode import (
    compute_theta_pnorm,
    compute_theta_softmin,
    compute_theta_statedep,
    ClosingODESystem,
)
from line_solver.solvers.solver_fld.utils.fcfs_approximation import (
    compute_service_moments,
    fit_coxian_2,
    update_service_process,
)


class TestComputeTheta(unittest.TestCase):
    """Test theta computation functions for closing method."""

    def setUp(self):
        """Set up test fixtures."""
        self.x = np.array([1.0, 2.0, 3.0])
        self.c = np.array([2.0, 2.0, 2.0])
        self.SQ = np.eye(3)

    def test_theta_pnorm_bounds(self):
        """Test that p-norm theta is bounded in [0, 1]."""
        theta = compute_theta_pnorm(self.x, self.c, self.SQ, p=20.0)

        self.assertTrue(np.all(theta >= 0.0))
        self.assertTrue(np.all(theta <= 1.0))

    def test_theta_pnorm_basic(self):
        """Test basic p-norm theta computation."""
        # When x << c, theta should be close to 1
        x_low = np.array([0.1, 0.1, 0.1])
        theta = compute_theta_pnorm(x_low, self.c, self.SQ, p=20.0)
        self.assertTrue(np.all(theta > 0.99))

        # When x >> c, theta should be small (approximately c/x)
        # For x=100, c=2: theta ≈ 2/100 = 0.02
        x_high = np.array([100.0, 100.0, 100.0])
        theta = compute_theta_pnorm(x_high, self.c, self.SQ, p=20.0)
        self.assertTrue(np.all(theta < 0.05))  # Allow some margin

    def test_theta_softmin_bounds(self):
        """Test that softmin theta is bounded."""
        theta = compute_theta_softmin(self.x, self.c, self.SQ, alpha=20.0)

        self.assertTrue(np.all(theta >= 0.0))
        self.assertTrue(np.all(theta <= 1.0))

    def test_theta_statedep_bounds(self):
        """Test that state-dependent theta is bounded."""
        theta = compute_theta_statedep(self.x, self.c, self.SQ)

        self.assertTrue(np.all(theta >= 0.0))
        self.assertTrue(np.all(theta <= 1.0))

    def test_theta_methods_consistency(self):
        """Test that all theta methods give reasonable results."""
        theta_pnorm = compute_theta_pnorm(self.x, self.c, self.SQ, p=20.0)
        theta_softmin = compute_theta_softmin(self.x, self.c, self.SQ, alpha=20.0)
        theta_statedep = compute_theta_statedep(self.x, self.c, self.SQ)

        # All should be within reasonable range
        self.assertTrue(np.all(theta_pnorm >= 0))
        self.assertTrue(np.all(theta_softmin >= 0))
        self.assertTrue(np.all(theta_statedep >= 0))


class TestServiceMoments(unittest.TestCase):
    """Test service moment computation."""

    def test_compute_service_moments_basic(self):
        """Test basic moment computation."""
        throughput = np.array([[0.5, 0.3], [0.4, 0.2]])
        utilization = np.array([[0.3, 0.2], [0.2, 0.1]])
        service_rate = np.array([[1.0, 1.0], [1.0, 1.0]])

        mean_S, var_S = compute_service_moments(throughput, utilization, service_rate)

        # Mean service time = utilization / throughput
        expected_mean = np.divide(utilization, throughput)
        np.testing.assert_array_almost_equal(mean_S, expected_mean)

        # Variance should be non-negative
        self.assertTrue(np.all(var_S >= 0))

    def test_moment_computation_stability(self):
        """Test numerical stability of moment computation."""
        # Test with very small values
        throughput = np.array([[1e-6, 1e-6]])
        utilization = np.array([[1e-8, 1e-8]])
        service_rate = np.array([[1.0, 1.0]])

        mean_S, var_S = compute_service_moments(throughput, utilization, service_rate)

        # Should handle small values gracefully
        self.assertTrue(np.all(np.isfinite(mean_S)))
        self.assertTrue(np.all(np.isfinite(var_S)))


class TestCoxianFitting(unittest.TestCase):
    """Test Coxian distribution fitting."""

    def test_fit_coxian_2_exponential(self):
        """Test fitting exponential (SCV=1)."""
        mean = 1.0
        variance = 1.0  # Exponential: var = mean²

        mu_1, mu_2, p = fit_coxian_2(mean, variance)

        # For exponential, p should be 0 (no second phase)
        self.assertEqual(p, 0.0)

    def test_fit_coxian_2_hyperexp(self):
        """Test fitting hyperexponential (SCV > 1)."""
        mean = 1.0
        variance = 5.0  # Hyperexponential

        mu_1, mu_2, p = fit_coxian_2(mean, variance)

        # For SCV > 1, p should be > 0 (use second phase)
        self.assertGreater(p, 0.0)
        self.assertLess(p, 1.0)

    def test_fit_coxian_2_bounds(self):
        """Test that Coxian parameters are bounded."""
        for variance in [0.5, 1.0, 2.0, 5.0, 10.0]:
            mu_1, mu_2, p = fit_coxian_2(1.0, variance)

            # Rates should be positive
            self.assertGreater(mu_1, 0.0)
            self.assertGreater(mu_2, 0.0)

            # Probability should be in [0, 1]
            self.assertGreaterEqual(p, 0.0)
            self.assertLessEqual(p, 1.0)


class TestUpdateServiceProcess(unittest.TestCase):
    """Test service process update."""

    def test_update_exponential_process(self):
        """Test update for exponential service."""
        mean_S = 1.0
        var_S = 1.0  # Exponential
        current_rate = 1.0

        process = update_service_process(mean_S, var_S, current_rate)

        self.assertEqual(process['type'], 'Exp')
        self.assertIn('mu', process['params'])

    def test_update_coxian_process(self):
        """Test update for general (Coxian) service."""
        mean_S = 1.0
        var_S = 2.5  # Slightly variable
        current_rate = 1.0

        process = update_service_process(mean_S, var_S, current_rate)

        self.assertIn(process['type'], ['Coxian2', 'HyperExp'])

    def test_update_process_metadata(self):
        """Test that update includes all metadata."""
        mean_S = 1.0
        var_S = 1.5
        current_rate = 1.0

        process = update_service_process(mean_S, var_S, current_rate)

        # Check required keys
        self.assertIn('type', process)
        self.assertIn('params', process)
        self.assertIn('mean', process)
        self.assertIn('variance', process)
        self.assertIn('scv', process)
        self.assertIn('rate', process)


class TestClosingODESystem(unittest.TestCase):
    """Test ClosingODESystem class."""

    def setUp(self):
        """Set up test fixtures."""
        self.W_T = -np.eye(2)
        self.SQ = np.eye(2)
        self.c = np.array([1.0, 1.0])
        self.A_Lambda = np.array([0.5, 0.3])

    def test_closing_ode_system_pnorm(self):
        """Test ClosingODESystem with p-norm."""
        system = ClosingODESystem(self.W_T, self.SQ, self.c, self.A_Lambda, method='pnorm', p=20.0)

        x = np.array([1.0, 1.0])
        dxdt = system(0, x)

        self.assertEqual(len(dxdt), 2)
        self.assertTrue(np.all(np.isfinite(dxdt)))

    def test_closing_ode_system_softmin(self):
        """Test ClosingODESystem with softmin."""
        system = ClosingODESystem(self.W_T, self.SQ, self.c, self.A_Lambda, method='softmin', alpha=20.0)

        x = np.array([1.0, 1.0])
        dxdt = system(0, x)

        self.assertEqual(len(dxdt), 2)
        self.assertTrue(np.all(np.isfinite(dxdt)))

    def test_closing_ode_system_statedep(self):
        """Test ClosingODESystem with state-dependent."""
        system = ClosingODESystem(self.W_T, self.SQ, self.c, self.A_Lambda, method='statedep')

        x = np.array([1.0, 1.0])
        dxdt = system(0, x)

        self.assertEqual(len(dxdt), 2)
        self.assertTrue(np.all(np.isfinite(dxdt)))


if __name__ == '__main__':
    unittest.main(verbosity=2)
