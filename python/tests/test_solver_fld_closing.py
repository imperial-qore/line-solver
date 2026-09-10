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
)
from line_solver.solvers.solver_fld.ode.closing_ode import (
    compute_theta_statedep,
    ClosingODESystem,
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
