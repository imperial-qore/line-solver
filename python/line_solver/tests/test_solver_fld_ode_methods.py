"""
Unit tests for alternative ODE methods: softmin and statedep.

Tests the softmin smoothing and state-dependent constraint implementations.
"""

import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver.solvers.solver_fld.ode.softmin import (
    softmin_smooth,
    softmin_smooth_vectorized,
    fluid_ode_softmin,
    SoftminODESystem,
)
from line_solver.solvers.solver_fld.ode.statedep import (
    statedep_constraint,
    statedep_constraint_vectorized,
    fluid_ode_statedep,
    StateDepODESystem,
)


class TestSoftminSmooth(unittest.TestCase):
    """Test softmin smoothing function."""

    def test_softmin_smooth_basic(self):
        """Test basic softmin smoothing computation."""
        # When x << c, softmin ≈ x, so g = softmin/x ≈ 1
        result = softmin_smooth(x=0.1, c=10.0, alpha=20.0)
        self.assertGreater(result, 0.99)

        # When x >> c, softmin ≈ c, so g = c/x << 1
        result = softmin_smooth(x=100.0, c=1.0, alpha=20.0)
        self.assertLess(result, 0.02)

        # When x = c, softmin should be between x and c
        result = softmin_smooth(x=1.0, c=1.0, alpha=20.0)
        self.assertGreater(result, 0.5)
        self.assertLess(result, 1.0)

    def test_softmin_smooth_bounds(self):
        """Test that softmin smoothing is bounded in (0, 1]."""
        for x in [0.01, 0.1, 1.0, 10.0, 100.0]:
            for c in [0.1, 1.0, 10.0]:
                result = softmin_smooth(x, c)
                self.assertGreaterEqual(result, 0.0)
                self.assertLessEqual(result, 1.0)

    def test_softmin_smooth_monotonicity(self):
        """Test that softmin is monotonically decreasing in x."""
        c = 10.0
        x_vals = [0.1, 1.0, 5.0, 10.0, 50.0, 100.0]
        results = []
        for x in x_vals:
            results.append(softmin_smooth(x, c))

        # Check monotonic decrease
        for i in range(len(results) - 1):
            self.assertGreaterEqual(results[i], results[i + 1],
                                  f"Not monotonic at x={x_vals[i]}")

    def test_softmin_smooth_alpha_effect(self):
        """Test that alpha parameter affects smoothing behavior."""
        x, c = 1.5, 1.0

        # Lower alpha: smoother transition
        g_low = softmin_smooth(x, c, alpha=2.0)
        # Higher alpha: sharper transition (more like hard min)
        g_high = softmin_smooth(x, c, alpha=100.0)

        # Higher alpha should give higher result (closer to hard min = 1/1.5)
        # when x/c > 1
        self.assertLess(g_low, g_high)


class TestSoftminVectorized(unittest.TestCase):
    """Test vectorized softmin."""

    def test_softmin_vectorized_vs_scalar(self):
        """Test that vectorized matches scalar."""
        x = np.array([0.1, 1.0, 10.0])
        c = np.array([1.0, 1.0, 1.0])

        result_vec = softmin_smooth_vectorized(x, c, 20.0)
        result_scalar = np.array([softmin_smooth(xi, ci, 20.0) for xi, ci in zip(x, c)])

        np.testing.assert_array_almost_equal(result_vec, result_scalar)


class TestStateDepConstraint(unittest.TestCase):
    """Test state-dependent (hard) constraint."""

    def test_statedep_constraint_basic(self):
        """Test basic state-dependent constraint."""
        # When x << c, g ≈ 1
        result = statedep_constraint(x=0.1, c=10.0)
        self.assertAlmostEqual(result, 1.0, places=1)

        # When x >> c, g ≈ c/x
        result = statedep_constraint(x=100.0, c=1.0)
        self.assertAlmostEqual(result, 0.01, places=3)

        # When x = c, g = 1
        result = statedep_constraint(x=1.0, c=1.0)
        self.assertAlmostEqual(result, 1.0)

    def test_statedep_constraint_bounds(self):
        """Test that constraint is bounded in (0, 1]."""
        for x in [0.01, 0.1, 1.0, 10.0, 100.0]:
            for c in [0.1, 1.0, 10.0]:
                result = statedep_constraint(x, c)
                self.assertGreaterEqual(result, 0.0)
                self.assertLessEqual(result, 1.0)

    def test_statedep_constraint_monotonicity(self):
        """Test that constraint decreases with x."""
        c = 10.0
        x_vals = [0.1, 1.0, 5.0, 10.0, 50.0, 100.0]
        results = []
        for x in x_vals:
            results.append(statedep_constraint(x, c))

        # Check monotonic decrease
        for i in range(len(results) - 1):
            self.assertGreaterEqual(results[i], results[i + 1])


class TestStateDepVectorized(unittest.TestCase):
    """Test vectorized state-dependent constraint."""

    def test_statedep_vectorized_vs_scalar(self):
        """Test that vectorized matches scalar."""
        x = np.array([0.1, 1.0, 10.0])
        c = np.array([1.0, 1.0, 1.0])

        result_vec = statedep_constraint_vectorized(x, c)
        result_scalar = np.array([statedep_constraint(xi, ci) for xi, ci in zip(x, c)])

        np.testing.assert_array_almost_equal(result_vec, result_scalar)


class TestFluidODESoftmin(unittest.TestCase):
    """Test fluid ODE with softmin."""

    def setUp(self):
        """Set up test fixtures."""
        self.W = np.array([[-1.0]])
        self.SQ = np.array([[1.0]])
        self.Sa = np.array([1.0])
        self.ALambda = np.array([0.5])

    def test_ode_softmin_stability(self):
        """Test ODE for numerical stability."""
        for x_val in [0.01, 0.1, 1.0, 10.0]:
            x = np.array([x_val])
            dxdt = fluid_ode_softmin(0, x, self.W, self.SQ, self.Sa, self.ALambda)

            self.assertTrue(np.all(np.isfinite(dxdt)))

    def test_ode_softmin_with_alpha(self):
        """Test ODE with different alpha values."""
        x = np.array([1.0])

        for alpha in [2.0, 5.0, 20.0, 100.0]:
            dxdt = fluid_ode_softmin(0, x, self.W, self.SQ, self.Sa, self.ALambda, alpha)
            self.assertTrue(np.isfinite(dxdt[0]))


class TestSoftminODESystem(unittest.TestCase):
    """Test SoftminODESystem class."""

    def setUp(self):
        """Set up test fixtures."""
        self.W = np.array([[-1.0]])
        self.SQ = np.array([[1.0]])
        self.Sa = np.array([1.0])
        self.ALambda = np.array([0.5])

    def test_softmin_system_init(self):
        """Test SoftminODESystem initialization."""
        system = SoftminODESystem(self.W, self.SQ, self.Sa, self.ALambda)
        self.assertAlmostEqual(system.alpha, 20.0)

    def test_softmin_system_call(self):
        """Test calling SoftminODESystem."""
        system = SoftminODESystem(self.W, self.SQ, self.Sa, self.ALambda, alpha=50.0)
        x = np.array([1.0])
        dxdt = system(0, x)

        self.assertEqual(len(dxdt), 1)
        self.assertTrue(np.isfinite(dxdt[0]))


class TestFluidODEStateDep(unittest.TestCase):
    """Test fluid ODE with state-dependent constraint."""

    def setUp(self):
        """Set up test fixtures."""
        self.W = np.array([[-1.0]])
        self.SQ = np.array([[1.0]])
        self.Sa = np.array([1.0])
        self.ALambda = np.array([0.5])

    def test_ode_statedep_stability(self):
        """Test ODE for numerical stability."""
        for x_val in [0.01, 0.1, 1.0, 10.0]:
            x = np.array([x_val])
            dxdt = fluid_ode_statedep(0, x, self.W, self.SQ, self.Sa, self.ALambda)

            self.assertTrue(np.all(np.isfinite(dxdt)))


class TestStateDepODESystem(unittest.TestCase):
    """Test StateDepODESystem class."""

    def setUp(self):
        """Set up test fixtures."""
        self.W = np.array([[-1.0]])
        self.SQ = np.array([[1.0]])
        self.Sa = np.array([1.0])
        self.ALambda = np.array([0.5])

    def test_statedep_system_init(self):
        """Test StateDepODESystem initialization."""
        system = StateDepODESystem(self.W, self.SQ, self.Sa, self.ALambda)
        self.assertIsNotNone(system)

    def test_statedep_system_call(self):
        """Test calling StateDepODESystem."""
        system = StateDepODESystem(self.W, self.SQ, self.Sa, self.ALambda)
        x = np.array([1.0])
        dxdt = system(0, x)

        self.assertEqual(len(dxdt), 1)
        self.assertTrue(np.isfinite(dxdt[0]))


class TestODEMethodComparison(unittest.TestCase):
    """Compare different ODE methods."""

    def setUp(self):
        """Set up test fixtures."""
        self.W = np.array([[-1.0]])
        self.SQ = np.array([[1.0]])
        self.Sa = np.array([1.0])
        self.ALambda = np.array([0.5])

    def test_methods_give_similar_results_at_low_x(self):
        """When x << c, all methods should give similar results."""
        x = np.array([0.1])

        from line_solver.solvers.solver_fld.ode.pnorm import fluid_ode_pnorm

        dxdt_pnorm = fluid_ode_pnorm(0, x, self.W, self.SQ, self.Sa, self.ALambda)
        dxdt_softmin = fluid_ode_softmin(0, x, self.W, self.SQ, self.Sa, self.ALambda)
        dxdt_statedep = fluid_ode_statedep(0, x, self.W, self.SQ, self.Sa, self.ALambda)

        # At low x, all should be close (within 5%)
        np.testing.assert_allclose([dxdt_pnorm], [dxdt_softmin], rtol=0.05)
        np.testing.assert_allclose([dxdt_pnorm], [dxdt_statedep], rtol=0.05)


if __name__ == '__main__':
    unittest.main(verbosity=2)
