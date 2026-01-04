"""
Unit tests for p-norm smoothing ODE system.

Tests the p-norm smoothed capacity constraint and ODE implementation.
"""

import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver.solvers.solver_fld.ode.pnorm import (
    pnorm_smooth,
    pnorm_smooth_vectorized,
    fluid_ode_pnorm,
    PNormODESystem,
)


class TestPNormSmooth(unittest.TestCase):
    """Test p-norm smoothing function."""

    def test_pnorm_smooth_basic(self):
        """Test basic p-norm smoothing computation."""
        # When x << c, g_hat ≈ 1
        result = pnorm_smooth(x=0.1, c=10.0, p=20.0)
        self.assertGreater(result, 0.99)

        # When x >> c, g_hat → (c/x)^(1/p) → 0 as x → ∞
        result = pnorm_smooth(x=100.0, c=1.0, p=20.0)
        self.assertLess(result, 0.1)

        # When x = c, g_hat = 1 / 2^(1/p)
        # For p=20, this is 1 / 2^(1/20) ≈ 0.966
        result = pnorm_smooth(x=1.0, c=1.0, p=20.0)
        self.assertGreater(result, 0.96)
        self.assertLess(result, 1.0)

    def test_pnorm_smooth_bounds(self):
        """Test that p-norm smoothing is bounded in (0, 1]."""
        for x in [0.01, 0.1, 1.0, 10.0, 100.0]:
            for c in [0.1, 1.0, 10.0]:
                result = pnorm_smooth(x, c)
                self.assertGreaterEqual(result, 0.0, f"g_hat({x},{c}) = {result} < 0")
                self.assertLessEqual(result, 1.0, f"g_hat({x},{c}) = {result} > 1")

    def test_pnorm_smooth_monotonicity(self):
        """Test that g_hat is decreasing in x."""
        c = 10.0
        x_vals = [0.1, 1.0, 5.0, 10.0, 50.0, 100.0]
        results = []
        for x in x_vals:
            results.append(pnorm_smooth(x, c))

        # Check monotonic decrease
        for i in range(len(results) - 1):
            self.assertGreaterEqual(results[i], results[i + 1],
                                  f"Not monotonic: g({x_vals[i]}) < g({x_vals[i+1]})")

    def test_pnorm_smooth_p_convergence(self):
        """Test convergence of p-norm as p increases.

        As p → ∞, g_hat(x,c,p) → min(1, c/x).
        For x/c > 1, this means g_hat → c/x < 1.
        """
        x, c = 1.5, 1.0  # x/c = 1.5

        # Lower p: smoother transition
        g_low = pnorm_smooth(x, c, p=2.0)
        # Higher p: sharper transition, should converge toward c/x = 1/1.5
        g_high = pnorm_smooth(x, c, p=100.0)

        # At x/c > 1, higher p should give higher g_hat (closer to c/x)
        self.assertGreater(g_high, g_low)

        # g_high should be close to c/x = 1/1.5 ≈ 0.667
        self.assertAlmostEqual(g_high, 1.0/1.5, places=2)

    def test_pnorm_smooth_edge_cases(self):
        """Test edge cases."""
        # x = 0
        self.assertEqual(pnorm_smooth(0, 1.0), 1.0)

        # c = 0 (invalid)
        self.assertEqual(pnorm_smooth(1.0, 0), 1.0)

        # p = 0 (invalid)
        self.assertEqual(pnorm_smooth(1.0, 1.0, p=0), 1.0)

        # Negative values
        self.assertEqual(pnorm_smooth(-1.0, 1.0), 1.0)
        self.assertEqual(pnorm_smooth(1.0, -1.0), 1.0)


class TestPNormSmoothVectorized(unittest.TestCase):
    """Test vectorized p-norm smoothing."""

    def test_vectorized_vs_scalar(self):
        """Test that vectorized version matches scalar version."""
        x = np.array([0.1, 1.0, 10.0])
        c = np.array([1.0, 1.0, 1.0])
        p = np.array([20.0, 20.0, 20.0])

        result_vec = pnorm_smooth_vectorized(x, c, p)
        result_scalar = np.array([pnorm_smooth(xi, ci, pi) for xi, ci, pi in zip(x, c, p)])

        np.testing.assert_array_almost_equal(result_vec, result_scalar)

    def test_vectorized_default_p(self):
        """Test vectorized with default p parameter."""
        x = np.array([0.1, 1.0, 10.0])
        c = np.array([1.0, 1.0, 1.0])

        result = pnorm_smooth_vectorized(x, c)

        self.assertEqual(len(result), 3)
        self.assertTrue(np.all((result >= 0) & (result <= 1)))


class TestFluidODEPNorm(unittest.TestCase):
    """Test fluid ODE with p-norm smoothing."""

    def setUp(self):
        """Set up test fixtures."""
        # Simple 1D ODE test case
        self.W = np.array([[-1.0]])  # Departure rate
        self.SQ = np.array([[1.0]])  # State-to-queue mapping
        self.Sa = np.array([1.0])    # Server capacity
        self.ALambda = np.array([0.5])  # Arrival rate

    def test_ode_initial_condition(self):
        """Test ODE at initial state."""
        x = np.array([0.0])

        dxdt = fluid_ode_pnorm(0, x, self.W, self.SQ, self.Sa, self.ALambda)

        # At x=0, arrivals should dominate
        self.assertGreater(dxdt[0], 0)

    def test_ode_stability(self):
        """Test ODE for numerical stability."""
        # Test at multiple time points
        for t in [0, 0.1, 1.0, 10.0]:
            for x_val in [0.01, 0.1, 1.0, 10.0]:
                x = np.array([x_val])
                dxdt = fluid_ode_pnorm(t, x, self.W, self.SQ, self.Sa, self.ALambda)

                # Derivative should be finite
                self.assertTrue(np.all(np.isfinite(dxdt)),
                              f"Non-finite derivative at t={t}, x={x_val}: {dxdt}")

    def test_ode_with_pstar(self):
        """Test ODE with explicit p-norm parameters."""
        x = np.array([1.0])
        pstar = np.array([20.0])

        dxdt = fluid_ode_pnorm(0, x, self.W, self.SQ, self.Sa, self.ALambda, pstar)

        self.assertEqual(len(dxdt), 1)
        self.assertTrue(np.isfinite(dxdt[0]))


class TestPNormODESystem(unittest.TestCase):
    """Test PNormODESystem class."""

    def setUp(self):
        """Set up test fixtures."""
        self.W = np.array([[-1.0]])
        self.SQ = np.array([[1.0]])
        self.Sa = np.array([1.0])
        self.ALambda = np.array([0.5])

    def test_pnorm_ode_system_init(self):
        """Test PNormODESystem initialization."""
        system = PNormODESystem(self.W, self.SQ, self.Sa, self.ALambda)

        self.assertIsNotNone(system.pstar)
        self.assertEqual(len(system.pstar), 1)

    def test_pnorm_ode_system_call(self):
        """Test calling ODE system."""
        system = PNormODESystem(self.W, self.SQ, self.Sa, self.ALambda)

        x = np.array([1.0])
        dxdt = system(0, x)

        self.assertEqual(len(dxdt), 1)
        self.assertTrue(np.isfinite(dxdt[0]))

    def test_pnorm_ode_system_custom_pstar(self):
        """Test PNormODESystem with custom p-norm parameter."""
        system = PNormODESystem(self.W, self.SQ, self.Sa, self.ALambda, pstar=50.0)

        self.assertAlmostEqual(system.pstar[0], 50.0)

    def test_pnorm_ode_system_steady_state(self):
        """Test steady state constraint method."""
        system = PNormODESystem(self.W, self.SQ, self.Sa, self.ALambda)
        residual = system.steady_state_constraint()

        # Residual should be callable
        x = np.array([1.0])
        result = residual(x)

        self.assertEqual(len(result), 1)


class TestPNormContinuity(unittest.TestCase):
    """Test continuity properties of p-norm smoothing."""

    def test_continuity_at_boundary(self):
        """Test that p-norm smoothing is continuous."""
        c = 1.0
        x_vals = np.linspace(0.01, 10.0, 100)

        results = [pnorm_smooth(x, c) for x in x_vals]

        # Check for smoothness: no large jumps
        diffs = np.diff(results)
        max_diff = np.max(np.abs(diffs))

        self.assertLess(max_diff, 0.1,
                       f"Large discontinuity detected: {max_diff}")

    def test_smoothness_parameter_variation(self):
        """Test that results vary monotonically as p increases.

        The function g_hat(x,c,p) is monotonic in p when x/c is fixed.
        """
        x, c = 1.5, 1.0

        p_vals = np.linspace(1, 100, 50)
        results = [pnorm_smooth(x, c, p) for p in p_vals]

        # Check monotonicity: as p increases from 1 to 100,
        # g_hat should increase monotonically (for x/c > 1)
        for i in range(len(results) - 1):
            self.assertGreaterEqual(results[i+1], results[i] - 1e-10,
                                   f"Not monotonic at p_vals[{i}]={p_vals[i]}")


if __name__ == '__main__':
    unittest.main(verbosity=2)
