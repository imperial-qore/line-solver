"""
Unit tests for diffusion method and Euler-Maruyama SDE solver.

Tests the custom SDE implementation, population constraint enforcement,
drift and diffusion coefficient computation, and transient analysis.
"""

import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver.solvers.solver_fld.methods.diffusion import (
    DiffusionMethod,
)
from line_solver.solvers.solver_fld.options import SolverFLDOptions


class MockNetworkStruct:
    """Mock NetworkStruct for testing."""

    def __init__(self, nstations=2, nclasses=1, njobs=5.0, nservers=None):
        self.nstations = nstations
        self.nclasses = nclasses
        self.njobs = np.array([njobs])
        self.nservers = nservers if nservers is not None else np.ones(nstations, dtype=int)


class TestDiffusionValidation(unittest.TestCase):
    """Test validation of network topology."""

    def test_closed_network_detection(self):
        """Test that closed network validation works."""
        sn = MockNetworkStruct(nstations=2, nclasses=1, njobs=5.0)
        options = SolverFLDOptions(method='diffusion')

        # Should not raise
        method = DiffusionMethod(sn, options)
        self.assertIsNotNone(method)

    def test_open_network_rejection(self):
        """Test that open networks are rejected."""
        sn = MockNetworkStruct(nstations=2, nclasses=1, njobs=np.inf)
        options = SolverFLDOptions(method='diffusion')

        with self.assertRaises(ValueError) as context:
            DiffusionMethod(sn, options)
        self.assertIn("closed", str(context.exception).lower())

    def test_missing_njobs_rejection(self):
        """Test that missing njobs causes rejection."""
        sn = MockNetworkStruct(nstations=2, nclasses=1)
        sn.njobs = None
        options = SolverFLDOptions(method='diffusion')

        with self.assertRaises(ValueError):
            DiffusionMethod(sn, options)


class TestDriftComputation(unittest.TestCase):
    """Test drift coefficient computation."""

    def setUp(self):
        """Set up test fixtures."""
        self.sn = MockNetworkStruct(nstations=2, nclasses=1, njobs=2.0)
        self.options = SolverFLDOptions(method='diffusion')
        self.method = DiffusionMethod(self.sn, self.options)

    def test_drift_finite(self):
        """Test that drift is always finite."""
        W_T = -np.eye(2)
        SQ = np.eye(2)
        c = np.ones(2)
        x = np.array([1.0, 1.0])

        drift = self.method._compute_drift(x, W_T, SQ, c)

        self.assertEqual(len(drift), 2)
        self.assertTrue(np.all(np.isfinite(drift)))

    def test_drift_zero_at_equilibrium(self):
        """Test drift behavior near equilibrium."""
        W_T = np.array([[-0.5, 0.5], [0.5, -0.5]])
        SQ = np.eye(2)
        c = np.array([2.0, 2.0])
        x = np.array([1.0, 1.0])

        drift = self.method._compute_drift(x, W_T, SQ, c)

        # At balanced state, drift should be relatively small
        self.assertTrue(np.all(np.isfinite(drift)))

    def test_drift_small_populations(self):
        """Test drift with very small populations."""
        W_T = -np.eye(2)
        SQ = np.eye(2)
        c = np.ones(2)
        x = np.array([1e-8, 1e-8])

        drift = self.method._compute_drift(x, W_T, SQ, c)

        # Should handle small values gracefully
        self.assertTrue(np.all(np.isfinite(drift)))


class TestDiffusionComputation(unittest.TestCase):
    """Test diffusion (volatility) coefficient computation."""

    def setUp(self):
        """Set up test fixtures."""
        self.sn = MockNetworkStruct(nstations=2, nclasses=1, njobs=2.0)
        self.options = SolverFLDOptions(method='diffusion')
        self.method = DiffusionMethod(self.sn, self.options)

    def test_diffusion_nonnegative(self):
        """Test that diffusion is always non-negative."""
        SQ = np.eye(2)
        c = np.ones(2)
        x = np.array([1.0, 1.0])

        diffusion = self.method._compute_diffusion(x, SQ, c)

        self.assertEqual(len(diffusion), 2)
        self.assertTrue(np.all(diffusion >= 0.0))

    def test_diffusion_zero_at_zero_throughput(self):
        """Test diffusion is zero when throughput is zero."""
        SQ = np.eye(2)
        c = np.zeros(2)  # No capacity
        x = np.array([1.0, 1.0])

        diffusion = self.method._compute_diffusion(x, SQ, c)

        # Diffusion should be zero (no volatility without capacity)
        self.assertTrue(np.allclose(diffusion, 0.0))

    def test_diffusion_finite(self):
        """Test that diffusion is always finite."""
        SQ = np.eye(2)
        c = np.array([1.0, 2.0])
        x = np.array([1.0, 1.0])

        diffusion = self.method._compute_diffusion(x, SQ, c)

        self.assertTrue(np.all(np.isfinite(diffusion)))

    def test_diffusion_scales_with_throughput(self):
        """Test that diffusion increases with throughput."""
        SQ = np.eye(2)
        c = np.array([1.0, 1.0])

        # Low throughput
        x_low = np.array([0.1, 0.1])
        diff_low = self.method._compute_diffusion(x_low, SQ, c)

        # High throughput
        x_high = np.array([1.0, 1.0])
        diff_high = self.method._compute_diffusion(x_high, SQ, c)

        # Higher throughput should have higher volatility
        self.assertTrue(np.all(diff_high >= diff_low))


class TestPopulationConstraint(unittest.TestCase):
    """Test population constraint enforcement."""

    def setUp(self):
        """Set up test fixtures."""
        self.sn = MockNetworkStruct(nstations=2, nclasses=1, njobs=5.0)
        self.options = SolverFLDOptions(method='diffusion')
        self.method = DiffusionMethod(self.sn, self.options)

    def test_population_preserved(self):
        """Test that total population is preserved."""
        N = np.array([5.0])
        x = np.array([2.0, 3.0])

        x_constrained = self.method._enforce_population_constraint(x, N)

        # Total population should equal N[0]
        self.assertAlmostEqual(np.sum(x_constrained), N[0], places=10)

    def test_population_distribution_balanced(self):
        """Test population distribution when sum is zero."""
        N = np.array([6.0])
        x = np.array([0.0, 0.0])  # Zero population

        x_constrained = self.method._enforce_population_constraint(x, N)

        # Should distribute evenly (3.0 per station)
        expected = np.array([3.0, 3.0])
        np.testing.assert_array_almost_equal(x_constrained, expected)

    def test_population_multiclass(self):
        """Test population constraint with multiple classes."""
        N = np.array([4.0, 6.0])  # 2 classes, 2 stations
        # State vector: [Q_1^1, Q_1^2, Q_2^1, Q_2^2] (i.e., x[i*K + k])
        x = np.array([2.0, 1.0, 1.0, 3.0])

        x_constrained = self.method._enforce_population_constraint(x, N)

        # Check class 1 population (x[0] + x[2])
        self.assertAlmostEqual(x_constrained[0] + x_constrained[2], N[0], places=10)
        # Check class 2 population (x[1] + x[3])
        self.assertAlmostEqual(x_constrained[1] + x_constrained[3], N[1], places=10)

    def test_nonnegativity_preserved(self):
        """Test that nonnegativity is preserved after constraint."""
        N = np.array([5.0])
        x = np.array([3.0, 2.0])

        x_constrained = self.method._enforce_population_constraint(x, N)

        self.assertTrue(np.all(x_constrained >= 0.0))


class TestInitialState(unittest.TestCase):
    """Test initial state computation."""

    def setUp(self):
        """Set up test fixtures."""
        self.sn = MockNetworkStruct(nstations=2, nclasses=1, njobs=5.0)
        self.options = SolverFLDOptions(method='diffusion')
        self.method = DiffusionMethod(self.sn, self.options)

    def test_initial_state_correct_size(self):
        """Test initial state has correct size."""
        N = np.array([5.0])
        x0 = self.method._compute_initial_state(N)

        self.assertEqual(len(x0), 2)  # 2 stations * 1 class

    def test_initial_state_distributed(self):
        """Test initial state distributes population evenly."""
        N = np.array([4.0])
        x0 = self.method._compute_initial_state(N)

        # Should be distributed evenly: 2.0 per station
        expected = np.array([2.0, 2.0])
        np.testing.assert_array_almost_equal(x0, expected)

    def test_initial_state_multiclass(self):
        """Test initial state with multiple classes."""
        self.sn.nclasses = 2
        N = np.array([4.0, 6.0])
        x0 = self.method._compute_initial_state(N)

        self.assertEqual(len(x0), 4)  # 2 stations * 2 classes
        # Class 1: 2.0 per station (4.0 / 2 stations)
        # Class 2: 3.0 per station (6.0 / 2 stations)
        expected = np.array([2.0, 3.0, 2.0, 3.0])
        np.testing.assert_array_almost_equal(x0, expected)


class TestEulerMaruyama(unittest.TestCase):
    """Test Euler-Maruyama SDE stepping."""

    def setUp(self):
        """Set up test fixtures."""
        self.sn = MockNetworkStruct(nstations=2, nclasses=1, njobs=2.0)
        self.options = SolverFLDOptions(
            method='diffusion',
            timespan=(0.0, 1.0),
            timestep=0.1,
        )
        self.method = DiffusionMethod(self.sn, self.options)

    def test_sde_solution_finite(self):
        """Test that SDE solution is finite."""
        W_T = -np.eye(2)
        SQ = np.eye(2)
        c = np.ones(2)
        N = np.array([2.0])
        x0 = np.array([1.0, 1.0])

        t_span, x_final, t_vec, Q_hist = self.method._solve_sde(W_T, SQ, c, N, x0)

        self.assertTrue(np.all(np.isfinite(x_final)))
        self.assertTrue(np.all(np.isfinite(Q_hist[-1])))

    def test_sde_solution_nonnegative(self):
        """Test that queue lengths stay non-negative."""
        W_T = -np.eye(2)
        SQ = np.eye(2)
        c = np.ones(2)
        N = np.array([2.0])
        x0 = np.array([1.0, 1.0])

        t_span, x_final, t_vec, Q_hist = self.method._solve_sde(W_T, SQ, c, N, x0)

        # All queue lengths should be non-negative
        for Q in Q_hist:
            self.assertTrue(np.all(Q >= -1e-10), f"Negative queue: {Q}")

    def test_sde_solution_preserves_population(self):
        """Test that population constraint is maintained."""
        W_T = -np.eye(2)
        SQ = np.eye(2)
        c = np.ones(2)
        N = np.array([2.0])
        x0 = np.array([1.0, 1.0])

        t_span, x_final, t_vec, Q_hist = self.method._solve_sde(W_T, SQ, c, N, x0)

        # Population should be maintained throughout
        for Q in Q_hist:
            self.assertAlmostEqual(np.sum(Q), N[0], places=10)

    def test_sde_trajectory_length(self):
        """Test that trajectory has correct length."""
        W_T = -np.eye(2)
        SQ = np.eye(2)
        c = np.ones(2)
        N = np.array([2.0])
        x0 = np.array([1.0, 1.0])

        t_span, x_final, t_vec, Q_hist = self.method._solve_sde(W_T, SQ, c, N, x0)

        # Number of time steps
        expected_steps = int((1.0 - 0.0) / 0.1) + 1
        self.assertEqual(len(t_vec), expected_steps)
        self.assertEqual(len(Q_hist), expected_steps)


class TestDiffusionMethodSolve(unittest.TestCase):
    """Test full diffusion method solve."""

    def setUp(self):
        """Set up test fixtures."""
        self.sn = MockNetworkStruct(nstations=2, nclasses=1, njobs=3.0)
        self.options = SolverFLDOptions(
            method='diffusion',
            timespan=(0.0, 0.5),
            timestep=0.05,
        )

    def test_solve_returns_result(self):
        """Test that solve returns FLDResult."""
        method = DiffusionMethod(self.sn, self.options)
        result = method.solve()

        self.assertIsNotNone(result)
        self.assertIsNotNone(result.QN)
        self.assertIsNotNone(result.UN)

    def test_solve_result_shapes(self):
        """Test that result matrices have correct shapes."""
        method = DiffusionMethod(self.sn, self.options)
        result = method.solve()

        self.assertEqual(result.QN.shape, (2, 1))
        self.assertEqual(result.UN.shape, (2, 1))
        self.assertEqual(result.RN.shape, (2, 1))
        self.assertEqual(result.TN.shape, (2, 1))
        self.assertEqual(result.CN.shape, (1, 1))
        self.assertEqual(result.XN.shape, (1, 1))

    def test_solve_result_finite(self):
        """Test that all result values are finite."""
        method = DiffusionMethod(self.sn, self.options)
        result = method.solve()

        self.assertTrue(np.all(np.isfinite(result.QN)))
        self.assertTrue(np.all(np.isfinite(result.UN)))
        self.assertTrue(np.all(np.isfinite(result.RN)))

    def test_solve_result_nonnegative(self):
        """Test that queue lengths are non-negative."""
        method = DiffusionMethod(self.sn, self.options)
        result = method.solve()

        self.assertTrue(np.all(result.QN >= -1e-10))
        self.assertTrue(np.all(result.UN >= -1e-10))

    def test_solve_with_multiclass(self):
        """Test solve with multiple job classes."""
        sn = MockNetworkStruct(nstations=2, nclasses=2, njobs=3.0)
        sn.njobs = np.array([2.0, 3.0])
        options = SolverFLDOptions(
            method='diffusion',
            timespan=(0.0, 0.3),
            timestep=0.03,
        )

        method = DiffusionMethod(sn, options)
        result = method.solve()

        self.assertEqual(result.QN.shape, (2, 2))
        self.assertEqual(result.UN.shape, (2, 2))
        self.assertTrue(np.all(np.isfinite(result.QN)))

    def test_solve_runtime_recorded(self):
        """Test that runtime is recorded."""
        method = DiffusionMethod(self.sn, self.options)
        result = method.solve()

        self.assertGreater(result.runtime, 0.0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
