"""
Unit tests for MFQ method and Erlang-C analytical solvers.

Tests the M/M/c analytical solution, Erlang-C formula,
parameter extraction, and single-queue validation.
"""

import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver.solvers.solver_fld.methods.mfq import (
    MFQMethod,
)
from line_solver.solvers.solver_fld.options import SolverFLDOptions


class MockNetworkStruct:
    """Mock NetworkStruct for testing."""

    def __init__(self, nstations=1, nclasses=1, lambda_arr=None, rates=None, nservers=None):
        self.nstations = nstations
        self.nclasses = nclasses
        self.lambda_arr = lambda_arr if lambda_arr is not None else np.array([0.5])
        self.rates = rates if rates is not None else np.array([[1.0]])
        self.nservers = nservers if nservers is not None else np.array([1])


class TestMFQValidation(unittest.TestCase):
    """Test validation of single-queue topology."""

    def test_single_queue_detection(self):
        """Test that single-queue topology is accepted."""
        sn = MockNetworkStruct(nstations=1, nclasses=1)
        options = SolverFLDOptions(method='mfq')

        # Should not raise
        method = MFQMethod(sn, options)
        self.assertIsNotNone(method)

    def test_multiqueue_rejection(self):
        """Test that multi-queue networks are rejected."""
        sn = MockNetworkStruct(nstations=2, nclasses=1)
        options = SolverFLDOptions(method='mfq')

        with self.assertRaises(ValueError) as context:
            MFQMethod(sn, options)
        self.assertIn("single", str(context.exception).lower())

    def test_missing_nstations_rejection(self):
        """Test that missing nstations causes rejection."""
        sn = MockNetworkStruct(nstations=1, nclasses=1)
        delattr(sn, 'nstations')
        options = SolverFLDOptions(method='mfq')

        with self.assertRaises(ValueError):
            MFQMethod(sn, options)


class TestParameterExtraction(unittest.TestCase):
    """Test extraction of queue parameters from NetworkStruct."""

    def setUp(self):
        """Set up test fixtures."""
        self.sn = MockNetworkStruct(nstations=1, nclasses=1)
        self.options = SolverFLDOptions(method='mfq')
        self.method = MFQMethod(self.sn, self.options)

    def test_extract_single_class(self):
        """Test parameter extraction with single class."""
        sn = MockNetworkStruct(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([0.5]),
            rates=np.array([[1.0]]),
            nservers=np.array([1]),
        )
        method = MFQMethod(sn, self.options)
        params = method._extract_queue_parameters()

        self.assertEqual(params['nclasses'], 1)
        self.assertEqual(params['nservers'], 1)
        self.assertAlmostEqual(params['lambda_arr'][0], 0.5)
        self.assertAlmostEqual(params['mu_arr'][0], 1.0)

    def test_extract_multiclass(self):
        """Test parameter extraction with multiple classes."""
        sn = MockNetworkStruct(
            nstations=1,
            nclasses=2,
            lambda_arr=np.array([0.5, 0.3]),
            rates=np.array([[1.0, 1.2]]),
            nservers=np.array([2]),
        )
        method = MFQMethod(sn, self.options)
        params = method._extract_queue_parameters()

        self.assertEqual(params['nclasses'], 2)
        self.assertEqual(params['nservers'], 2)
        self.assertEqual(len(params['lambda_arr']), 2)
        self.assertEqual(len(params['mu_arr']), 2)

    def test_extract_defaults(self):
        """Test that defaults are used when parameters missing."""
        sn = MockNetworkStruct(nstations=1)
        sn.lambda_arr = None
        sn.rates = None
        sn.nservers = None
        method = MFQMethod(sn, self.options)
        params = method._extract_queue_parameters()

        self.assertIsNotNone(params['lambda_arr'])
        self.assertIsNotNone(params['mu_arr'])
        self.assertEqual(params['nservers'], 1)


class TestErlangC(unittest.TestCase):
    """Test Erlang-C formula implementation."""

    def setUp(self):
        """Set up test fixtures."""
        self.sn = MockNetworkStruct(nstations=1, nclasses=1)
        self.options = SolverFLDOptions(method='mfq')
        self.method = MFQMethod(self.sn, self.options)

    def test_erlang_c_bounds(self):
        """Test that Erlang-C is bounded in [0, 1]."""
        for lambda_val in [0.1, 0.5, 0.9]:
            pw = self.method._erlang_c(lambda_val, 1.0, 1)
            self.assertTrue(0.0 <= pw <= 1.0)

    def test_erlang_c_exponential(self):
        """Test Erlang-C for M/M/1 (single server)."""
        # For M/M/1: Pw = rho where rho = lambda/mu
        lambda_val = 0.5
        mu = 1.0
        c = 1
        pw = self.method._erlang_c(lambda_val, mu, c)

        # Expected: rho = 0.5
        self.assertAlmostEqual(pw, 0.5, places=5)

    def test_erlang_c_multiserver(self):
        """Test Erlang-C for M/M/c (multiple servers)."""
        # For M/M/2 with low load
        lambda_val = 0.5
        mu = 1.0
        c = 2
        pw = self.method._erlang_c(lambda_val, mu, c)

        # Should be less than single server
        pw_single = self.method._erlang_c(lambda_val, mu, 1)
        self.assertLess(pw, pw_single)

    def test_erlang_c_increasing_in_load(self):
        """Test that Erlang-C increases with traffic intensity."""
        mu = 1.0
        c = 1

        pw_low = self.method._erlang_c(0.3, mu, c)
        pw_mid = self.method._erlang_c(0.5, mu, c)
        pw_high = self.method._erlang_c(0.7, mu, c)

        self.assertLess(pw_low, pw_mid)
        self.assertLess(pw_mid, pw_high)

    def test_erlang_c_decreasing_in_servers(self):
        """Test that Erlang-C decreases with more servers."""
        lambda_val = 1.0
        mu = 1.0

        pw_1 = self.method._erlang_c(lambda_val, mu, 1)
        pw_2 = self.method._erlang_c(lambda_val, mu, 2)
        pw_3 = self.method._erlang_c(lambda_val, mu, 3)

        self.assertLess(pw_2, pw_1)
        self.assertLess(pw_3, pw_2)

    def test_erlang_c_unstable(self):
        """Test Erlang-C behavior at/above stability limit."""
        # When A >= c, system is unstable
        pw = self.method._erlang_c(2.0, 1.0, 1)  # A=2 >= c=1
        self.assertEqual(pw, 1.0)


class TestMMQueueSolver(unittest.TestCase):
    """Test M/M/c queue analytical solver."""

    def setUp(self):
        """Set up test fixtures."""
        self.sn = MockNetworkStruct(nstations=1, nclasses=1)
        self.options = SolverFLDOptions(method='mfq')
        self.method = MFQMethod(self.sn, self.options)

    def test_mm1_basic(self):
        """Test M/M/1 queue solution."""
        params = {
            'lambda_arr': [0.5],
            'mu_arr': [1.0],
            'nservers': 1,
            'nclasses': 1,
        }

        QN, UN, RN, TN = self.method._solve_mm_queue(params)

        # Check shapes
        self.assertEqual(QN.shape, (1, 1))
        self.assertEqual(UN.shape, (1, 1))
        self.assertEqual(RN.shape, (1, 1))
        self.assertEqual(TN.shape, (1, 1))

        # For M/M/1 with rho=0.5:
        # L = rho / (1 - rho) = 0.5 / 0.5 = 1
        # W = 1 / (mu - lambda) = 1 / 0.5 = 2
        expected_q = 0.5 / (1.0 - 0.5)  # ~1.0
        expected_r = 1.0 / (1.0 - 0.5)  # ~2.0

        self.assertAlmostEqual(QN[0, 0], expected_q, places=2)
        self.assertAlmostEqual(RN[0, 0], expected_r, places=2)
        self.assertAlmostEqual(TN[0, 0], 0.5, places=5)

    def test_mmc_multiserver(self):
        """Test M/M/c queue with multiple servers."""
        params = {
            'lambda_arr': [1.0],
            'mu_arr': [1.0],
            'nservers': 2,
            'nclasses': 1,
        }

        QN, UN, RN, TN = self.method._solve_mm_queue(params)

        # With 2 servers and A=1: rho = 1/2 = 0.5
        # Should be stable
        self.assertFalse(np.isinf(QN[0, 0]))
        self.assertFalse(np.isinf(RN[0, 0]))

    def test_mmc_unstable(self):
        """Test unstable M/M/c queue."""
        params = {
            'lambda_arr': [2.0],
            'mu_arr': [1.0],
            'nservers': 1,  # Insufficient capacity
            'nclasses': 1,
        }

        QN, UN, RN, TN = self.method._solve_mm_queue(params)

        # Should be unstable
        self.assertTrue(np.isinf(QN[0, 0]))
        self.assertTrue(np.isinf(RN[0, 0]))

    def test_mmc_multiclass(self):
        """Test M/M/c with multiple job classes."""
        params = {
            'lambda_arr': [0.5, 0.3],
            'mu_arr': [1.0, 1.0],
            'nservers': 1,
            'nclasses': 2,
        }

        QN, UN, RN, TN = self.method._solve_mm_queue(params)

        # Check shapes
        self.assertEqual(QN.shape, (1, 2))

        # Both should be stable (total lambda = 0.8 < mu = 1.0)
        self.assertFalse(np.any(np.isinf(QN)))

    def test_utilization_bounds(self):
        """Test that utilization is bounded in [0, 1]."""
        params = {
            'lambda_arr': [0.5],
            'mu_arr': [1.0],
            'nservers': 1,
            'nclasses': 1,
        }

        QN, UN, RN, TN = self.method._solve_mm_queue(params)

        self.assertGreaterEqual(UN[0, 0], 0.0)
        self.assertLessEqual(UN[0, 0], 1.0)

    def test_queue_length_nonnegative(self):
        """Test that queue lengths are non-negative."""
        params = {
            'lambda_arr': [0.3],
            'mu_arr': [1.0],
            'nservers': 2,
            'nclasses': 1,
        }

        QN, UN, RN, TN = self.method._solve_mm_queue(params)

        self.assertGreaterEqual(QN[0, 0], 0.0)

    def test_throughput_preserved(self):
        """Test that throughput equals arrival rate."""
        params = {
            'lambda_arr': [0.5],
            'mu_arr': [1.0],
            'nservers': 1,
            'nclasses': 1,
        }

        QN, UN, RN, TN = self.method._solve_mm_queue(params)

        self.assertAlmostEqual(TN[0, 0], 0.5, places=5)


class TestMFQMethodSolve(unittest.TestCase):
    """Test full MFQ method solve."""

    def setUp(self):
        """Set up test fixtures."""
        self.sn = MockNetworkStruct(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([0.5]),
            rates=np.array([[1.0]]),
            nservers=np.array([1]),
        )
        self.options = SolverFLDOptions(method='mfq')

    def test_solve_returns_result(self):
        """Test that solve returns FLDResult."""
        method = MFQMethod(self.sn, self.options)
        result = method.solve()

        self.assertIsNotNone(result)
        self.assertIsNotNone(result.QN)
        self.assertIsNotNone(result.UN)

    def test_solve_result_shapes(self):
        """Test that result matrices have correct shapes."""
        method = MFQMethod(self.sn, self.options)
        result = method.solve()

        self.assertEqual(result.QN.shape, (1, 1))
        self.assertEqual(result.UN.shape, (1, 1))
        self.assertEqual(result.RN.shape, (1, 1))
        self.assertEqual(result.TN.shape, (1, 1))
        self.assertEqual(result.CN.shape, (1, 1))
        self.assertEqual(result.XN.shape, (1, 1))

    def test_solve_result_finite(self):
        """Test that result values are finite."""
        method = MFQMethod(self.sn, self.options)
        result = method.solve()

        self.assertTrue(np.all(np.isfinite(result.QN)))
        self.assertTrue(np.all(np.isfinite(result.UN)))
        self.assertTrue(np.all(np.isfinite(result.TN)))

    def test_solve_with_multiclass(self):
        """Test solve with multiple job classes."""
        sn = MockNetworkStruct(
            nstations=1,
            nclasses=2,
            lambda_arr=np.array([0.5, 0.3]),
            rates=np.array([[1.0, 1.0]]),
            nservers=np.array([1]),
        )
        options = SolverFLDOptions(method='mfq')

        method = MFQMethod(sn, options)
        result = method.solve()

        self.assertEqual(result.QN.shape, (1, 2))
        self.assertTrue(np.all(np.isfinite(result.QN)))

    def test_solve_with_multiserver(self):
        """Test solve with multiple servers."""
        sn = MockNetworkStruct(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([1.0]),
            rates=np.array([[1.0]]),
            nservers=np.array([2]),
        )
        options = SolverFLDOptions(method='mfq')

        method = MFQMethod(sn, options)
        result = method.solve()

        self.assertEqual(result.QN.shape, (1, 1))
        self.assertFalse(np.isinf(result.QN[0, 0]))

    def test_solve_runtime_recorded(self):
        """Test that runtime is recorded."""
        method = MFQMethod(self.sn, self.options)
        result = method.solve()

        self.assertGreater(result.runtime, 0.0)

    def test_solve_method_identifier(self):
        """Test that method identifier is set correctly."""
        method = MFQMethod(self.sn, self.options)
        result = method.solve()

        self.assertEqual(result.method, 'mfq')


if __name__ == '__main__':
    unittest.main(verbosity=2)
