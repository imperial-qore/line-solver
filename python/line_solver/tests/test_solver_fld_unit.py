"""
Unit tests for SolverFLD core infrastructure and data structures.

This file contains unit tests for:
- SolverFLD initialization and configuration
- FLDResult and SolverFLDOptions data structures
- Static methods (listValidMethods, supports, getFeatureSet)
- Result accessor methods

Comprehensive integration tests are in test_solver_fld_integration.py
"""

import unittest
import numpy as np
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver.solvers.solver_fld import SolverFLD
from line_solver.solvers.solver_fld.options import SolverFLDOptions, FLDResult


class SimpleNetworkStruct:
    """Minimal mock NetworkStruct for unit tests."""

    def __init__(self, nstations=2, nclasses=1):
        self.nstations = nstations
        self.nclasses = nclasses
        self.njobs = []
        self.nodetype = np.zeros(nstations, dtype=int)
        self.stationToNode = np.arange(nstations, dtype=int)
        self.rates = np.ones((nstations, nclasses))
        self.scv = np.ones((nstations, nclasses))
        self.lambda_arr = np.ones(nclasses) * 0.5

    def compileStruct(self):
        """Return self for testing."""
        return self


class TestFLDResultStructure(unittest.TestCase):
    """Test FLDResult data structure."""

    def test_fldresult_creation(self):
        """Test FLDResult can be created with required fields."""
        M, K = 3, 2
        QN = np.random.rand(M, K)
        UN = np.random.rand(M, K)
        RN = np.random.rand(M, K)
        TN = np.random.rand(M, K)
        CN = np.random.rand(1, K)
        XN = np.random.rand(1, K)

        result = FLDResult(
            QN=QN, UN=UN, RN=RN, TN=TN,
            CN=CN, XN=XN,
            iterations=10, method='matrix', runtime=0.5
        )

        self.assertEqual(result.method, 'matrix')
        self.assertEqual(result.iterations, 10)
        self.assertEqual(result.runtime, 0.5)
        np.testing.assert_array_equal(result.QN, QN)
        np.testing.assert_array_equal(result.UN, UN)

    def test_fldresult_optional_fields(self):
        """Test FLDResult with optional fields."""
        M, K = 2, 1
        result = FLDResult(
            QN=np.zeros((M, K)),
            UN=np.zeros((M, K)),
            RN=np.zeros((M, K)),
            TN=np.zeros((M, K)),
            CN=np.zeros((1, K)),
            XN=np.zeros((1, K))
        )

        self.assertIsNone(result.t)
        self.assertIsNone(result.xvec)
        self.assertEqual(len(result.QNt), 0)
        self.assertEqual(len(result.UNt), 0)


class TestSolverFLDOptions(unittest.TestCase):
    """Test SolverFLDOptions data structure."""

    def test_default_options(self):
        """Test default option values."""
        opts = SolverFLDOptions()

        self.assertEqual(opts.method, 'default')
        self.assertEqual(opts.tol, 1e-4)
        self.assertEqual(opts.iter_max, 200)
        self.assertEqual(opts.iter_tol, 1e-4)
        self.assertTrue(opts.stiff)
        self.assertFalse(opts.verbose)
        self.assertEqual(opts.pstar, 20.0)
        self.assertEqual(opts.softmin_alpha, 20.0)

    def test_custom_options(self):
        """Test custom option values."""
        opts = SolverFLDOptions(
            method='matrix',
            tol=1e-6,
            verbose=True,
            pstar=50.0
        )

        self.assertEqual(opts.method, 'matrix')
        self.assertEqual(opts.tol, 1e-6)
        self.assertTrue(opts.verbose)
        self.assertEqual(opts.pstar, 50.0)


class TestSolverFLDInitialization(unittest.TestCase):
    """Test SolverFLD initialization and configuration."""

    def test_solver_initialization_default(self):
        """Test SolverFLD initialization with defaults."""
        sn = SimpleNetworkStruct(nstations=2, nclasses=1)
        solver = SolverFLD(sn)

        self.assertIsNotNone(solver)
        self.assertEqual(solver.options.method, 'default')
        self.assertIsNone(solver.result)

    def test_solver_with_method_param(self):
        """Test SolverFLD with method parameter."""
        sn = SimpleNetworkStruct(nstations=2, nclasses=1)
        solver = SolverFLD(sn, method='matrix')

        self.assertEqual(solver.options.method, 'matrix')

    def test_solver_with_custom_options(self):
        """Test SolverFLD with custom options."""
        sn = SimpleNetworkStruct(nstations=2, nclasses=1)
        opts = SolverFLDOptions(method='matrix', tol=1e-5, verbose=True)
        solver = SolverFLD(sn, options=opts)

        self.assertEqual(solver.options.tol, 1e-5)
        self.assertTrue(solver.options.verbose)

    def test_method_parameter_overrides_options(self):
        """Test that method parameter overrides options method."""
        sn = SimpleNetworkStruct(nstations=2, nclasses=1)
        opts = SolverFLDOptions(method='softmin')
        solver = SolverFLD(sn, method='matrix', options=opts)

        self.assertEqual(solver.options.method, 'matrix')


class TestSolverFLDStaticMethods(unittest.TestCase):
    """Test SolverFLD static methods."""

    def test_list_valid_methods(self):
        """Test listValidMethods returns all supported methods."""
        methods = SolverFLD.listValidMethods()

        self.assertIn('default', methods)
        self.assertIn('matrix', methods)
        self.assertIn('pnorm', methods)
        self.assertIn('softmin', methods)
        self.assertIn('statedep', methods)
        self.assertIn('closing', methods)
        self.assertIn('diffusion', methods)
        self.assertIn('mfq', methods)
        self.assertIn('butools', methods)

    def test_supports_valid_method(self):
        """Test supports() for valid methods."""
        sn = SimpleNetworkStruct(nstations=2, nclasses=1)

        can_solve, reason = SolverFLD.supports(sn, 'matrix')
        self.assertTrue(can_solve)
        self.assertIsNone(reason)

    def test_supports_unknown_method(self):
        """Test supports() for unknown methods."""
        sn = SimpleNetworkStruct(nstations=2, nclasses=1)

        can_solve, reason = SolverFLD.supports(sn, 'nonexistent')
        self.assertFalse(can_solve)
        self.assertIsNotNone(reason)

    def test_get_feature_set(self):
        """Test getFeatureSet returns expected features."""
        features = SolverFLD.getFeatureSet()

        self.assertIn('open_networks', features)
        self.assertIn('closed_networks', features)
        self.assertIn('mixed_networks', features)
        self.assertTrue(features['open_networks'])
        self.assertTrue(features['closed_networks'])
        self.assertTrue(features['mixed_networks'])
        self.assertIn('FCFS', features['scheduling'])
        self.assertIn('PS', features['scheduling'])

    def test_default_options(self):
        """Test defaultOptions returns SolverFLDOptions."""
        opts = SolverFLD.defaultOptions()

        self.assertIsInstance(opts, SolverFLDOptions)
        self.assertEqual(opts.method, 'default')
        self.assertGreater(opts.iter_max, 0)


class TestResultAccessorMethods(unittest.TestCase):
    """Test result accessor methods."""

    def setUp(self):
        """Set up test fixtures."""
        self.sn = SimpleNetworkStruct(nstations=3, nclasses=2)
        self.solver = SolverFLD(self.sn)

        # Create dummy result
        M, K = 3, 2
        self.result = FLDResult(
            QN=np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]),
            UN=np.array([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]]),
            RN=np.array([[10.0, 20.0], [30.0, 40.0], [50.0, 60.0]]),
            TN=np.array([[0.01, 0.02], [0.03, 0.04], [0.05, 0.06]]),
            CN=np.array([[90.0, 120.0]]),
            XN=np.array([[0.03, 0.04]]),
            method='matrix'
        )
        self.solver.result = self.result

    def test_get_avg_qlen(self):
        """Test getAvgQLen returns correct values."""
        qlen = self.solver.getAvgQLen()

        expected = np.mean(self.result.QN, axis=1)
        np.testing.assert_array_almost_equal(qlen, expected)

    def test_get_avg_util(self):
        """Test getAvgUtil returns correct values."""
        util = self.solver.getAvgUtil()

        expected = np.mean(self.result.UN, axis=1)
        np.testing.assert_array_almost_equal(util, expected)

    def test_get_avg_respt(self):
        """Test getAvgRespT returns correct values."""
        rt = self.solver.getAvgRespT()

        expected = np.mean(self.result.RN, axis=1)
        np.testing.assert_array_almost_equal(rt, expected)

    def test_get_tput(self):
        """Test getTput returns correct values."""
        tput = self.solver.getTput()

        expected = np.mean(self.result.TN, axis=0)
        np.testing.assert_array_almost_equal(tput, expected)

    def test_get_avg_sys_resp_t(self):
        """Test getAvgSysRespT returns sum of response times."""
        sys_rt = self.solver.getAvgSysRespT()

        expected = np.sum(self.result.RN, axis=0)
        np.testing.assert_array_almost_equal(sys_rt, expected)

    def test_get_avg_sys_tput(self):
        """Test getAvgSysTput returns scalar."""
        sys_tput = self.solver.getAvgSysTput()

        expected = np.mean(self.result.XN)
        self.assertAlmostEqual(sys_tput, expected)

    def test_accessor_without_result_raises_error(self):
        """Test that accessors raise RuntimeError if result is None."""
        solver = SolverFLD(self.sn)
        solver.result = None

        with self.assertRaises(RuntimeError):
            solver.getAvgQLen()

        with self.assertRaises(RuntimeError):
            solver.getAvgUtil()


class TestMethodResolution(unittest.TestCase):
    """Test method name resolution."""

    def test_resolve_default_method(self):
        """Test that 'default' resolves to 'matrix'."""
        sn = SimpleNetworkStruct(nstations=2, nclasses=1)
        solver = SolverFLD(sn, method='default')

        resolved = solver._resolve_method()
        self.assertEqual(resolved, 'matrix')

    def test_resolve_alias_methods(self):
        """Test that method aliases resolve correctly."""
        sn = SimpleNetworkStruct(nstations=2, nclasses=1)

        # Test fluid.* aliases
        solver = SolverFLD(sn, method='fluid.matrix')
        self.assertEqual(solver._resolve_method(), 'matrix')

        solver = SolverFLD(sn, method='fluid.mfq')
        self.assertEqual(solver._resolve_method(), 'mfq')


if __name__ == '__main__':
    unittest.main(verbosity=2)
