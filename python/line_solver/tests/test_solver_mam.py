"""
Unit tests for SolverMAM data structures and static methods.

This file contains unit tests for:
- MAMResult data structure validation
- SolverMAM initialization and configuration
- Static methods (getFeatureSet, listValidMethods, defaultOptions)
- Closed network detection

Comprehensive end-to-end algorithm testing is performed in:
- test_solver_mam_integration.py (32 integration tests)
- test_solver_mam_matlab_parity.py (14 parity validation tests)
"""

import unittest
import numpy as np
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver.solvers.solver_mam import SolverMAM, SolverMAMOptions
from line_solver.solvers.solver_mam.algorithms import MAMResult
from line_solver.solvers.solver_mam.utils.network_adapter import check_closed_network


class SimpleNetworkStruct:
    """Minimal mock NetworkStruct for unit tests."""

    def __init__(self, nstations=2, nclasses=1, njobs=None):
        self.nstations = nstations
        self.nclasses = nclasses
        self.njobs = njobs if njobs is not None else []
        self.nodetype = np.zeros(nstations, dtype=int)  # NodeType.QUEUE=2

    def compileStruct(self):
        """Return self for testing."""
        return self


class TestMAMResult(unittest.TestCase):
    """Test MAMResult data structure."""

    def test_mamresult_creation(self):
        """Test MAMResult can be created with required fields."""
        QN = np.array([[1.0, 2.0], [3.0, 4.0]])
        UN = np.array([[0.5, 0.6], [0.7, 0.8]])
        RN = np.array([[1.0, 1.5], [2.0, 2.5]])
        TN = np.array([[0.8, 0.9]])

        result = MAMResult(
            QN=QN, UN=UN, RN=RN, TN=TN,
            totiter=10, method='dec.source', runtime=0.5
        )

        self.assertEqual(result.method, 'dec.source')
        self.assertEqual(result.totiter, 10)
        self.assertEqual(result.runtime, 0.5)
        np.testing.assert_array_equal(result.QN, QN)
        np.testing.assert_array_equal(result.UN, UN)


class TestNetworkAdapter(unittest.TestCase):
    """Test network parameter extraction utilities."""

    def test_check_closed_network(self):
        """Test closed network detection."""
        sn = SimpleNetworkStruct(nstations=2, nclasses=1)
        sn.njobs = [5]  # Closed: fixed jobs

        is_closed = check_closed_network(sn)

        self.assertTrue(is_closed)


class TestSolverMAMConfiguration(unittest.TestCase):
    """Test SolverMAM configuration and static methods."""

    def test_solver_initialization(self):
        """Test SolverMAM initialization."""
        sn = SimpleNetworkStruct(nstations=2, nclasses=1)

        solver = SolverMAM(sn, method='dec.source')

        self.assertIsNotNone(solver)
        self.assertEqual(solver.options.method, 'dec.source')

    def test_solver_with_options(self):
        """Test SolverMAM with custom options."""
        sn = SimpleNetworkStruct(nstations=2, nclasses=1)
        opts = SolverMAMOptions(method='mna_open', tol=1e-8, max_iter=200)

        solver = SolverMAM(sn, options=opts)

        self.assertEqual(solver.options.tol, 1e-8)
        self.assertEqual(solver.options.max_iter, 200)

    def test_list_valid_methods(self):
        """Test listValidMethods static method."""
        methods = SolverMAM.listValidMethods()

        self.assertIn('default', methods)
        self.assertIn('dec.source', methods)
        self.assertIn('mna_open', methods)
        self.assertIn('inap', methods)
        self.assertIn('fj', methods)

    def test_supports_method(self):
        """Test supports static method for valid methods."""
        sn = SimpleNetworkStruct(nstations=2, nclasses=1)

        # dec.source should support this network
        can_solve, reason = SolverMAM.supports(sn, 'dec.source')
        self.assertTrue(can_solve)

        # Unknown method should not be supported
        can_solve, reason = SolverMAM.supports(sn, 'nonexistent')
        self.assertFalse(can_solve)

    def test_get_feature_set(self):
        """Test getFeatureSet static method."""
        features = SolverMAM.getFeatureSet()

        self.assertIn('open_networks', features)
        self.assertIn('closed_networks', features)
        self.assertIn('scheduling', features)
        self.assertTrue(features['open_networks'])
        self.assertTrue(features['closed_networks'])

    def test_default_options(self):
        """Test defaultOptions static method."""
        opts = SolverMAM.defaultOptions()

        self.assertIsInstance(opts, SolverMAMOptions)
        self.assertEqual(opts.method, 'default')
        self.assertGreater(opts.max_iter, 0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
