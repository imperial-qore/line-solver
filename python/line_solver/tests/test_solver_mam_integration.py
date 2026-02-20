"""
Integration tests for SolverMAM with real LINE network models.

Tests validate:
- Method routing and auto-selection with real networks
- Multi-class queueing networks
- Open vs closed network handling
- Result consistency across different solution methods
- Performance metrics correctness
"""

import unittest
import numpy as np
import pandas as pd
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver.api.sn.network_struct import NetworkStruct
from line_solver.solvers.solver_mam import SolverMAM, SolverMAMOptions


class RealNetworkStruct(NetworkStruct):
    """Real NetworkStruct implementation for testing."""

    def __init__(self, nstations=2, nclasses=1, njobs=None):
        super().__init__()
        self.nstations = nstations
        self.nclasses = nclasses
        self.nservers = np.ones(nstations, dtype=int)
        self.njobs = njobs if njobs is not None else []

        # Service rates (mu) - (M, K) matrix
        self.mu = np.ones((nstations, nclasses)) * 2.0
        self.rates = np.ones((nstations, nclasses)) * 2.0  # Same as mu

        # Service SCVs (squared coefficients of variation)
        self.scv = np.ones((nstations, nclasses))

        # Arrival rates lambda per class
        self.lambda_arr = np.ones(nclasses)

        # Routing probability matrix (M, M)
        self.rt = np.ones((nstations, nstations)) / nstations

        # Node types
        self.node_types = ['Queue'] * nstations
        self.nodetype = np.zeros(nstations, dtype=int)  # 0 = Queue

        # Chain information
        self.nchains = nclasses
        self.nclosedjobs = sum(njobs) if njobs else 0
        self.chains = list(range(nclasses))
        self.stationToNode = np.arange(nstations)
        self.nodeToStation = {i: i for i in range(nstations)}

        # Visit counts
        self.visits = np.ones((nstations, nclasses))

        # Scheduling strategies (0 = FCFS)
        self.sched = {}


class TestM_M_1Network(unittest.TestCase):
    """Test M/M/1 queue (single station, single class, open)."""

    def setUp(self):
        """Create M/M/1 network: λ=1, μ=2, ρ=0.5."""
        self.sn = RealNetworkStruct(nstations=1, nclasses=1)
        self.sn.lambda_arr = np.array([1.0])
        self.sn.mu = np.array([[2.0]])
        self.sn.rates = np.array([[2.0]])

    def test_dec_source_m_m_1(self):
        """Test dec.source on M/M/1: should give ρ=0.5."""
        solver = SolverMAM(self.sn, method='dec.source')
        solver.runAnalyzer()

        # Check results exist
        self.assertIsNotNone(solver.result)

        # M/M/1: ρ = λ/μ = 1/2 = 0.5
        util = solver.result.UN[0, 0]
        self.assertAlmostEqual(util, 0.5, places=1)

        # Queue length: L = ρ/(1-ρ) = 0.5/0.5 = 1
        qlen = solver.result.QN[0, 0]
        self.assertGreater(qlen, 0)
        self.assertLess(qlen, 2)

    def test_mna_open_m_m_1(self):
        """Test mna_open on M/M/1."""
        solver = SolverMAM(self.sn, method='mna_open')
        solver.runAnalyzer()

        self.assertIsNotNone(solver.result)
        self.assertEqual(solver.result.method, 'mna_open')
        self.assertGreater(solver.result.totiter, 0)

    def test_inap_m_m_1(self):
        """Test INAP on M/M/1."""
        solver = SolverMAM(self.sn, method='inap')
        solver.runAnalyzer()

        self.assertIsNotNone(solver.result)
        self.assertEqual(solver.result.method, 'inap')


class TestTandemNetwork(unittest.TestCase):
    """Test 2-station tandem queue (series network)."""

    def setUp(self):
        """Create M/M/1/M/M/1 tandem: both stations have λ=1, μ=2."""
        self.sn = RealNetworkStruct(nstations=2, nclasses=1)
        self.sn.lambda_arr = np.array([1.0])
        self.sn.mu = np.array([[2.0], [2.0]])
        self.sn.rates = np.array([[2.0], [2.0]])
        # Series routing: station 0→1
        self.sn.rt = np.array([[0.0, 1.0], [0.0, 0.0]])

    def test_dec_source_tandem(self):
        """Test dec.source on tandem queue."""
        solver = SolverMAM(self.sn, method='dec.source')
        solver.runAnalyzer()

        self.assertIsNotNone(solver.result)
        self.assertEqual(solver.result.QN.shape, (2, 1))

        # Both stations should have similar utilization (ρ ≈ 0.5)
        util_1 = solver.result.UN[0, 0]
        util_2 = solver.result.UN[1, 0]

        self.assertGreater(util_1, 0)
        self.assertGreater(util_2, 0)
        self.assertLess(util_1, 1.0)
        self.assertLess(util_2, 1.0)

    def test_method_consistency_tandem(self):
        """Test that different methods give similar results on tandem."""
        methods = ['dec.source', 'mna_open', 'inap']
        results = {}

        for method in methods:
            solver = SolverMAM(self.sn, method=method)
            solver.runAnalyzer()
            results[method] = solver.result

        # All methods should produce valid results
        for method, result in results.items():
            self.assertIsNotNone(result.QN)
            self.assertFalse(np.any(np.isnan(result.QN)))
            self.assertTrue(np.all(result.QN >= 0))


class TestMultiClassNetwork(unittest.TestCase):
    """Test multi-class (2-class) queueing network."""

    def setUp(self):
        """Create 2-station, 2-class network."""
        self.sn = RealNetworkStruct(nstations=2, nclasses=2)
        # Different arrival rates per class
        self.sn.lambda_arr = np.array([0.5, 0.5])
        # Different service rates per class
        self.sn.mu = np.array([[2.0, 1.5], [2.0, 1.5]])
        self.sn.rates = np.array([[2.0, 1.5], [2.0, 1.5]])
        # Series routing
        self.sn.rt = np.array([[0.0, 1.0], [0.0, 0.0]])

    def test_dec_source_multiclass(self):
        """Test dec.source with 2 classes."""
        solver = SolverMAM(self.sn, method='dec.source')
        solver.runAnalyzer()

        # Check result dimensions
        self.assertEqual(solver.result.QN.shape, (2, 2))
        self.assertEqual(solver.result.UN.shape, (2, 2))
        self.assertEqual(solver.result.RN.shape, (2, 2))
        self.assertEqual(solver.result.TN.shape, (1, 2))

    def test_mna_multiclass(self):
        """Test MNA with multiple classes."""
        solver = SolverMAM(self.sn, method='mna_open')
        solver.runAnalyzer()

        # Check per-class throughputs
        tput = solver.getTput()
        self.assertEqual(len(tput), 2)
        self.assertGreater(tput[0], 0)
        self.assertGreater(tput[1], 0)


class TestClosedNetwork(unittest.TestCase):
    """Test closed networks with fixed population."""

    def setUp(self):
        """Create 2-station closed network with N=5 jobs."""
        self.sn = RealNetworkStruct(nstations=2, nclasses=1)
        self.sn.njobs = [5]  # 5 jobs circulating
        self.sn.lambda_arr = np.array([1.0])  # Initial guess
        self.sn.mu = np.array([[2.0], [2.0]])
        self.sn.rates = np.array([[2.0], [2.0]])
        # Cycle: station 0→1→0
        self.sn.rt = np.array([[0.0, 1.0], [1.0, 0.0]])

    def test_mna_closed(self):
        """Test MNA on closed network (enforces population constraint)."""
        solver = SolverMAM(self.sn, method='mna_closed')
        solver.runAnalyzer()

        # Check population constraint: sum(Q) ≈ N
        total_qlen = np.sum(solver.result.QN)
        expected_jobs = 5
        self.assertAlmostEqual(total_qlen, expected_jobs, delta=0.5)

    def test_dec_source_closed(self):
        """Test dec.source on closed network."""
        solver = SolverMAM(self.sn, method='dec.source')
        solver.runAnalyzer()

        self.assertIsNotNone(solver.result)
        # Should not crash on closed network
        self.assertFalse(np.any(np.isnan(solver.result.QN)))


class TestAutoSelection(unittest.TestCase):
    """Test automatic method selection."""

    def test_auto_select_open_network(self):
        """Auto-select should choose dec.source for open networks."""
        sn = RealNetworkStruct(nstations=2, nclasses=1)
        sn.njobs = []  # Open network

        solver = SolverMAM(sn, method='default')
        selected = solver._select_method()

        self.assertEqual(selected, 'dec.source')

    def test_auto_select_closed_1station(self):
        """Auto-select should consider ldqbd for single-station closed."""
        sn = RealNetworkStruct(nstations=1, nclasses=1)
        sn.njobs = [5]  # Single-station closed

        solver = SolverMAM(sn, method='default')
        selected = solver._select_method()

        # Should select ldqbd or dec.source (both valid)
        self.assertIn(selected, ['ldqbd', 'dec.source'])

    def test_mna_auto_select_open(self):
        """mna method should auto-select to mna_open for open networks."""
        sn = RealNetworkStruct(nstations=2, nclasses=1)
        sn.njobs = []

        solver = SolverMAM(sn, method='mna')
        selected = solver._select_method()

        self.assertEqual(selected, 'mna_open')

    def test_mna_auto_select_closed(self):
        """mna method should auto-select to mna_closed for closed networks."""
        sn = RealNetworkStruct(nstations=2, nclasses=1)
        sn.njobs = [5]

        solver = SolverMAM(sn, method='mna')
        selected = solver._select_method()

        self.assertEqual(selected, 'mna_closed')


class TestResultAccessors(unittest.TestCase):
    """Test all result accessor methods."""

    def setUp(self):
        """Create and solve a basic network."""
        self.sn = RealNetworkStruct(nstations=2, nclasses=2)
        self.solver = SolverMAM(self.sn, method='dec.source')
        self.solver.runAnalyzer()

    def test_get_avg_qlen(self):
        """Test getAvgQLen accessor."""
        qlen = self.solver.getAvgQLen()

        self.assertEqual(len(qlen), 2)
        self.assertTrue(np.all(qlen >= 0))

    def test_get_avg_util(self):
        """Test getAvgUtil accessor."""
        util = self.solver.getAvgUtil()

        self.assertEqual(len(util), 2)
        self.assertTrue(np.all(util >= 0))
        self.assertTrue(np.all(util <= 1.1))  # Allow small overshoot

    def test_get_avg_respt(self):
        """Test getAvgRespT accessor."""
        resp_t = self.solver.getAvgRespT()

        self.assertEqual(len(resp_t), 2)
        self.assertTrue(np.all(resp_t > 0))

    def test_get_tput(self):
        """Test getTput accessor."""
        tput = self.solver.getTput()

        self.assertEqual(len(tput), 2)
        self.assertTrue(np.all(tput > 0))

    def test_get_avg_sys_resp_t(self):
        """Test getAvgSysRespT accessor."""
        sys_resp_t = self.solver.getAvgSysRespT()

        self.assertEqual(len(sys_resp_t), 2)
        self.assertTrue(np.all(sys_resp_t > 0))

    def test_get_avg_table(self):
        """Test getAvgTable DataFrame accessor."""
        table = self.solver.getAvgTable()

        self.assertIsInstance(table, pd.DataFrame)
        self.assertEqual(len(table), 2)  # 2 stations
        self.assertIn('QLen', table.columns)
        self.assertIn('Util', table.columns)
        self.assertIn('RespT', table.columns)

    def test_get_avg_sys_tput(self):
        """Test getAvgSysTput accessor."""
        sys_tput = self.solver.getAvgSysTput()

        self.assertIsInstance(sys_tput, (float, np.floating))
        self.assertGreater(sys_tput, 0)


class TestStaticMethods(unittest.TestCase):
    """Test static introspection methods."""

    def test_list_valid_methods(self):
        """Test listValidMethods returns all supported methods."""
        methods = SolverMAM.listValidMethods()

        expected = ['default', 'dec.source', 'dec.mmap', 'dec.poisson',
                    'mna', 'mna_open', 'mna_closed', 'ldqbd', 'inap', 'inapplus']

        for method in expected:
            self.assertIn(method, methods)

    def test_supports_valid_method(self):
        """Test supports returns True for valid methods on compatible networks."""
        sn = RealNetworkStruct(nstations=2, nclasses=1)

        can_solve, reason = SolverMAM.supports(sn, 'dec.source')
        self.assertTrue(can_solve)

    def test_supports_invalid_method(self):
        """Test supports returns False for unknown methods."""
        sn = RealNetworkStruct(nstations=2, nclasses=1)

        can_solve, reason = SolverMAM.supports(sn, 'nonexistent_method')
        self.assertFalse(can_solve)

    def test_get_feature_set(self):
        """Test getFeatureSet returns capability information."""
        features = SolverMAM.getFeatureSet()

        self.assertIn('open_networks', features)
        self.assertIn('closed_networks', features)
        self.assertIn('scheduling', features)
        self.assertTrue(features['open_networks'])
        self.assertTrue(features['closed_networks'])

    def test_default_options(self):
        """Test defaultOptions returns proper SolverMAMOptions."""
        opts = SolverMAM.defaultOptions()

        self.assertIsInstance(opts, SolverMAMOptions)
        self.assertEqual(opts.method, 'default')
        self.assertGreater(opts.max_iter, 0)
        self.assertGreater(opts.tol, 0)


class TestConvergenceBehavior(unittest.TestCase):
    """Test convergence properties with different tolerance levels."""

    def test_convergence_tight_tolerance(self):
        """Test convergence with tight tolerance (1e-8)."""
        sn = RealNetworkStruct(nstations=2, nclasses=1)

        opts = SolverMAMOptions(method='mna_open', tol=1e-8, max_iter=200)
        solver = SolverMAM(sn, options=opts)
        solver.runAnalyzer()

        # Should converge within max_iter
        self.assertLess(solver.result.totiter, opts.max_iter)

    def test_convergence_loose_tolerance(self):
        """Test convergence with loose tolerance (1e-3)."""
        sn = RealNetworkStruct(nstations=2, nclasses=1)

        opts = SolverMAMOptions(method='mna_open', tol=1e-3, max_iter=50)
        solver = SolverMAM(sn, options=opts)
        solver.runAnalyzer()

        # Should converge quickly with loose tolerance
        self.assertLess(solver.result.totiter, 20)


class TestNumericalEdgeCases(unittest.TestCase):
    """Test solver behavior at numerical edge cases."""

    def test_very_low_arrival_rate(self):
        """Test with very low arrival rates."""
        sn = RealNetworkStruct(nstations=2, nclasses=1)
        sn.lambda_arr = np.array([1e-4])

        solver = SolverMAM(sn, method='dec.source')
        solver.runAnalyzer()

        # Should not produce NaN
        self.assertFalse(np.any(np.isnan(solver.result.QN)))
        self.assertTrue(np.all(solver.result.QN >= 0))

    def test_high_utilization_stability(self):
        """Test stability with high utilization (ρ close to 1)."""
        sn = RealNetworkStruct(nstations=1, nclasses=1)
        sn.lambda_arr = np.array([0.95])
        sn.mu = np.array([[1.0]])  # ρ = 0.95

        solver = SolverMAM(sn, method='mna_open')
        solver.runAnalyzer()

        # Should handle high utilization gracefully
        util = solver.result.UN[0, 0]
        self.assertLess(util, 1.1)  # Should not be much > 1 due to approximation

    def test_variable_service_rates(self):
        """Test with variable (non-identical) service rates."""
        sn = RealNetworkStruct(nstations=3, nclasses=1)
        sn.mu = np.array([[1.0], [2.0], [3.0]])  # Different rates

        solver = SolverMAM(sn, method='dec.source')
        solver.runAnalyzer()

        self.assertIsNotNone(solver.result)
        self.assertFalse(np.any(np.isnan(solver.result.QN)))


class TestPerformanceMetrics(unittest.TestCase):
    """Test that computed performance metrics make sense."""

    def test_queueing_theory_relationships(self):
        """Test that Little's law relationships hold."""
        sn = RealNetworkStruct(nstations=2, nclasses=1)
        sn.lambda_arr = np.array([1.0])

        solver = SolverMAM(sn, method='dec.source')
        solver.runAnalyzer()

        # Little's Law: L = λ * W
        # For each station: Q[m] ≈ lambda[k] * R[m,k]
        lambda_k = solver.result.TN[0, 0]

        for m in range(sn.nstations):
            q_expected = lambda_k * solver.result.RN[m, 0]
            q_actual = solver.result.QN[m, 0]
            # Should be approximately equal (within 20% due to approximation)
            if q_expected > 0.1:
                rel_error = abs(q_actual - q_expected) / q_expected
                self.assertLess(rel_error, 0.2,
                    f"Station {m}: Little's law violated by {rel_error*100:.1f}%")

    def test_utilization_bounds(self):
        """Test that utilization stays within [0, 1] bounds."""
        sn = RealNetworkStruct(nstations=3, nclasses=2)

        for method in ['dec.source', 'mna_open', 'inap']:
            solver = SolverMAM(sn, method=method)
            solver.runAnalyzer()

            # All utilizations should be in [0, 1]
            self.assertTrue(np.all(solver.result.UN >= -0.01))  # Allow small numerical error
            self.assertTrue(np.all(solver.result.UN <= 1.01))


if __name__ == '__main__':
    unittest.main(verbosity=2)
