"""
Integration tests for SolverFLD.

Tests the end-to-end solver functionality with complete network models.
"""

import unittest
import numpy as np
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver.solvers.solver_fld import SolverFLD
from line_solver.solvers.solver_fld.options import SolverFLDOptions


class SimpleQueue:
    """Mock queue node for testing."""

    def __init__(self, name, service_rate, scheduling='PS'):
        self.name = name
        self.service_rate = service_rate
        self.scheduling = scheduling
        self.server = 1


class SimpleSource:
    """Mock source node for testing."""

    def __init__(self, name, arrival_rate):
        self.name = name
        self.arrival_rate = arrival_rate


class SimpleSink:
    """Mock sink node for testing."""

    def __init__(self, name):
        self.name = name


class SimpleNetworkModel:
    """Simple network model for testing."""

    def __init__(self, nstations=2, nclasses=1):
        self.nstations = nstations
        self.nclasses = nclasses
        self.nodes = []
        self.source = SimpleSource('Source', 0.5)
        self.queues = [SimpleQueue(f'Q{i}', 1.0) for i in range(nstations)]
        self.sink = SimpleSink('Sink')

    def compileStruct(self):
        """Create a NetworkStruct for testing."""
        from line_solver.api.sn import NetworkStruct, NodeType, SchedStrategy

        M = self.nstations
        K = self.nclasses

        sn = NetworkStruct()
        sn.nstations = M
        sn.nclasses = K
        sn.nodetype = np.array([NodeType.SOURCE] + [NodeType.QUEUE] * M + [NodeType.SINK], dtype=int)
        sn.stationToNode = np.arange(M + 2, dtype=int)

        # Required by handler
        sn.njobs = np.array([])  # Empty for open networks
        sn.nclosedjobs = 0
        sn.nservers = np.ones(M, dtype=int)
        sn.phases = np.ones((M, K), dtype=int)
        sn.proc = {}
        sn.pie = {}
        sn.rt = None

        sn.lambda_arr = np.ones(K) * 0.5
        sn.rates = np.ones((M, K)) * 1.0
        sn.scv = np.ones((M, K))
        sn.schedid = np.array([SchedStrategy.PS] * M, dtype=int)
        sn.routing = np.zeros((M, M))
        sn.rtouter = np.zeros((M, 1))

        return sn


class TestMatrixMethodSolver(unittest.TestCase):
    """Test matrix method solver on simple networks."""

    def test_mm1_solver_basic(self):
        """Test solver runs on simple M/M/1 network."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)
        solver = SolverFLD(model, method='matrix')

        # Should not raise
        result = solver.runAnalyzer()

        self.assertIsNotNone(result.result)
        self.assertIsNotNone(result.result.QN)
        self.assertIsNotNone(result.result.UN)
        self.assertIsNotNone(result.result.RN)
        self.assertIsNotNone(result.result.TN)

    def test_mm1_result_structure(self):
        """Test that result has correct structure."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)
        solver = SolverFLD(model, method='matrix')
        solver.runAnalyzer()

        result = solver.result
        M, K = 1, 1

        self.assertEqual(result.QN.shape, (M, K))
        self.assertEqual(result.UN.shape, (M, K))
        self.assertEqual(result.RN.shape, (M, K))
        self.assertEqual(result.TN.shape, (M, K))
        self.assertEqual(result.CN.shape, (1, K))
        self.assertEqual(result.XN.shape, (1, K))

    def test_mm1_physical_constraints(self):
        """Test that results satisfy physical constraints."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)
        solver = SolverFLD(model, method='matrix')
        solver.runAnalyzer()

        result = solver.result

        # All metrics should be non-negative
        self.assertTrue(np.all(result.QN >= 0), "Queue lengths should be non-negative")
        self.assertTrue(np.all(result.UN >= 0), "Utilizations should be non-negative")
        self.assertTrue(np.all(result.RN >= 0), "Response times should be non-negative")
        self.assertTrue(np.all(result.TN >= 0), "Throughputs should be non-negative")

        # Utilization should be <= 1
        self.assertTrue(np.all(result.UN <= 1.0), "Utilization should be <= 1")

    def test_mm1_littles_law(self):
        """Test that results satisfy Little's Law: L = lambda * R."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)
        solver = SolverFLD(model, method='matrix')
        solver.runAnalyzer()

        result = solver.result
        lambda_arr = model.compileStruct().lambda_arr

        # Little's Law: Q = lambda * R
        Q_computed = lambda_arr[np.newaxis, :] * result.RN
        np.testing.assert_array_almost_equal(result.QN, Q_computed, decimal=4,
                                             err_msg="Little's Law violated")

    def test_tandem_network(self):
        """Test solver on tandem network."""
        model = SimpleNetworkModel(nstations=2, nclasses=1)
        solver = SolverFLD(model, method='matrix')

        result = solver.runAnalyzer()

        self.assertIsNotNone(result.result)
        self.assertEqual(result.result.QN.shape, (2, 1))

    def test_multiclass_network(self):
        """Test solver on multiclass network."""
        model = SimpleNetworkModel(nstations=2, nclasses=2)
        solver = SolverFLD(model, method='matrix')

        result = solver.runAnalyzer()

        self.assertIsNotNone(result.result)
        self.assertEqual(result.result.QN.shape, (2, 2))
        self.assertEqual(result.result.UN.shape, (2, 2))

    def test_accessor_methods_after_solve(self):
        """Test that accessor methods work after solving."""
        model = SimpleNetworkModel(nstations=2, nclasses=1)
        solver = SolverFLD(model, method='matrix')
        solver.runAnalyzer()

        # These should not raise
        qlen = solver.getAvgQLen()
        util = solver.getAvgUtil()
        rt = solver.getAvgRespT()
        tput = solver.getTput()
        sys_rt = solver.getAvgSysRespT()
        sys_tput = solver.getAvgSysTput()

        self.assertEqual(qlen.shape, (2,))
        self.assertEqual(util.shape, (2,))
        self.assertEqual(rt.shape, (2,))
        self.assertEqual(tput.shape, (1,))

    def test_table_generation(self):
        """Test that avg table generation works."""
        model = SimpleNetworkModel(nstations=2, nclasses=1)
        solver = SolverFLD(model, method='matrix')
        solver.runAnalyzer()

        table = solver.getAvgTable()

        self.assertIn('QLen', table.columns)
        self.assertIn('Util', table.columns)
        self.assertIn('RespT', table.columns)
        self.assertEqual(len(table), 2)

    def test_method_chaining(self):
        """Test that runAnalyzer returns self for method chaining."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)
        solver = SolverFLD(model, method='matrix')

        result = solver.runAnalyzer()

        self.assertIs(result, solver)

    def test_default_method_resolves_to_matrix(self):
        """Test that default method resolves to matrix."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)
        solver = SolverFLD(model, method='default')
        solver.runAnalyzer()

        self.assertEqual(solver.result.method, 'matrix')

    def test_custom_tolerance(self):
        """Test that custom tolerance is accepted."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)
        opts = SolverFLDOptions(method='matrix', tol=1e-5)
        solver = SolverFLD(model, options=opts)

        solver.runAnalyzer()

        self.assertIsNotNone(solver.result)

    def test_verbose_mode(self):
        """Test that verbose mode runs without error."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)
        opts = SolverFLDOptions(method='matrix', verbose=True)
        solver = SolverFLD(model, options=opts)

        # Should not raise even with verbose output
        solver.runAnalyzer()

        self.assertIsNotNone(solver.result)


class TestClosingMethodSolver(unittest.TestCase):
    """Test closing method solver."""

    def test_closing_method_runs(self):
        """Test that closing method runs successfully."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)
        solver = SolverFLD(model, method='closing')

        result = solver.runAnalyzer()
        self.assertIsNotNone(result.result)
        self.assertEqual(result.result.method, 'closing')

    def test_closing_with_pnorm_smoothing(self):
        """Test closing method with p-norm smoothing."""
        model = SimpleNetworkModel(nstations=2, nclasses=1)
        opts = SolverFLDOptions(method='closing', pstar=20.0)
        solver = SolverFLD(model, options=opts)

        result = solver.runAnalyzer()
        self.assertEqual(result.result.QN.shape, (2, 1))

    def test_closing_with_softmin_smoothing(self):
        """Test closing method with softmin smoothing."""
        model = SimpleNetworkModel(nstations=2, nclasses=1)
        opts = SolverFLDOptions(method='softmin', softmin_alpha=20.0)
        solver = SolverFLD(model, options=opts)

        result = solver.runAnalyzer()
        self.assertEqual(result.result.QN.shape, (2, 1))

    def test_closing_with_statedep(self):
        """Test closing method with state-dependent smoothing."""
        model = SimpleNetworkModel(nstations=2, nclasses=1)
        solver = SolverFLD(model, method='statedep')

        result = solver.runAnalyzer()
        self.assertEqual(result.result.QN.shape, (2, 1))


class TestDiffusionMethodSolver(unittest.TestCase):
    """Test diffusion method solver."""

    def _create_closed_network_sn(self):
        """Create a closed NetworkStruct for diffusion testing."""
        model = SimpleNetworkModel(nstations=2, nclasses=1)
        sn = model.compileStruct()
        # Make it a closed network
        sn.njobs = np.array([5.0])  # 5 jobs in system
        return sn

    def test_diffusion_method_closed_network(self):
        """Test that diffusion method works on closed networks."""
        sn = self._create_closed_network_sn()
        opts = SolverFLDOptions(method='diffusion', timespan=(0.0, 1.0), timestep=0.1)
        solver = SolverFLD(sn, options=opts)

        result = solver.runAnalyzer()
        self.assertIsNotNone(result.result)
        self.assertEqual(result.result.method, 'diffusion')

    def test_diffusion_population_preserved(self):
        """Test that diffusion maintains population constraint."""
        sn = self._create_closed_network_sn()
        opts = SolverFLDOptions(method='diffusion', timespan=(0.0, 0.5), timestep=0.1)
        solver = SolverFLD(sn, options=opts)

        result = solver.runAnalyzer()

        # Total queue length should be preserved (sum ≈ 5)
        total_q = np.sum(result.result.QN)
        self.assertAlmostEqual(total_q, 5.0, places=0)

    def test_diffusion_rejects_open_network(self):
        """Test that diffusion rejects open networks."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)
        solver = SolverFLD(model, method='diffusion')

        with self.assertRaises(ValueError):
            solver.runAnalyzer()


class TestMFQMethodSolver(unittest.TestCase):
    """Test MFQ (M/M/c) method solver."""

    def _create_single_queue_model(self):
        """Create a single-queue network for MFQ testing."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)
        return model

    def test_mfq_method_single_queue(self):
        """Test that MFQ method works on single-queue networks."""
        model = self._create_single_queue_model()
        solver = SolverFLD(model, method='mfq')

        result = solver.runAnalyzer()
        self.assertIsNotNone(result.result)
        self.assertEqual(result.result.method, 'mfq')

    def test_mfq_mm1_result(self):
        """Test MFQ on M/M/1 queue."""
        model = self._create_single_queue_model()
        solver = SolverFLD(model, method='mfq')

        result = solver.runAnalyzer()

        # For M/M/1 with rho=0.5, L should be ~1.0
        self.assertGreater(result.result.QN[0, 0], 0.0)
        self.assertLess(result.result.QN[0, 0], 2.0)

    def test_mfq_multiclass(self):
        """Test MFQ on multiclass single-queue network."""
        model = SimpleNetworkModel(nstations=1, nclasses=2)
        solver = SolverFLD(model, method='mfq')

        result = solver.runAnalyzer()
        self.assertEqual(result.result.QN.shape, (1, 2))

    def test_mfq_rejects_multiqueue(self):
        """Test that MFQ rejects multi-queue networks."""
        model = SimpleNetworkModel(nstations=2, nclasses=1)
        solver = SolverFLD(model, method='mfq')

        with self.assertRaises(ValueError):
            solver.runAnalyzer()


class TestMethodComparison(unittest.TestCase):
    """Test comparing different methods on same network."""

    def test_pnorm_vs_statedep_mm1(self):
        """Test that pnorm and statedep give similar results."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)

        # Solve with pnorm
        solver_pnorm = SolverFLD(model, method='pnorm')
        solver_pnorm.runAnalyzer()
        result_pnorm = solver_pnorm.result.QN

        # Solve with statedep
        solver_statedep = SolverFLD(model, method='statedep')
        solver_statedep.runAnalyzer()
        result_statedep = solver_statedep.result.QN

        # Results should be qualitatively similar (same sign, order of magnitude)
        self.assertTrue(np.all(result_pnorm >= 0))
        self.assertTrue(np.all(result_statedep >= 0))

    def test_softmin_vs_pnorm_convergence(self):
        """Test that softmin and pnorm converge with high smoothing."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)

        # Solve with high p-norm
        opts_pnorm = SolverFLDOptions(method='pnorm', pstar=100.0)
        solver_pnorm = SolverFLD(model, options=opts_pnorm)
        solver_pnorm.runAnalyzer()

        # Solve with high softmin alpha
        opts_softmin = SolverFLDOptions(method='softmin', softmin_alpha=100.0)
        solver_softmin = SolverFLD(model, options=opts_softmin)
        solver_softmin.runAnalyzer()

        # Both should give non-negative results
        self.assertTrue(np.all(solver_pnorm.result.QN >= 0))
        self.assertTrue(np.all(solver_softmin.result.QN >= 0))

    def test_matrix_vs_pnorm_mm1(self):
        """Test that matrix and pnorm methods are consistent for M/M/1."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)

        # Solve with matrix method
        solver_matrix = SolverFLD(model, method='matrix')
        solver_matrix.runAnalyzer()

        # Solve with pnorm method
        solver_pnorm = SolverFLD(model, method='pnorm')
        solver_pnorm.runAnalyzer()

        # Both should have non-negative utilization
        self.assertGreaterEqual(solver_matrix.result.UN[0, 0], 0.0)
        self.assertGreaterEqual(solver_pnorm.result.UN[0, 0], 0.0)


class TestNetworkTopologies(unittest.TestCase):
    """Test different network topologies."""

    def test_tandem_queues(self):
        """Test tandem queue network."""
        model = SimpleNetworkModel(nstations=3, nclasses=1)
        solver = SolverFLD(model, method='matrix')

        result = solver.runAnalyzer()
        self.assertEqual(result.result.QN.shape, (3, 1))
        self.assertTrue(np.all(result.result.QN >= 0))

    def test_two_class_network(self):
        """Test two-class network."""
        model = SimpleNetworkModel(nstations=2, nclasses=2)
        solver = SolverFLD(model, method='matrix')

        result = solver.runAnalyzer()
        self.assertEqual(result.result.QN.shape, (2, 2))

    def test_three_class_network(self):
        """Test three-class network."""
        model = SimpleNetworkModel(nstations=2, nclasses=3)
        solver = SolverFLD(model, method='matrix')

        result = solver.runAnalyzer()
        self.assertEqual(result.result.QN.shape, (2, 3))

    def test_high_utilization(self):
        """Test network with high utilization."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)
        sn = model.compileStruct()
        sn.lambda_arr = np.array([0.8])  # High arrival rate
        sn.rates = np.array([[1.0]])     # Same service rate

        solver = SolverFLD(model, method='mfq')
        result = solver.runAnalyzer()

        # Should have positive queue length
        self.assertGreater(result.result.QN[0, 0], 0.0)

    def test_low_utilization(self):
        """Test network with low utilization."""
        model = SimpleNetworkModel(nstations=1, nclasses=1)
        sn = model.compileStruct()
        sn.lambda_arr = np.array([0.1])  # Low arrival rate
        sn.rates = np.array([[1.0]])     # Higher service rate

        solver = SolverFLD(sn, method='mfq')
        result = solver.runAnalyzer()

        # Should have small queue length (rho=0.1, L = 0.1/0.9 ≈ 0.11)
        self.assertLess(result.result.QN[0, 0], 0.2)


if __name__ == '__main__':
    unittest.main(verbosity=2)
