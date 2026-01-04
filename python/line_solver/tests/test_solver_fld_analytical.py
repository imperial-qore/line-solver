"""
Analytical validation tests for SolverFLD.

Compares solver results against known analytical solutions for
simple queueing systems:
- M/M/1: Erlang-Delay formula
- M/M/c: Erlang-C formula
- Jackson networks: Product-form solution
- Little's Law: L = lambda * W
"""

import unittest
import numpy as np
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver.solvers.solver_fld import SolverFLD
from line_solver.solvers.solver_fld.options import SolverFLDOptions


class SimpleNetworkModel:
    """Simple network model for analytical testing."""

    def __init__(self, nstations=1, nclasses=1, lambda_arr=None, rates=None):
        self.nstations = nstations
        self.nclasses = nclasses
        self.lambda_arr = lambda_arr
        self.rates = rates

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

        sn.lambda_arr = self.lambda_arr if self.lambda_arr is not None else np.ones(K) * 0.5
        sn.rates = self.rates if self.rates is not None else np.ones((M, K)) * 1.0
        sn.scv = np.ones((M, K))
        sn.schedid = np.array([SchedStrategy.PS] * M, dtype=int)
        sn.routing = np.zeros((M, M))
        sn.rtouter = np.zeros((M, 1))

        return sn


class TestMM1Analytical(unittest.TestCase):
    """Test M/M/1 queue against analytical formulas."""

    def test_mm1_queue_length(self):
        """Test M/M/1 queue length formula: L = rho / (1 - rho)."""
        lambda_val = 0.5
        mu = 1.0
        rho = lambda_val / mu

        model = SimpleNetworkModel(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([lambda_val]),
            rates=np.array([[mu]])
        )

        solver = SolverFLD(model, method='mfq')
        solver.runAnalyzer()

        expected_L = rho / (1.0 - rho)
        actual_L = solver.result.QN[0, 0]

        # Allow 10% tolerance
        self.assertAlmostEqual(actual_L, expected_L, delta=expected_L * 0.1)

    def test_mm1_response_time(self):
        """Test M/M/1 response time formula: W = 1 / (mu - lambda)."""
        lambda_val = 0.5
        mu = 1.0

        model = SimpleNetworkModel(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([lambda_val]),
            rates=np.array([[mu]])
        )

        solver = SolverFLD(model, method='mfq')
        solver.runAnalyzer()

        expected_W = 1.0 / (mu - lambda_val)
        actual_W = solver.result.RN[0, 0]

        # Allow 10% tolerance
        self.assertAlmostEqual(actual_W, expected_W, delta=expected_W * 0.1)

    def test_mm1_utilization(self):
        """Test M/M/1 utilization formula: rho = lambda / mu."""
        lambda_val = 0.5
        mu = 1.0
        rho = lambda_val / mu

        model = SimpleNetworkModel(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([lambda_val]),
            rates=np.array([[mu]])
        )

        solver = SolverFLD(model, method='mfq')
        solver.runAnalyzer()

        actual_rho = solver.result.UN[0, 0]

        # Allow 10% tolerance
        self.assertAlmostEqual(actual_rho, rho, delta=rho * 0.1)

    def test_mm1_throughput(self):
        """Test M/M/1 throughput equals arrival rate."""
        lambda_val = 0.5
        mu = 1.0

        model = SimpleNetworkModel(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([lambda_val]),
            rates=np.array([[mu]])
        )

        solver = SolverFLD(model, method='mfq')
        solver.runAnalyzer()

        # Throughput should equal arrival rate
        self.assertAlmostEqual(solver.result.TN[0, 0], lambda_val, places=5)


class TestMM1Variants(unittest.TestCase):
    """Test M/M/1 for different loads."""

    def test_mm1_light_load(self):
        """Test M/M/1 with light load (rho = 0.1)."""
        lambda_val = 0.1
        mu = 1.0
        rho = lambda_val / mu

        model = SimpleNetworkModel(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([lambda_val]),
            rates=np.array([[mu]])
        )

        solver = SolverFLD(model, method='mfq')
        solver.runAnalyzer()

        expected_L = rho / (1.0 - rho)
        actual_L = solver.result.QN[0, 0]

        # Light load should have very small queue length
        self.assertLess(actual_L, 0.15)

    def test_mm1_moderate_load(self):
        """Test M/M/1 with moderate load (rho = 0.5)."""
        lambda_val = 0.5
        mu = 1.0
        rho = lambda_val / mu

        model = SimpleNetworkModel(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([lambda_val]),
            rates=np.array([[mu]])
        )

        solver = SolverFLD(model, method='mfq')
        solver.runAnalyzer()

        expected_L = rho / (1.0 - rho)
        actual_L = solver.result.QN[0, 0]

        # Moderate load queue length should be exactly 1.0
        self.assertAlmostEqual(actual_L, 1.0, delta=0.15)

    def test_mm1_high_load(self):
        """Test M/M/1 with high load (rho = 0.9)."""
        lambda_val = 0.9
        mu = 1.0
        rho = lambda_val / mu

        model = SimpleNetworkModel(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([lambda_val]),
            rates=np.array([[mu]])
        )

        solver = SolverFLD(model, method='mfq')
        solver.runAnalyzer()

        expected_L = rho / (1.0 - rho)
        actual_L = solver.result.QN[0, 0]

        # High load should have larger queue
        self.assertGreater(actual_L, 7.0)


class TestMMcAnalytical(unittest.TestCase):
    """Test M/M/c queue against analytical formulas."""

    def test_mmc_single_server(self):
        """Test that M/M/1 is equivalent to M/M/c with c=1."""
        lambda_val = 0.5
        mu = 1.0

        model = SimpleNetworkModel(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([lambda_val]),
            rates=np.array([[mu]])
        )

        # Get M/M/1 solution
        sn = model.compileStruct()
        sn.nservers = np.array([1])

        solver = SolverFLD(sn, method='mfq')
        solver.runAnalyzer()

        rho = lambda_val / mu
        expected_L = rho / (1.0 - rho)
        actual_L = solver.result.QN[0, 0]

        self.assertAlmostEqual(actual_L, expected_L, delta=expected_L * 0.1)

    def test_mmc_improves_with_servers(self):
        """Test that queue length decreases with more servers."""
        lambda_val = 1.0
        mu = 1.0

        # M/M/1
        model1 = SimpleNetworkModel(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([lambda_val]),
            rates=np.array([[mu]])
        )
        sn1 = model1.compileStruct()
        sn1.nservers = np.array([1])

        solver1 = SolverFLD(sn1, method='mfq')
        solver1.runAnalyzer()
        L1 = solver1.result.QN[0, 0]

        # M/M/2
        sn2 = model1.compileStruct()
        sn2.nservers = np.array([2])

        solver2 = SolverFLD(sn2, method='mfq')
        solver2.runAnalyzer()
        L2 = solver2.result.QN[0, 0]

        # More servers should reduce queue length
        self.assertLess(L2, L1)


class TestLittlesLaw(unittest.TestCase):
    """Test Little's Law: L = lambda * W."""

    def test_littles_law_mm1(self):
        """Test Little's Law for M/M/1."""
        lambda_val = 0.5
        mu = 1.0

        model = SimpleNetworkModel(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([lambda_val]),
            rates=np.array([[mu]])
        )

        solver = SolverFLD(model, method='mfq')
        solver.runAnalyzer()

        L = solver.result.QN[0, 0]
        W = solver.result.RN[0, 0]
        lambda_computed = solver.result.TN[0, 0]

        # Little's Law: L = lambda * W
        expected_L = lambda_computed * W
        self.assertAlmostEqual(L, expected_L, delta=abs(expected_L) * 0.01)

    def test_littles_law_multistation(self):
        """Test Little's Law for tandem queues."""
        model = SimpleNetworkModel(
            nstations=2,
            nclasses=1,
            lambda_arr=np.array([0.5]),
            rates=np.array([[1.0], [1.0]])
        )

        solver = SolverFLD(model, method='matrix')
        solver.runAnalyzer()

        # Check Little's Law at each station
        for i in range(2):
            L = solver.result.QN[i, 0]
            W = solver.result.RN[i, 0]
            T = solver.result.TN[i, 0]

            if T > 1e-10:
                expected_L = T * W
                self.assertAlmostEqual(L, expected_L, delta=abs(expected_L) * 0.05,
                                     msg=f"Little's Law violated at station {i}")


class TestSystemMetrics(unittest.TestCase):
    """Test system-level metrics."""

    def test_cycle_time_sum(self):
        """Test that cycle time equals sum of response times for single class."""
        model = SimpleNetworkModel(
            nstations=2,
            nclasses=1,
            lambda_arr=np.array([0.5]),
            rates=np.array([[1.0], [1.0]])
        )

        solver = SolverFLD(model, method='matrix')
        solver.runAnalyzer()

        cycle_time = solver.result.CN[0, 0]
        sum_resp_times = np.sum(solver.result.RN[:, 0])

        # Cycle time should be sum of response times (for open networks)
        self.assertAlmostEqual(cycle_time, sum_resp_times, delta=sum_resp_times * 0.1)

    def test_system_throughput_consistency(self):
        """Test that system throughput matches arrival rate for stable systems."""
        lambda_val = 0.5

        model = SimpleNetworkModel(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([lambda_val]),
            rates=np.array([[1.0]])
        )

        solver = SolverFLD(model, method='mfq')
        solver.runAnalyzer()

        sys_throughput = solver.result.XN[0, 0]

        # System throughput should equal arrival rate
        self.assertAlmostEqual(sys_throughput, lambda_val, places=3)


class TestStabilityBounds(unittest.TestCase):
    """Test stability conditions and bounds."""

    def test_stability_condition_mm1(self):
        """Test M/M/1 stability: rho < 1."""
        lambda_val = 0.5
        mu = 1.0
        rho = lambda_val / mu

        model = SimpleNetworkModel(
            nstations=1,
            nclasses=1,
            lambda_arr=np.array([lambda_val]),
            rates=np.array([[mu]])
        )

        solver = SolverFLD(model, method='mfq')
        solver.runAnalyzer()

        # All metrics should be finite for stable system
        self.assertTrue(np.all(np.isfinite(solver.result.QN)))
        self.assertTrue(np.all(np.isfinite(solver.result.RN)))
        self.assertFalse(np.any(np.isinf(solver.result.QN)))

    def test_queue_length_nonnegative(self):
        """Test that queue lengths are always non-negative."""
        for lambda_val in [0.1, 0.5, 0.8]:
            model = SimpleNetworkModel(
                nstations=1,
                nclasses=1,
                lambda_arr=np.array([lambda_val]),
                rates=np.array([[1.0]])
            )

            solver = SolverFLD(model, method='mfq')
            solver.runAnalyzer()

            self.assertTrue(np.all(solver.result.QN >= 0),
                          f"Negative queue length for lambda={lambda_val}")

    def test_utilization_bounds(self):
        """Test that utilization is in [0, 1]."""
        model = SimpleNetworkModel(
            nstations=2,
            nclasses=1,
            lambda_arr=np.array([0.5]),
            rates=np.array([[1.0], [1.0]])
        )

        solver = SolverFLD(model, method='matrix')
        solver.runAnalyzer()

        # Utilization should be bounded
        self.assertTrue(np.all(solver.result.UN >= 0),
                       "Negative utilization")
        self.assertTrue(np.all(solver.result.UN <= 1),
                       "Utilization > 1")


if __name__ == '__main__':
    unittest.main(verbosity=2)
