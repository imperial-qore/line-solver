"""Regression tests for FLD AoI parity and MFQ dispatch."""

import os
import sys
import unittest
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver.lib.thirdparty.aoi import solve_bufferless, solve_singlebuffer
from line_solver.api.sn import NetworkStruct, NodeType, SchedStrategy
from line_solver.solvers.solver_fld import SolverFLD


def _exp_dist(rate: float):
    return SimpleNamespace(mean=1.0 / rate, cv=1.0)


def build_aoi_struct(capacity: int) -> NetworkStruct:
    sn = NetworkStruct()
    sn.nstations = 2
    sn.nclasses = 1
    sn.nnodes = 0
    sn.nclosedjobs = 0
    sn.njobs = np.array([np.inf])
    sn.nodetype = np.array([NodeType.SOURCE, NodeType.QUEUE, NodeType.SINK], dtype=int)
    sn.nodeToStation = np.array([0, 1, -1], dtype=int)
    sn.stationToNode = np.array([0, 1], dtype=int)
    sn.nservers = np.array([1, 1], dtype=int)
    sn.cap = np.array([np.inf, capacity], dtype=float)
    sn.sched = {0: SchedStrategy.EXT, 1: SchedStrategy.FCFS}
    sn.schedid = np.array([SchedStrategy.EXT, SchedStrategy.FCFS], dtype=int)
    sn.lambda_arr = np.array([0.8], dtype=float)
    sn.rates = np.array([[0.8], [2.0]], dtype=float)
    sn.scv = np.ones((2, 1))
    sn.proc = {
        0: {0: _exp_dist(0.8)},
        1: {0: _exp_dist(2.0)},
    }
    sn.rt = np.zeros((2, 2))
    return sn


class TestAoISolvers(unittest.TestCase):
    def test_bufferless_exact_solver(self):
        result = solve_bufferless(
            tau=np.array([1.0]),
            T=np.array([[-0.8]]),
            sigma=np.array([1.0]),
            S=np.array([[-2.0]]),
            p=0.0,
        )

        self.assertEqual(result['status'], 'success')
        self.assertEqual(result['systemType'], 'bufferless')
        self.assertTrue(np.isfinite(result['AoI_mean']))
        self.assertTrue(np.isfinite(result['PAoI_mean']))
        self.assertGreater(result['AoI_mean'], 0.0)
        self.assertGreater(result['PAoI_mean'], 0.0)
        self.assertEqual(result['AoI_A'].shape[0], result['AoI_A'].shape[1])

    def test_singlebuffer_exact_solver(self):
        result = solve_singlebuffer(
            lambda_rate=0.8,
            sigma=np.array([1.0]),
            S=np.array([[-2.0]]),
            r=0.0,
        )

        self.assertEqual(result['status'], 'success')
        self.assertEqual(result['systemType'], 'singlebuffer')
        self.assertTrue(np.isfinite(result['AoI_mean']))
        self.assertTrue(np.isfinite(result['PAoI_mean']))
        self.assertGreater(result['AoI_mean'], 0.0)
        self.assertGreater(result['PAoI_mean'], 0.0)
        self.assertEqual(result['PAoI_A'].shape[0], result['PAoI_A'].shape[1])


class TestSolverFLDAoI(unittest.TestCase):
    def test_mfq_dispatch_exposes_aoi_accessors(self):
        for capacity, expected_type in ((1, 'bufferless'), (2, 'singlebuffer')):
            with self.subTest(capacity=capacity):
                solver = SolverFLD(build_aoi_struct(capacity), method='mfq')
                solver.runAnalyzer()

                self.assertEqual(solver.result.method, 'mfq')
                self.assertEqual(solver.result.aoiResults['systemType'], expected_type)

                aoi, paoi, table = solver.getAvgAoI()
                self.assertAlmostEqual(aoi['mean'], solver.result.aoiResults['AoI_mean'])
                self.assertAlmostEqual(paoi['mean'], solver.result.aoiResults['PAoI_mean'])
                self.assertEqual(list(table['Metric']), ['AoI', 'Peak AoI'])

                aoi_cdf, paoi_cdf = solver.getCdfAoI()
                self.assertEqual(aoi_cdf.shape[1], 2)
                self.assertEqual(paoi_cdf.shape[1], 2)
                self.assertTrue(np.all(np.diff(aoi_cdf[:, 1]) >= 0))
                self.assertTrue(np.all(np.diff(paoi_cdf[:, 1]) >= 0))
                self.assertTrue(np.all(np.diff(aoi_cdf[:, 0]) >= -1e-8))
                self.assertTrue(np.all(np.diff(paoi_cdf[:, 0]) >= -1e-8))


if __name__ == '__main__':
    unittest.main()
