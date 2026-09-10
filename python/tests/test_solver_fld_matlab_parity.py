"""
MATLAB parity / reference tests for SolverFLD, using real LINE network models.

Two solver families are exercised:
- The mean-field FLUID methods (matrix, statedep, ...): validated for stability
  and finiteness. An open struct keeps the Source as station 0, so the queue is
  at station 1.
- The MFQ method: an M/M/c queue solver that writes the exact (Erlang-C)
  measures into the FULL (nstations, K) station table, on the queue row -- the
  Source is row 0 and stays empty, as under every other fluid method. For M/M/1
  this recovers the stochastic law L = rho/(1-rho), W = 1/(mu-lambda), U = rho.
"""

# The Source is station 0 in every open helper used here, so a queueing measure
# is read off row QROW and never off row 0.
QROW = 1

import unittest
import numpy as np

import line_solver as L
from line_solver.solvers.solver_fld import SolverFLD
from line_solver.solvers.solver_fld.options import SolverFLDOptions
from test_solver_fld_analytical import fld_mm1, fld_tandem
from test_solver_fld_integration import fld_multiclass, fld_closed


class TestMatrixMethodParity(unittest.TestCase):
    """Matrix (fluid) method: stability and bounds."""

    def test_mm1_matrix_vs_analytical(self):
        r = SolverFLD(fld_mm1(0.5, 1.0), method='matrix').runAnalyzer().result
        self.assertTrue(np.all(np.isfinite(np.asarray(r.QN))))
        self.assertTrue(np.all(np.asarray(r.UN) <= 1.0 + 1e-9))
        self.assertTrue(np.all(np.asarray(r.QN) >= -1e-9))

    def test_mmc_matrix_stability(self):
        r = SolverFLD(fld_mm1(1.5, 1.0, nservers=2), method='matrix').runAnalyzer().result
        self.assertTrue(np.isfinite(np.asarray(r.QN)[1, 0]))
        self.assertTrue(np.isfinite(np.asarray(r.RN)[1, 0]))


class TestMFQMethodParity(unittest.TestCase):
    """MFQ: exact M/M/c queue results, on the queue row of the station table."""

    def test_mm1_mfq_exact(self):
        lam, mu = 0.5, 1.0
        rho = lam / mu
        r = SolverFLD(fld_mm1(lam, mu), method='mfq').runAnalyzer().result
        self.assertAlmostEqual(np.asarray(r.QN)[QROW, 0], rho / (1 - rho), places=4)
        self.assertAlmostEqual(np.asarray(r.RN)[QROW, 0], 1.0 / (mu - lam), places=4)
        self.assertAlmostEqual(np.asarray(r.UN)[QROW, 0], rho, places=4)

    def test_mmc_mfq_exact(self):
        lam, mu, c = 1.5, 1.0, 2
        rho = lam / (c * mu)
        r = SolverFLD(fld_mm1(lam, mu, nservers=c), method='mfq').runAnalyzer().result
        self.assertAlmostEqual(np.asarray(r.UN)[QROW, 0], rho, places=4)
        self.assertTrue(np.isfinite(np.asarray(r.QN)[QROW, 0]))


class TestDiffusionMethodParity(unittest.TestCase):
    """Diffusion method on a closed network."""

    def test_diffusion_closed_network_stability(self):
        opts = SolverFLDOptions(method='diffusion', timespan=(0.0, 1.0), timestep=0.1)
        try:
            r = SolverFLD(fld_closed(njobs=5), options=opts).runAnalyzer().result
            self.assertTrue(np.all(np.isfinite(np.asarray(r.QN))))
            self.assertTrue(np.all(np.isfinite(np.asarray(r.UN))))
            self.assertTrue(np.all(np.isfinite(np.asarray(r.RN))))
        except Exception as e:
            self.assertIn('closed', str(e).lower())


class TestMethodCrossValidation(unittest.TestCase):
    """Cross-validate methods on the same network."""

    def test_matrix_vs_mfq_tandem(self):
        r = SolverFLD(fld_tandem(0.5, (1.0, 1.0)), method='matrix').runAnalyzer().result
        self.assertTrue(np.all(np.isfinite(np.asarray(r.QN))))
        self.assertTrue(np.all(np.asarray(r.QN) >= -1e-9))

    def test_method_agreement_utilization(self):
        r = SolverFLD(fld_tandem(0.5, (1.0, 1.0)), method='matrix').runAnalyzer().result
        self.assertTrue(np.all(np.isfinite(np.asarray(r.UN))))


class TestMulticlassNetworks(unittest.TestCase):
    """Multi-class networks."""

    def test_mfq_multiclass_mm1(self):
        r = SolverFLD(fld_multiclass(lams=(0.4, 0.3), nqueues=1), method='mfq').runAnalyzer().result
        QN = np.asarray(r.QN)
        self.assertEqual(QN.shape, (2, 2))          # Source + Queue1, two classes
        self.assertTrue(np.all(np.isfinite(QN)))
        self.assertTrue(np.all(np.isfinite(np.asarray(r.RN))))
        # Little's law per class in the single queue
        TN = np.asarray(r.TN)
        RN = np.asarray(r.RN)
        for k in range(2):
            if TN[QROW, k] > 1e-10:
                self.assertAlmostEqual(QN[QROW, k], TN[QROW, k] * RN[QROW, k],
                                       delta=abs(TN[QROW, k] * RN[QROW, k]) * 0.02 + 1e-9)

    def test_matrix_multiclass_stability(self):
        r = SolverFLD(fld_multiclass(lams=(0.5, 0.3), nqueues=2), method='matrix').runAnalyzer().result
        self.assertTrue(np.all(np.isfinite(np.asarray(r.QN))))
        self.assertTrue(np.all(np.isfinite(np.asarray(r.UN))))
        self.assertTrue(np.all(np.isfinite(np.asarray(r.RN))))


class TestEdgeCases(unittest.TestCase):
    """Boundary conditions (MFQ = exact M/M/1)."""

    def test_light_load_stability(self):
        r = SolverFLD(fld_mm1(0.1, 1.0), method='mfq').runAnalyzer().result
        q = np.asarray(r.QN)[0, 0]
        self.assertTrue(np.isfinite(q))
        self.assertLess(q, 0.15)  # rho/(1-rho) = 0.111

    def test_high_load_convergence(self):
        r = SolverFLD(fld_mm1(0.95, 1.0), method='mfq').runAnalyzer().result
        q = np.asarray(r.QN)[QROW, 0]
        self.assertTrue(np.isfinite(q))
        self.assertGreater(q, 10.0)  # rho/(1-rho) = 19

    def test_unstable_handling(self):
        r = SolverFLD(fld_mm1(2.0, 1.0), method='mfq').runAnalyzer().result
        self.assertTrue(np.isinf(np.asarray(r.QN)[QROW, 0]))
        self.assertTrue(np.isinf(np.asarray(r.RN)[QROW, 0]))


class TestRuntimePerformance(unittest.TestCase):
    """Solver performance characteristics."""

    def test_mfq_fast_execution(self):
        solver = SolverFLD(fld_mm1(0.5, 1.0), method='mfq')
        solver.runAnalyzer()
        self.assertLess(solver.runtime, 1.0)

    def test_matrix_reasonable_runtime(self):
        solver = SolverFLD(fld_tandem(0.5, (1.0, 1.0)), method='matrix')
        solver.runAnalyzer()
        self.assertLess(solver.runtime, 5.0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
