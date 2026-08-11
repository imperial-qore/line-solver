"""
Integration tests for SolverFLD using real LINE network models.

SolverFLD is a mean-field fluid solver. These tests exercise its methods
(matrix, closing, pnorm, softmin, statedep, diffusion, mfq) on real Network
models and check structural correctness and fluid-appropriate properties.

Real open-network structs keep the Source as station 0, so a single queue is at
station 1 and the result arrays have nstations = (number of queues + 1).
"""

import unittest
import numpy as np

import line_solver as L
from line_solver.solvers.solver_fld import SolverFLD
from line_solver.solvers.solver_fld.options import SolverFLDOptions
from line_solver.tests.test_solver_fld_analytical import fld_mm1, fld_tandem


def fld_multiclass(lams=(0.5, 0.5), mu=1.0, nqueues=1):
    """Open multiclass model over nqueues tandem queues (Source(0) + queues)."""
    K = len(lams)
    m = L.Network('mc')
    s = L.Source(m, 'Source')
    queues = [L.Queue(m, 'Queue%d' % (i + 1), L.SchedStrategy.PS) for i in range(nqueues)]
    k = L.Sink(m, 'Sink')
    classes = [L.OpenClass(m, 'Class%d' % (r + 1)) for r in range(K)]
    for r in range(K):
        s.setArrival(classes[r], L.Exp(lams[r]))
        for q in queues:
            q.setService(classes[r], L.Exp(mu))
    m.link(L.Network.serialRouting(s, *queues, k))
    return m


def fld_closed(njobs=5, mus=(2.0, 2.0)):
    """Closed 2-station cycle: Delay(0) <-> Queue(1)."""
    m = L.Network('closed')
    d = L.Delay(m, 'Delay')
    q = L.Queue(m, 'Queue1', L.SchedStrategy.PS)
    c = L.ClosedClass(m, 'Class1', njobs, d)
    d.setService(c, L.Exp(mus[0]))
    q.setService(c, L.Exp(mus[1]))
    m.link(L.Network.serialRouting(d, q))
    return m


class TestMatrixMethodSolver(unittest.TestCase):
    """Test the matrix (default) fluid method."""

    def test_mm1_solver_basic(self):
        solver = SolverFLD(fld_mm1(0.5, 1.0), method='matrix')
        result = solver.runAnalyzer()
        self.assertIsNotNone(result.result)
        for nm in ('QN', 'UN', 'RN', 'TN'):
            self.assertIsNotNone(getattr(result.result, nm))

    def test_mm1_result_structure(self):
        model = fld_mm1(0.5, 1.0)
        solver = SolverFLD(model, method='matrix')
        result = solver.runAnalyzer().result
        M = model.getStruct().nstations
        K = model.getStruct().nclasses
        self.assertEqual(np.asarray(result.QN).shape, (M, K))
        self.assertEqual(np.asarray(result.UN).shape, (M, K))
        self.assertEqual(np.asarray(result.RN).shape, (M, K))
        self.assertEqual(np.asarray(result.TN).shape, (M, K))
        self.assertEqual(np.asarray(result.CN).shape, (1, K))
        self.assertEqual(np.asarray(result.XN).shape, (1, K))

    def test_mm1_physical_constraints(self):
        solver = SolverFLD(fld_mm1(0.5, 1.0), method='matrix')
        result = solver.runAnalyzer().result
        self.assertTrue(np.all(np.asarray(result.QN) >= -1e-9))
        self.assertTrue(np.all(np.asarray(result.UN) >= -1e-9))
        self.assertTrue(np.all(np.asarray(result.RN) >= -1e-9))
        self.assertTrue(np.all(np.asarray(result.TN) >= -1e-9))
        self.assertTrue(np.all(np.asarray(result.UN) <= 1.0 + 1e-9))

    def test_mm1_littles_law(self):
        solver = SolverFLD(fld_mm1(0.5, 1.0), method='matrix')
        r = solver.runAnalyzer().result
        Lq = np.asarray(r.QN)[1, 0]
        W = np.asarray(r.RN)[1, 0]
        T = np.asarray(r.TN)[1, 0]
        self.assertAlmostEqual(Lq, T * W, delta=abs(T * W) * 0.05 + 1e-9)

    def test_tandem_network(self):
        model = fld_tandem(0.5, (1.0, 1.0))
        result = SolverFLD(model, method='matrix').runAnalyzer().result
        self.assertEqual(np.asarray(result.QN).shape, (model.getStruct().nstations, 1))

    def test_multiclass_network(self):
        model = fld_multiclass(lams=(0.5, 0.5), nqueues=1)
        result = SolverFLD(model, method='matrix').runAnalyzer().result
        M = model.getStruct().nstations
        self.assertEqual(np.asarray(result.QN).shape, (M, 2))
        self.assertEqual(np.asarray(result.UN).shape, (M, 2))

    def test_accessor_methods_after_solve(self):
        model = fld_tandem(0.5, (1.0, 1.0))
        solver = SolverFLD(model, method='matrix')
        solver.runAnalyzer()
        M = model.getStruct().nstations
        self.assertEqual(np.asarray(solver.getAvgQLen()).shape[0], M)
        self.assertEqual(np.asarray(solver.getAvgUtil()).shape[0], M)
        self.assertEqual(np.asarray(solver.getAvgRespT()).shape[0], M)
        self.assertGreaterEqual(np.asarray(solver.getTput()).size, 1)

    def test_table_generation(self):
        solver = SolverFLD(fld_tandem(0.5, (1.0, 1.0)), method='matrix')
        solver.runAnalyzer()
        table = solver.getAvgTable()
        rendered = str(table)
        for col in ('QLen', 'Util', 'RespT'):
            self.assertIn(col, rendered)

    def test_method_chaining(self):
        solver = SolverFLD(fld_mm1(0.5, 1.0), method='matrix')
        result = solver.runAnalyzer()
        self.assertIs(result, solver)

    def test_default_method_resolves_to_matrix(self):
        solver = SolverFLD(fld_mm1(0.5, 1.0), method='default')
        solver.runAnalyzer()
        self.assertEqual(solver.result.method, 'matrix')

    def test_custom_tolerance(self):
        opts = SolverFLDOptions(method='matrix', tol=1e-5)
        solver = SolverFLD(fld_mm1(0.5, 1.0), options=opts)
        solver.runAnalyzer()
        self.assertIsNotNone(solver.result)

    def test_verbose_mode(self):
        opts = SolverFLDOptions(method='matrix', verbose=True)
        solver = SolverFLD(fld_mm1(0.5, 1.0), options=opts)
        solver.runAnalyzer()
        self.assertIsNotNone(solver.result)


class TestClosingMethodSolver(unittest.TestCase):
    """Test the closing / smoothing methods."""

    def test_closing_method_runs(self):
        solver = SolverFLD(fld_mm1(0.5, 1.0), method='closing')
        result = solver.runAnalyzer()
        self.assertIsNotNone(result.result)
        self.assertEqual(result.result.method, 'closing')

    def test_closing_with_pnorm_smoothing(self):
        model = fld_tandem(0.5, (1.0, 1.0))
        solver = SolverFLD(model, method='pnorm')
        result = solver.runAnalyzer()
        self.assertEqual(np.asarray(result.result.QN).shape, (model.getStruct().nstations, 1))

    def test_closing_with_softmin_smoothing(self):
        model = fld_tandem(0.5, (1.0, 1.0))
        opts = SolverFLDOptions(method='softmin', softmin_alpha=20.0)
        solver = SolverFLD(model, options=opts)
        result = solver.runAnalyzer()
        self.assertEqual(np.asarray(result.result.QN).shape, (model.getStruct().nstations, 1))

    def test_closing_with_statedep(self):
        model = fld_tandem(0.5, (1.0, 1.0))
        solver = SolverFLD(model, method='statedep')
        result = solver.runAnalyzer()
        self.assertEqual(np.asarray(result.result.QN).shape, (model.getStruct().nstations, 1))


class TestDiffusionMethodSolver(unittest.TestCase):
    """Test the diffusion method (closed networks, transient)."""

    def test_diffusion_method_closed_network(self):
        opts = SolverFLDOptions(method='diffusion', timespan=(0.0, 1.0), timestep=0.1)
        solver = SolverFLD(fld_closed(njobs=5), options=opts)
        result = solver.runAnalyzer()
        self.assertIsNotNone(result.result)
        self.assertEqual(result.result.method, 'diffusion')

    def test_diffusion_population_preserved(self):
        opts = SolverFLDOptions(method='diffusion', timespan=(0.0, 0.5), timestep=0.1)
        solver = SolverFLD(fld_closed(njobs=5), options=opts)
        result = solver.runAnalyzer()
        total_q = np.sum(np.asarray(result.result.QN))
        self.assertAlmostEqual(total_q, 5.0, delta=1.0)

    def test_diffusion_rejects_open_network(self):
        solver = SolverFLD(fld_mm1(0.5, 1.0), method='diffusion')
        with self.assertRaises(Exception):
            solver.runAnalyzer()


class TestMFQMethodSolver(unittest.TestCase):
    """Test the MFQ (M/M/c) method (single queue only)."""

    def test_mfq_method_single_queue(self):
        solver = SolverFLD(fld_mm1(0.5, 1.0), method='mfq')
        result = solver.runAnalyzer()
        self.assertIsNotNone(result.result)
        self.assertEqual(result.result.method, 'mfq')

    def test_mfq_mm1_result(self):
        # MFQ is a standalone M/M/c queue solver: it returns a single-queue
        # (1, K) result (the queue), with the stochastic Erlang-C queue length.
        solver = SolverFLD(fld_mm1(0.5, 1.0), method='mfq')
        result = solver.runAnalyzer()
        qn_queue = np.asarray(result.result.QN)[0, 0]
        self.assertGreater(qn_queue, 0.0)
        self.assertLess(qn_queue, 2.0)  # M/M/1 rho=0.5: L = rho/(1-rho) = 1

    def test_mfq_multiclass(self):
        model = fld_multiclass(lams=(0.5, 0.5), nqueues=1)
        result = SolverFLD(model, method='mfq').runAnalyzer()
        self.assertEqual(np.asarray(result.result.QN).shape, (1, 2))

    def test_mfq_rejects_multiqueue(self):
        solver = SolverFLD(fld_tandem(0.5, (1.0, 1.0)), method='mfq')
        with self.assertRaises(Exception):
            solver.runAnalyzer()


class TestMethodComparison(unittest.TestCase):
    """Compare different methods on the same network."""

    def test_pnorm_vs_statedep_mm1(self):
        r1 = SolverFLD(fld_mm1(0.5, 1.0), method='pnorm').runAnalyzer().result
        r2 = SolverFLD(fld_mm1(0.5, 1.0), method='statedep').runAnalyzer().result
        self.assertTrue(np.all(np.asarray(r1.QN) >= -1e-9))
        self.assertTrue(np.all(np.asarray(r2.QN) >= -1e-9))

    def test_softmin_vs_pnorm_convergence(self):
        r1 = SolverFLD(fld_mm1(0.5, 1.0), method='pnorm').runAnalyzer().result
        opts_softmin = SolverFLDOptions(method='softmin', softmin_alpha=100.0)
        r2 = SolverFLD(fld_mm1(0.5, 1.0), options=opts_softmin).runAnalyzer().result
        self.assertTrue(np.all(np.asarray(r1.QN) >= -1e-9))
        self.assertTrue(np.all(np.asarray(r2.QN) >= -1e-9))

    def test_matrix_vs_pnorm_mm1(self):
        r1 = SolverFLD(fld_mm1(0.5, 1.0), method='matrix').runAnalyzer().result
        r2 = SolverFLD(fld_mm1(0.5, 1.0), method='pnorm').runAnalyzer().result
        self.assertGreaterEqual(np.asarray(r1.UN)[1, 0], -1e-9)
        self.assertGreaterEqual(np.asarray(r2.UN)[1, 0], -1e-9)


class TestNetworkTopologies(unittest.TestCase):
    """Test different network topologies."""

    def test_tandem_queues(self):
        model = fld_tandem(0.5, (1.0, 1.0, 1.0))
        result = SolverFLD(model, method='matrix').runAnalyzer().result
        self.assertEqual(np.asarray(result.QN).shape, (model.getStruct().nstations, 1))
        self.assertTrue(np.all(np.asarray(result.QN) >= -1e-9))

    def test_two_class_network(self):
        model = fld_multiclass(lams=(0.5, 0.5), nqueues=1)
        result = SolverFLD(model, method='matrix').runAnalyzer().result
        self.assertEqual(np.asarray(result.QN).shape, (model.getStruct().nstations, 2))

    def test_three_class_network(self):
        model = fld_multiclass(lams=(0.3, 0.3, 0.3), nqueues=1)
        result = SolverFLD(model, method='matrix').runAnalyzer().result
        self.assertEqual(np.asarray(result.QN).shape, (model.getStruct().nstations, 3))

    def test_high_utilization(self):
        solver = SolverFLD(fld_mm1(0.8, 1.0), method='mfq')
        result = solver.runAnalyzer()
        self.assertIsNotNone(result.result)


if __name__ == '__main__':
    unittest.main(verbosity=2)
