"""
Analytical validation tests for SolverFLD, using real LINE network models.

SolverFLD is a mean-field FLUID solver, so its steady state is the fluid fixed
point (lambda = mu * min(x, c)), not the stochastic M/M/1 law. For a stable open
single-server queue the fluid results are therefore:
    Util   = rho = lambda/mu
    QLen   = rho          (not rho/(1-rho))
    RespT  = 1/mu         (service time; no stochastic waiting in the fluid limit)
    Tput   = lambda
Little's law (QLen = Tput * RespT) holds exactly in the fluid limit.

Real open-network structs keep the Source as station 0, so the single queue is
station 1 (and tandem queues are stations 1..M).
"""

import unittest
import numpy as np

import line_solver as L
from line_solver.solvers.solver_fld import SolverFLD


def fld_mm1(lam=0.5, mu=1.0, nservers=1, sched=L.SchedStrategy.PS):
    """Open M/M/1 (or M/M/c): Source(0) -> Queue(1) -> Sink."""
    m = L.Network('mm1')
    s = L.Source(m, 'Source')
    q = L.Queue(m, 'Queue1', sched)
    if nservers != 1:
        q.setNumberOfServers(nservers)
    k = L.Sink(m, 'Sink')
    c = L.OpenClass(m, 'Class1')
    s.setArrival(c, L.Exp(lam))
    q.setService(c, L.Exp(mu))
    m.link(L.Network.serialRouting(s, q, k))
    return m


def fld_tandem(lam=0.5, mus=(1.0, 1.0), sched=L.SchedStrategy.PS):
    """Open tandem: Source(0) -> Q1(1) -> ... -> QM(M) -> Sink."""
    m = L.Network('tandem')
    s = L.Source(m, 'Source')
    queues = [L.Queue(m, 'Queue%d' % (i + 1), sched) for i in range(len(mus))]
    k = L.Sink(m, 'Sink')
    c = L.OpenClass(m, 'Class1')
    s.setArrival(c, L.Exp(lam))
    for i, q in enumerate(queues):
        q.setService(c, L.Exp(mus[i]))
    m.link(L.Network.serialRouting(s, *queues, k))
    return m


class TestMM1Analytical(unittest.TestCase):
    """M/M/1 fluid steady state (lam=0.5, mu=1.0, rho=0.5)."""

    def _solve(self, method='matrix'):
        solver = SolverFLD(fld_mm1(0.5, 1.0), method=method)
        solver.runAnalyzer()
        return solver.result

    def test_mm1_queue_length(self):
        # fluid QLen at the queue = rho
        self.assertAlmostEqual(np.asarray(self._solve().QN)[1, 0], 0.5, delta=0.1)

    def test_mm1_response_time(self):
        # fluid RespT at the queue = 1/mu
        self.assertAlmostEqual(np.asarray(self._solve().RN)[1, 0], 1.0, delta=0.1)

    def test_mm1_utilization(self):
        self.assertAlmostEqual(np.asarray(self._solve().UN)[1, 0], 0.5, delta=0.05)

    def test_mm1_throughput(self):
        self.assertAlmostEqual(np.asarray(self._solve().TN)[1, 0], 0.5, places=4)


class TestMM1Variants(unittest.TestCase):
    """M/M/1 fluid QLen = rho across loads."""

    def _qlen(self, lam):
        solver = SolverFLD(fld_mm1(lam, 1.0), method='matrix')
        solver.runAnalyzer()
        return np.asarray(solver.result.QN)[1, 0]

    def test_mm1_light_load(self):
        self.assertLess(self._qlen(0.1), 0.15)  # fluid QLen = rho = 0.1

    def test_mm1_moderate_load(self):
        self.assertAlmostEqual(self._qlen(0.5), 0.5, delta=0.15)  # fluid QLen = rho

    def test_mm1_high_load(self):
        self.assertAlmostEqual(self._qlen(0.9), 0.9, delta=0.15)  # fluid QLen = rho


class TestMMcAnalytical(unittest.TestCase):
    """M/M/c fluid behavior."""

    def test_mmc_single_server(self):
        solver = SolverFLD(fld_mm1(0.5, 1.0, nservers=1), method='matrix')
        solver.runAnalyzer()
        # fluid QLen = rho for c=1
        self.assertAlmostEqual(np.asarray(solver.result.QN)[1, 0], 0.5, delta=0.1)

    def test_mmc_improves_with_servers(self):
        s1 = SolverFLD(fld_mm1(1.0, 1.0, nservers=1), method='matrix'); s1.runAnalyzer()
        s2 = SolverFLD(fld_mm1(1.0, 1.0, nservers=2), method='matrix'); s2.runAnalyzer()
        L1 = np.asarray(s1.result.QN)[1, 0]
        L2 = np.asarray(s2.result.QN)[1, 0]
        self.assertLessEqual(L2, L1 + 1e-9)


class TestLittlesLaw(unittest.TestCase):
    """Little's Law (QLen = Tput * RespT) holds exactly in the fluid limit."""

    def test_littles_law_mm1(self):
        solver = SolverFLD(fld_mm1(0.5, 1.0), method='matrix'); solver.runAnalyzer()
        r = solver.result
        Lq = np.asarray(r.QN)[1, 0]
        W = np.asarray(r.RN)[1, 0]
        T = np.asarray(r.TN)[1, 0]
        self.assertAlmostEqual(Lq, T * W, delta=abs(T * W) * 0.02 + 1e-9)

    def test_littles_law_multistation(self):
        solver = SolverFLD(fld_tandem(0.5, (1.0, 1.0)), method='matrix'); solver.runAnalyzer()
        r = solver.result
        for i in (1, 2):  # queue stations
            Lq = np.asarray(r.QN)[i, 0]
            W = np.asarray(r.RN)[i, 0]
            T = np.asarray(r.TN)[i, 0]
            if T > 1e-10:
                self.assertAlmostEqual(Lq, T * W, delta=abs(T * W) * 0.05 + 1e-9,
                                       msg="Little's Law violated at station %d" % i)


class TestSystemMetrics(unittest.TestCase):
    """System-level fluid metrics."""

    def test_cycle_time_sum(self):
        solver = SolverFLD(fld_tandem(0.5, (1.0, 1.0)), method='matrix'); solver.runAnalyzer()
        r = solver.result
        cycle_time = np.asarray(r.CN)[0, 0]
        sum_resp = np.sum(np.asarray(r.RN)[[1, 2], 0])
        self.assertAlmostEqual(cycle_time, sum_resp, delta=abs(sum_resp) * 0.1 + 1e-9)

    def test_system_throughput_consistency(self):
        solver = SolverFLD(fld_mm1(0.5, 1.0), method='matrix'); solver.runAnalyzer()
        self.assertAlmostEqual(np.asarray(solver.result.XN)[0, 0], 0.5, places=3)


class TestStabilityBounds(unittest.TestCase):
    """Stability conditions and bounds."""

    def test_stability_condition_mm1(self):
        solver = SolverFLD(fld_mm1(0.5, 1.0), method='matrix'); solver.runAnalyzer()
        self.assertTrue(np.all(np.isfinite(np.asarray(solver.result.QN))))
        self.assertTrue(np.all(np.isfinite(np.asarray(solver.result.RN))))

    def test_queue_length_nonnegative(self):
        for lam in [0.1, 0.5, 0.8]:
            solver = SolverFLD(fld_mm1(lam, 1.0), method='matrix'); solver.runAnalyzer()
            self.assertTrue(np.all(np.asarray(solver.result.QN) >= -1e-9),
                            "Negative queue length for lambda=%s" % lam)

    def test_utilization_bounds(self):
        solver = SolverFLD(fld_tandem(0.5, (1.0, 1.0)), method='matrix'); solver.runAnalyzer()
        UN = np.asarray(solver.result.UN)
        self.assertTrue(np.all(UN >= -1e-9))
        self.assertTrue(np.all(UN <= 1.0 + 1e-9))


if __name__ == '__main__':
    unittest.main(verbosity=2)
