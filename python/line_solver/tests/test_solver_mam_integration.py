"""
Integration tests for SolverMAM with real LINE network models.

Tests validate:
- Method routing and auto-selection with real networks
- Multi-class queueing networks
- Open vs closed network handling
- Result consistency across different solution methods
- Performance metrics correctness

These use real line_solver Network models (not hand-built NetworkStruct doubles),
so the result arrays follow the real struct convention: an open network's Source
is station 0 and the queues follow at stations 1..M.
"""

import unittest
import numpy as np
import pandas as pd
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import line_solver as L
from line_solver.solvers.solver_mam import SolverMAM, SolverMAMOptions


def mm1_open(lam=1.0, mu=2.0):
    """Open M/M/1: Source -> Queue -> Sink. Stations: [Source(0), Queue(1)]."""
    m = L.Network('mm1')
    s = L.Source(m, 'Source')
    q = L.Queue(m, 'Queue1', L.SchedStrategy.FCFS)
    k = L.Sink(m, 'Sink')
    c = L.OpenClass(m, 'Class1')
    s.setArrival(c, L.Exp(lam))
    q.setService(c, L.Exp(mu))
    m.link(L.Network.serialRouting(s, q, k))
    return m


def tandem_open(lam=1.0, mus=(2.0, 2.0)):
    """Open tandem: Source -> Q1 -> ... -> QM -> Sink.
    Stations: [Source(0), Q1(1), ..., QM(M)]."""
    m = L.Network('tandem')
    s = L.Source(m, 'Source')
    queues = [L.Queue(m, 'Queue%d' % (i + 1), L.SchedStrategy.FCFS) for i in range(len(mus))]
    k = L.Sink(m, 'Sink')
    c = L.OpenClass(m, 'Class1')
    s.setArrival(c, L.Exp(lam))
    for i, q in enumerate(queues):
        q.setService(c, L.Exp(mus[i]))
    m.link(L.Network.serialRouting(s, *queues, k))
    return m


def multiclass_tandem(lams=(0.5, 0.5), mus=((2.0, 1.5), (2.0, 1.5))):
    """Open 2-class tandem over len(mus) queues. mus[i] = (rate_class0, rate_class1).
    Stations: [Source(0), Q1(1), ..., QM(M)]."""
    K = len(lams)
    m = L.Network('mc')
    s = L.Source(m, 'Source')
    queues = [L.Queue(m, 'Queue%d' % (i + 1), L.SchedStrategy.FCFS) for i in range(len(mus))]
    k = L.Sink(m, 'Sink')
    classes = [L.OpenClass(m, 'Class%d' % (r + 1)) for r in range(K)]
    for r in range(K):
        s.setArrival(classes[r], L.Exp(lams[r]))
        for i, q in enumerate(queues):
            q.setService(classes[r], L.Exp(mus[i][r]))
    m.link(L.Network.serialRouting(s, *queues, k))
    return m


def closed_cycle(njobs=5, mus=(2.0, 2.0)):
    """Closed 2-station cycle: Delay <-> Queue. Stations: [Delay(0), Queue(1)]."""
    m = L.Network('closed')
    d = L.Delay(m, 'Delay')
    q = L.Queue(m, 'Queue1', L.SchedStrategy.FCFS)
    c = L.ClosedClass(m, 'Class1', njobs, d)
    d.setService(c, L.Exp(mus[0]))
    q.setService(c, L.Exp(mus[1]))
    m.link(L.Network.serialRouting(d, q))
    return m


def closed_single(njobs=5, mu=2.0, think=1.0):
    """Single-queue closed model: Delay <-> Queue with one queueing station."""
    return closed_cycle(njobs=njobs, mus=(1.0 / think, mu))


def _queue_rows(sn):
    """Indices of queueing stations (exclude the Source of an open model)."""
    from line_solver.constants import NodeType
    rows = []
    for i in range(sn.nstations):
        nd = int(sn.stationToNode[i])
        if sn.nodetype[nd] != NodeType.Source:
            rows.append(i)
    return rows


class TestM_M_1Network(unittest.TestCase):
    """Test M/M/1 queue (Source -> Queue -> Sink), lam=1, mu=2, rho=0.5."""

    def setUp(self):
        self.model = mm1_open(lam=1.0, mu=2.0)

    def test_dec_source_m_m_1(self):
        solver = SolverMAM(self.model, method='dec.source')
        solver.runAnalyzer()
        self.assertIsNotNone(solver.result)
        # Queue is station 1 (station 0 is the Source): rho = 0.5
        util = np.asarray(solver.result.UN)[1, 0]
        self.assertAlmostEqual(util, 0.5, places=1)
        qlen = np.asarray(solver.result.QN)[1, 0]
        self.assertGreater(qlen, 0)
        self.assertLess(qlen, 2)

    def test_mna_open_m_m_1(self):
        solver = SolverMAM(self.model, method='mna_open')
        solver.runAnalyzer()
        self.assertIsNotNone(solver.result)
        self.assertEqual(solver.result.method, 'mna_open')
        self.assertGreater(solver.result.totiter, 0)

    def test_inap_m_m_1(self):
        solver = SolverMAM(self.model, method='inap')
        solver.runAnalyzer()
        self.assertIsNotNone(solver.result)
        self.assertEqual(solver.result.method, 'inap')


class TestTandemNetwork(unittest.TestCase):
    """Test 2-station tandem queue: Source -> Q1 -> Q2 -> Sink."""

    def setUp(self):
        self.model = tandem_open(lam=1.0, mus=(2.0, 2.0))

    def test_dec_source_tandem(self):
        solver = SolverMAM(self.model, method='dec.source')
        solver.runAnalyzer()
        self.assertIsNotNone(solver.result)
        # 3 stations (Source, Q1, Q2), 1 class
        self.assertEqual(np.asarray(solver.result.QN).shape, (3, 1))
        UN = np.asarray(solver.result.UN)
        for i in (1, 2):  # the two queues
            self.assertGreater(UN[i, 0], 0)
            self.assertLess(UN[i, 0], 1.0)

    def test_method_consistency_tandem(self):
        for method in ['dec.source', 'mna_open', 'inap']:
            solver = SolverMAM(self.model, method=method)
            solver.runAnalyzer()
            QN = np.asarray(solver.result.QN)
            self.assertFalse(np.any(np.isnan(QN)))
            self.assertTrue(np.all(QN >= 0))


class TestMultiClassNetwork(unittest.TestCase):
    """Test multi-class (2-class) open tandem network."""

    def setUp(self):
        self.model = multiclass_tandem(lams=(0.5, 0.5), mus=((2.0, 1.5), (2.0, 1.5)))

    def test_dec_source_multiclass(self):
        solver = SolverMAM(self.model, method='dec.source')
        solver.runAnalyzer()
        # 3 stations (Source, Q1, Q2), 2 classes
        self.assertEqual(np.asarray(solver.result.QN).shape, (3, 2))
        self.assertEqual(np.asarray(solver.result.UN).shape, (3, 2))
        self.assertEqual(np.asarray(solver.result.RN).shape, (3, 2))
        self.assertEqual(np.asarray(solver.result.TN).shape, (3, 2))

    def test_mna_multiclass(self):
        solver = SolverMAM(self.model, method='mna_open')
        solver.runAnalyzer()
        tput = np.asarray(solver.getTput()).ravel()
        self.assertGreaterEqual(tput.size, 2)  # at least per-class
        self.assertTrue(np.all(tput >= -1e-9))


class TestClosedNetwork(unittest.TestCase):
    """Test closed networks with fixed population."""

    def setUp(self):
        self.model = closed_cycle(njobs=5, mus=(2.0, 2.0))

    def test_mna_closed(self):
        solver = SolverMAM(self.model, method='mna_closed')
        solver.runAnalyzer()
        total_qlen = np.sum(np.asarray(solver.result.QN))
        self.assertAlmostEqual(total_qlen, 5, delta=0.5)

    def test_dec_source_closed(self):
        solver = SolverMAM(self.model, method='dec.source')
        solver.runAnalyzer()
        self.assertIsNotNone(solver.result)
        self.assertFalse(np.any(np.isnan(np.asarray(solver.result.QN))))


class TestAutoSelection(unittest.TestCase):
    """Test automatic method selection."""

    def test_auto_select_open_network(self):
        solver = SolverMAM(tandem_open(), method='default')
        self.assertEqual(solver._select_method(), 'dec.source')

    def test_auto_select_closed_1station(self):
        solver = SolverMAM(closed_single(njobs=5), method='default')
        self.assertIn(solver._select_method(), ['ldqbd', 'dec.source'])

    def test_mna_auto_select_open(self):
        solver = SolverMAM(tandem_open(), method='mna')
        self.assertEqual(solver._select_method(), 'mna_open')

    def test_mna_auto_select_closed(self):
        solver = SolverMAM(closed_cycle(njobs=5), method='mna')
        self.assertEqual(solver._select_method(), 'mna_closed')


class TestResultAccessors(unittest.TestCase):
    """Test all result accessor methods (2-class open tandem, 3 stations)."""

    def setUp(self):
        self.model = multiclass_tandem(lams=(0.5, 0.5), mus=((2.0, 1.5), (2.0, 1.5)))
        self.solver = SolverMAM(self.model, method='dec.source')
        self.solver.runAnalyzer()
        self.nst = 3

    def test_get_avg_qlen(self):
        qlen = np.asarray(self.solver.getAvgQLen())
        self.assertEqual(qlen.shape[0], self.nst)
        self.assertTrue(np.all(qlen >= 0))

    def test_get_avg_util(self):
        util = np.asarray(self.solver.getAvgUtil())
        self.assertEqual(util.shape[0], self.nst)
        self.assertTrue(np.all(util >= 0))
        self.assertTrue(np.all(util <= 1.1))

    def test_get_avg_respt(self):
        resp_t = np.asarray(self.solver.getAvgRespT())
        self.assertEqual(resp_t.shape[0], self.nst)
        self.assertTrue(np.all(resp_t >= 0))

    def test_get_tput(self):
        tput = np.asarray(self.solver.getTput()).ravel()
        self.assertGreaterEqual(tput.size, self.nst)
        self.assertTrue(np.all(tput >= -1e-9))

    def test_get_avg_table(self):
        table = self.solver.getAvgTable()
        self.assertIsNotNone(table)
        rendered = str(table)
        for col in ('QLen', 'Util', 'RespT'):
            self.assertIn(col, rendered)


class TestStaticMethods(unittest.TestCase):
    """Test static introspection methods."""

    def test_list_valid_methods(self):
        methods = SolverMAM.listValidMethods()
        for method in ['default', 'dec.source', 'dec.mmap', 'dec.poisson',
                       'mna', 'mna_open', 'mna_closed', 'ldqbd', 'inap', 'inapplus']:
            self.assertIn(method, methods)

    def test_supports_valid_method(self):
        sn = tandem_open().getStruct()
        can_solve, reason = SolverMAM.supports(sn, 'dec.source')
        self.assertTrue(can_solve)

    def test_supports_invalid_method(self):
        sn = tandem_open().getStruct()
        can_solve, reason = SolverMAM.supports(sn, 'nonexistent_method')
        self.assertFalse(can_solve)

    def test_get_feature_set(self):
        features = SolverMAM.getFeatureSet()
        self.assertIn('OpenClass', features)
        self.assertIn('ClosedClass', features)
        self.assertIn('SchedStrategy_FCFS', features)
        self.assertIn('MAP', features)
        self.assertIn('Retrial', features)

    def test_default_options(self):
        opts = SolverMAM.defaultOptions()
        self.assertIsInstance(opts, SolverMAMOptions)
        self.assertEqual(opts.method, 'default')
        self.assertGreater(opts.max_iter, 0)
        self.assertGreater(opts.tol, 0)


class TestConvergenceBehavior(unittest.TestCase):
    """Test convergence properties with different tolerance levels."""

    def test_convergence_tight_tolerance(self):
        opts = SolverMAMOptions(method='mna_open', tol=1e-8, max_iter=200)
        solver = SolverMAM(tandem_open(), options=opts)
        solver.runAnalyzer()
        self.assertLessEqual(solver.result.totiter, opts.max_iter)

    def test_convergence_loose_tolerance(self):
        opts = SolverMAMOptions(method='mna_open', tol=1e-3, max_iter=50)
        solver = SolverMAM(tandem_open(), options=opts)
        solver.runAnalyzer()
        self.assertLessEqual(solver.result.totiter, opts.max_iter)


class TestNumericalEdgeCases(unittest.TestCase):
    """Test solver behavior at numerical edge cases."""

    def test_very_low_arrival_rate(self):
        solver = SolverMAM(tandem_open(lam=1e-4, mus=(2.0, 2.0)), method='dec.source')
        solver.runAnalyzer()
        QN = np.asarray(solver.result.QN)
        self.assertFalse(np.any(np.isnan(QN)))
        self.assertTrue(np.all(QN >= 0))

    def test_high_utilization_stability(self):
        solver = SolverMAM(mm1_open(lam=0.95, mu=1.0), method='mna_open')
        solver.runAnalyzer()
        util = np.asarray(solver.result.UN)[1, 0]
        self.assertLess(util, 1.1)

    def test_variable_service_rates(self):
        solver = SolverMAM(tandem_open(lam=0.9, mus=(1.0, 2.0, 3.0)), method='dec.source')
        solver.runAnalyzer()
        self.assertIsNotNone(solver.result)
        self.assertFalse(np.any(np.isnan(np.asarray(solver.result.QN))))


class TestPerformanceMetrics(unittest.TestCase):
    """Test that computed performance metrics make sense."""

    def test_queueing_theory_relationships(self):
        model = tandem_open(lam=1.0, mus=(2.0, 2.0))
        solver = SolverMAM(model, method='dec.source')
        solver.runAnalyzer()
        sn = model.getStruct()
        QN = np.asarray(solver.result.QN)
        RN = np.asarray(solver.result.RN)
        TN = np.asarray(solver.result.TN)
        for m in _queue_rows(sn):
            q_expected = TN[m, 0] * RN[m, 0]  # Little's law per station
            q_actual = QN[m, 0]
            if q_expected > 0.1:
                rel_error = abs(q_actual - q_expected) / q_expected
                self.assertLess(rel_error, 0.2,
                                "Station %d: Little's law violated by %.1f%%" % (m, rel_error * 100))

    def test_utilization_bounds(self):
        model = multiclass_tandem(lams=(0.5, 0.5), mus=((2.0, 1.5), (2.0, 1.5), (2.0, 1.5)))
        for method in ['dec.source', 'mna_open', 'inap']:
            solver = SolverMAM(model, method=method)
            solver.runAnalyzer()
            UN = np.asarray(solver.result.UN)
            self.assertTrue(np.all(UN >= -0.01))
            self.assertTrue(np.all(UN <= 1.01))


if __name__ == '__main__':
    unittest.main(verbosity=2)
