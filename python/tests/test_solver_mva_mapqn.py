"""SolverMVA method 'amva.mapqn': the horizontal-cut mean value analysis (Casale-Smirni
DSN 2009 balances closed by the arrival theorem) for a closed model of one
exponential delay and one FCFS single-server queue with a MAP service per class.

The oracle is SolverCTMC on the same model, the accuracy asserted being the one
the method has (a few percent), plus the reference values of the prototype,
which every codebase must reproduce to solver precision.
"""
import unittest

import numpy as np

from line_solver import (Network, Queue, Delay, ClosedClass, SchedStrategy, Exp,
                         MAP, SolverMVA, SolverCTMC)
from line_solver.api.mapqn import mapqn_amva

D0a = np.array([[-3.0, 0.5], [0.2, -0.4]]); D1a = np.array([[2.5, 0.0], [0.2, 0.0]])
D0b = np.array([[-0.3, 0.3], [0.6, -0.6]]) - np.diag([1.5, 0.3]); D1b = np.diag([1.5, 0.3])


def _model(think, services, njobs, sched=SchedStrategy.FCFS):
    model = Network('mapqn')
    delay = Delay(model, 'Delay')
    queue = Queue(model, 'Queue', sched)
    for r, (z, svc, n) in enumerate(zip(think, services, njobs)):
        c = ClosedClass(model, 'C%d' % (r + 1), n, delay)
        delay.setService(c, Exp(1.0 / z))
        queue.setService(c, svc)
    model.link(Network.serialRouting(delay, queue))
    return model


def _cols(solver):
    t = solver.getAvgTable()
    return {c: t[c].to_numpy().astype(float) for c in ('QLen', 'Util', 'RespT', 'Tput')}


class TestAgainstTheReferenceValues(unittest.TestCase):
    """The prototype's numbers on the two models of the chapter; the delay row comes first."""

    def test_two_classes(self):
        m = _model([2.0, 8.0], [MAP(D0a, D1a), MAP(D0b, D1b)], [2, 2])
        c = _cols(SolverMVA(m, 'amva.mapqn'))
        np.testing.assert_allclose(c['Tput'][:2], [0.565971964, 0.1934693276], rtol=1e-7)
        np.testing.assert_allclose(c['QLen'][2:], [0.8680560719, 0.4522453793], rtol=1e-7)
        np.testing.assert_allclose(c['QLen'][:2], [1.1319439281, 1.5477546207], rtol=1e-7)

    def test_single_class(self):
        m = _model([2.0], [MAP(D0a, D1a)], [4])
        c = _cols(SolverMVA(m, 'amva.mapqn'))
        np.testing.assert_allclose(c['Tput'][0], 0.9717108597, rtol=1e-7)
        np.testing.assert_allclose(c['QLen'][1], 2.0565782806, rtol=1e-7)

    def test_direct_api(self):
        r = mapqn_amva([0.5, 0.125], [D0a, D0b], [D1a, D1b], [2, 2])
        np.testing.assert_allclose(r.X, [0.565971964, 0.1934693276], rtol=1e-7)
        np.testing.assert_allclose(r.Qq, [0.8680560719, 0.4522453793], rtol=1e-7)


class TestAgainstTheExactChain(unittest.TestCase):

    def test_within_the_method_accuracy(self):
        m = _model([2.0, 8.0], [MAP(D0a, D1a), MAP(D0b, D1b)], [2, 2])
        a = _cols(SolverMVA(m, 'amva.mapqn')); e = _cols(SolverCTMC(m))
        for col in ('QLen', 'Util', 'RespT', 'Tput'):
            np.testing.assert_allclose(a[col], e[col], rtol=0.10, err_msg=col)

    def test_exponential_service_is_exact_mva(self):
        m = _model([2.0, 8.0], [Exp(1.2), Exp(1.2)], [3, 3])
        a = _cols(SolverMVA(m, 'amva.mapqn')); e = _cols(SolverMVA(m, 'exact'))
        for col in ('QLen', 'Util', 'RespT', 'Tput'):
            np.testing.assert_allclose(a[col], e[col], rtol=1e-9, err_msg=col)


class TestTheGate(unittest.TestCase):

    def test_offered_and_supported_on_its_shape(self):
        s = SolverMVA(_model([2.0, 8.0], [MAP(D0a, D1a), MAP(D0b, D1b)], [2, 2]), 'amva.mapqn')
        self.assertIn('amva.mapqn', s.listValidMethods())
        ok, reason = s.supportsModelMethod('amva.mapqn')
        self.assertTrue(ok, reason)

    def test_withheld_off_its_shape(self):
        ps = SolverMVA(_model([2.0], [MAP(D0a, D1a)], [2], sched=SchedStrategy.PS), 'amva.mapqn')
        self.assertNotIn('amva.mapqn', ps.listValidMethods())
        ok, reason = ps.supportsModelMethod('amva.mapqn')
        self.assertFalse(ok)
        self.assertTrue('FCFS' in reason or 'SchedStrategy_PS' in reason, reason)
        with self.assertRaises(Exception):
            ps.getAvgTable()

    def test_the_default_method_still_takes_the_environment_route(self):
        m = _model([2.0], [MAP(D0a, D1a)], [4])
        s = SolverMVA(m)
        self.assertTrue(s.needsMapEnv(s.options))
        s2 = SolverMVA(m, 'amva.mapqn')
        self.assertFalse(s2.needsMapEnv(s2.options))
        s2.getAvgTable()
        self.assertEqual(s2._result['method'], 'amva.mapqn')

    def test_citation(self):
        s = SolverMVA(_model([2.0], [MAP(D0a, D1a)], [2]), 'amva.mapqn')
        keys = [c['key'] if isinstance(c, dict) else getattr(c, 'key', c) for c in s.citations()]
        self.assertTrue(any('CasSmi09' in str(k) for k in keys), keys)


if __name__ == '__main__':
    unittest.main()
