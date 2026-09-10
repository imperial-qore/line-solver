"""Cross-language regression for the two Bertsimas open-network bounds.

``npfqn_bnd_bpt`` is the first-order LP relaxation of the achievable region
(Bertsimas-Paschalidis-Tsitsiklis, Ann. Appl. Prob. 4(1), 1994) and
``npfqn_bnd_bgt`` is the piecewise-linear Lyapunov bound (Bertsimas-Gamarnik-
Tsitsiklis, Ann. Appl. Prob. 11(4), 2001). The reference values are MATLAB's
(matlab/src/api/npfqn/) and are reproduced by the JAR and the C++ ports to the
digits asserted here; see _kb/06-solver-catalog.md under ``bpt.lower`` and
``bgt.upper``.
"""

import unittest

import numpy as np

from line_solver import (Exp, Network, OpenClass, SchedStrategy, Queue, Sink,
                         SolverBA, Source)
from line_solver.api.npfqn import npfqn_bnd_bgt, npfqn_bnd_bpt


class TestNpfqnBndBpt(unittest.TestCase):

    def test_exact_on_mm1(self):
        """The relaxation is tight on M/M/1 at every load."""
        for rho in (0.3, 0.5, 0.7, 0.9):
            z = npfqn_bnd_bpt([rho], [1.0], np.zeros((1, 1)), [0], [1.0]).zlb
            self.assertAlmostEqual(z, 1.0 / (1.0 - rho), places=9)

    def test_two_classes_one_station(self):
        """The bound sits below the cmu-optimal achievable point, 3.6875."""
        P = np.zeros((2, 2))
        l0, mu, st = [0.4, 0.4], [2.0, 1.0], [0, 0]
        self.assertAlmostEqual(npfqn_bnd_bpt(l0, mu, P, st, [1, 1]).zlb, 3.4375, places=9)
        self.assertAlmostEqual(npfqn_bnd_bpt(l0, mu, P, st, [1, 0]).zlb, 0.625, places=9)
        self.assertAlmostEqual(npfqn_bnd_bpt(l0, mu, P, st, [0, 1]).zlb, 1.0 / 0.6, places=9)
        self.assertLessEqual(npfqn_bnd_bpt(l0, mu, P, st, [1, 1]).zlb, 3.6875 + 1e-12)

    def test_tandem_first_station_is_exact(self):
        """Station 1 is an M/M/1 in isolation; station 2 falls back to 1/mu."""
        P = np.zeros((2, 2))
        P[0, 1] = 1.0
        l0, mu, st = [0.5, 0.0], [1.0, 1.0], [0, 1]
        self.assertAlmostEqual(npfqn_bnd_bpt(l0, mu, P, st, [1, 0]).zlb, 2.0, places=9)
        self.assertAlmostEqual(npfqn_bnd_bpt(l0, mu, P, st, [0, 1]).zlb, 1.0, places=9)

    def test_saturated_station_is_refused(self):
        with self.assertRaises(ValueError):
            npfqn_bnd_bpt([1.2], [1.0], np.zeros((1, 1)), [0], [1.0])

    def test_solver_ba_reproduces_the_exact_mm1(self):
        model = Network('mm1')
        src = Source(model, 'Src')
        q = Queue(model, 'Q1', SchedStrategy.FCFS)
        snk = Sink(model, 'Snk')
        c1 = OpenClass(model, 'C1')
        src.setArrival(c1, Exp(0.7))
        q.setService(c1, Exp(1.0))
        model.link(Network.serialRouting(src, q, snk))
        Q, U, R, T = SolverBA(model, 'bpt.lower').getAvg()[:4]
        self.assertAlmostEqual(float(R[1, 0]), 1.0 / 0.3, places=8)
        self.assertAlmostEqual(float(Q[1, 0]), 0.7 / 0.3, places=8)


class TestNpfqnBndBgt(unittest.TestCase):

    def test_mm1_gamma_has_a_closed_form(self):
        """With J = 1 the uniformized drift of L = 1 gives (mu-lambda)/(lambda+mu)."""
        qref = {0.3: 23.9340659340659, 0.5: 32.6666666666667,
                0.7: 53.6862745098039, 0.9: 160.105263157895}
        for rho, qub in qref.items():
            r = npfqn_bnd_bgt([rho], [[1.0]], [[0]], 1)
            self.assertAlmostEqual(r.gamma, (1 - rho) / (1 + rho), places=9)
            self.assertAlmostEqual(float(r.Qub[0][0]), qub, places=7)
            # it IS an upper bound on the exact mean queue length
            self.assertGreaterEqual(float(r.Qub[0][0]), rho / (1 - rho))

    def test_tandem(self):
        r = npfqn_bnd_bgt([0.5], [[1.0, 1.0]], [[0, 1]], 2)
        self.assertAlmostEqual(r.gamma, 0.2, places=9)
        self.assertAlmostEqual(r.Lmax, 1.0, places=9)
        self.assertAlmostEqual(r.B, 5529.6, places=6)
        self.assertAlmostEqual(r.U, 5578.0, places=6)
        self.assertAlmostEqual(r.tail_ratio, 1.1 / 1.15, places=9)

    def test_globally_unstable_lu_kumar_is_refused(self):
        """Every station is at 0.7, yet rho_2 + rho_4 = 1.2 > 1 (Rybko-Stolyar).

        A per-station load test would accept this network; GLP[dm] must not.
        """
        with self.assertRaises(ValueError):
            npfqn_bnd_bgt([1.0], [[1 / 0.1, 1 / 0.6, 1 / 0.1, 1 / 0.6]], [[0, 1, 1, 0]], 2)

    def test_stable_lu_kumar(self):
        r = npfqn_bnd_bgt([1.0], [[1 / 0.3, 1 / 0.6, 1 / 0.3, 1 / 0.1]], [[0, 1, 1, 0]], 2)
        self.assertAlmostEqual(r.gamma, 0.00574712643678161, places=12)
        self.assertAlmostEqual(float(r.rho_station[0]), 0.4, places=12)
        self.assertAlmostEqual(float(r.rho_station[1]), 0.9, places=12)

    def test_solver_ba_brackets_the_tandem(self):
        """bgt.upper and bpt.lower bracket the exact solution, loosely above."""
        model = Network('tandem')
        src = Source(model, 'Src')
        q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
        q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
        snk = Sink(model, 'Snk')
        c1 = OpenClass(model, 'C1')
        src.setArrival(c1, Exp(0.5))
        q1.setService(c1, Exp(1.0))
        q2.setService(c1, Exp(1.0))
        model.link(Network.serialRouting(src, q1, q2, snk))
        Qup = SolverBA(model, 'bgt.upper').getAvg()[0]
        Qlo = SolverBA(model, 'bpt.lower').getAvg()[0]
        for i in (1, 2):
            self.assertLessEqual(float(Qlo[i, 0]), 1.0 + 1e-9)
            self.assertGreaterEqual(float(Qup[i, 0]), 1.0)
        self.assertAlmostEqual(float(Qup[1, 0]), 5578.0, places=6)

    def test_probabilistic_split_is_refused_by_name(self):
        model = Network('split')
        src = Source(model, 'Src')
        q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
        q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
        snk = Sink(model, 'Snk')
        c1 = OpenClass(model, 'C1')
        src.setArrival(c1, Exp(0.5))
        q1.setService(c1, Exp(1.0))
        q2.setService(c1, Exp(1.0))
        P = model.initRoutingMatrix()
        P.set(c1, c1, src, q1, 0.5)
        P.set(c1, c1, src, q2, 0.5)
        P.set(c1, c1, q1, snk, 1.0)
        P.set(c1, c1, q2, snk, 1.0)
        model.link(P)
        with self.assertRaises(ValueError):
            SolverBA(model, 'bgt.upper').getAvg()


if __name__ == '__main__':
    unittest.main()
