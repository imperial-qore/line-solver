"""Regression tests for the qsys solver paths added in the Hillier-Yu benchmark cases.

Each test builds a single-class open Source-Queue-Sink network at rho=0.99
and asserts the expected Lq from the corresponding LINE solver path. Tolerance
is 5e-4 absolute (values are stable to >=4 decimal places across all three
LINE codebases). Reference values were verified independently against the
analytical / matrix-geometric formulas during implementation.
"""

import os
import sys
import unittest
import warnings

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from line_solver import (
    Network,
    Source,
    Queue,
    Sink,
    OpenClass,
    SchedStrategy,
    Exp,
    Det,
    Erlang,
    SolverMVA,
    SolverMAM,
)


def _build(name, c, build_arrival, build_service):
    m = Network(name)
    src = Source(m, "Source")
    q = Queue(m, "Queue1", SchedStrategy.FCFS)
    q.setNumberOfServers(c)
    Sink(m, "Sink")
    oc = OpenClass(m, "Class1")
    build_arrival(src, oc)
    build_service(q, oc)
    P = m.initRoutingMatrix()
    P.set(oc, oc, [[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    m.link(P)
    return m


def _lq(solver, c):
    Q = np.asarray(solver.getAvgQLen()).flatten()
    U = np.asarray(solver.getAvgUtil()).flatten()
    return float(Q[1] - c * U[1])


class TestHillierYuQsys(unittest.TestCase):
    """End-to-end Lq regression at rho=0.99 across the new dispatch paths."""

    TOL = 5e-4

    def setUp(self):
        # Suppress the noisy "Reducible network topology" warning emitted on
        # every model with a Sink reachable only through a single queue.
        warnings.filterwarnings("ignore", message=".*[Rr]educible.*")
        warnings.filterwarnings("ignore", message=".*matlib is.*")

    # ---------------- M/M/k via MVA exact ---------------------------------
    def test_mm4_mva_exact(self):
        m = _build(
            "MM4", 4,
            lambda src, oc: src.setArrival(oc, Exp(3.96)),
            lambda q, oc: q.setService(oc, Exp(1.0)),
        )
        s = SolverMVA(m, method="exact")
        self.assertAlmostEqual(_lq(s, 4), 96.8126, delta=self.TOL)

    # ---------------- M/G/1 (M/E3/1) via MVA exact (PK) -------------------
    def test_me3_1_mva_exact(self):
        m = _build(
            "ME3_1", 1,
            lambda src, oc: src.setArrival(oc, Exp(0.99)),
            lambda q, oc: q.setService(oc, Erlang.fitMeanAndOrder(1.0, 3)),
        )
        s = SolverMVA(m, method="exact")
        self.assertAlmostEqual(_lq(s, 1), 65.3400, delta=self.TOL)

    # ---------------- M/E3/c via MAM (rate-scaling + surrogate) -----------
    def test_me3_3_mam(self):
        m = _build(
            "ME3_3", 3,
            lambda src, oc: src.setArrival(oc, Exp(2.97)),
            lambda q, oc: q.setService(oc, Erlang.fitMeanAndOrder(1.0, 3)),
        )
        s = SolverMAM(m)
        # 64.8036799734 is the EXACT M/E3/3 Lq, from the MAP/PH/c QBD of
        # qsys_mapphc in C++ (double and real:50) and in MATLAB, all agreeing
        # with SolverMAM to 11 digits and with Little's law on Wq. The 65.3400
        # this asserted until 2026-08-16 is the c=1 value of the test above:
        # Pollaczek-Khinchine at lambda=0.99 with E[S^2]=4/3 gives exactly
        # 0.9801*(4/3)/(2*0.01) = 65.34, a different queue from this one.
        self.assertAlmostEqual(_lq(s, 3), 64.8036799734, delta=self.TOL)

    # ---------------- M/D/c via MAM (Crommelin embedded DTMC) -------------
    def test_md_3_mam_crommelin(self):
        m = _build(
            "MD_3", 3,
            lambda src, oc: src.setArrival(oc, Exp(2.97)),
            lambda q, oc: q.setService(oc, Det(1.0)),
        )
        s = SolverMAM(m)
        self.assertAlmostEqual(_lq(s, 3), 48.6193, delta=self.TOL)

    # ---------------- D/M/c via MAM (Smith embedded DTMC) -----------------
    def test_dm_4_mam_smith(self):
        m = _build(
            "DM_4", 4,
            lambda src, oc: src.setArrival(oc, Det(1.0 / 3.96)),
            lambda q, oc: q.setService(oc, Exp(1.0)),
        )
        s = SolverMAM(m)
        self.assertAlmostEqual(_lq(s, 4), 47.8380, delta=self.TOL)

    # ---------------- PH/M/c via MAM (matrix-geometric) -------------------
    def _erlang_m_c(self, c, expected_lq):
        m = _build(
            f"E3M{c}", c,
            lambda src, oc: src.setArrival(oc, Erlang.fitMeanAndOrder(1.0 / (c * 0.99), 3)),
            lambda q, oc: q.setService(oc, Exp(1.0)),
        )
        s = SolverMAM(m)
        self.assertAlmostEqual(_lq(s, c), expected_lq, delta=self.TOL)

    def test_e3_m_1(self):  self._erlang_m_c(1, 65.1206)
    def test_e3_m_2(self):  self._erlang_m_c(2, 64.7159)
    def test_e3_m_3(self):  self._erlang_m_c(3, 64.4035)
    def test_e3_m_4(self):  self._erlang_m_c(4, 64.1398)
    def test_e3_m_5(self):  self._erlang_m_c(5, 63.9075)
    def test_e3_m_8(self):  self._erlang_m_c(8, 63.3256)
    def test_e3_m_10(self): self._erlang_m_c(10, 62.9986)

    # ---------------- E3/M/1 via MVA exact (qsys_phm1 path) ---------------
    def test_e3_m_1_mva_exact(self):
        m = _build(
            "E3M1_MVA", 1,
            lambda src, oc: src.setArrival(oc, Erlang.fitMeanAndOrder(1.0 / 0.99, 3)),
            lambda q, oc: q.setService(oc, Exp(1.0)),
        )
        s = SolverMVA(m, method="exact")
        self.assertAlmostEqual(_lq(s, 1), 65.1206, delta=self.TOL)


if __name__ == "__main__":
    unittest.main()
