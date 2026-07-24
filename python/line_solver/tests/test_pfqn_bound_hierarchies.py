"""Tests for the hierarchical / multiserver / LD bound functions
(pfqn_pbh/pbk/bjbk/cbh/mcub/ssd/sib/ldbcmp) and single-class Marie, ported from
the MATLAB SolverBA bound suite. Bounds must bracket the exact solution.
"""

import unittest
import numpy as np

from line_solver.api.pfqn.bound_hierarchies import (
    pfqn_pbh, pfqn_pbk, pfqn_bjbk, pfqn_cbh, pfqn_mcub, pfqn_ssd,
    pfqn_sib, pfqn_ldbcmp,
)
from line_solver.api.pfqn.marie import pfqn_marie
from line_solver.api.pfqn.mva import pfqn_mva


def _exact_X(L, N, Z):
    XN = pfqn_mva(np.asarray(L, dtype=float).reshape(-1, 1),
                  np.array([N]), np.array([Z]))[0]
    return float(np.ravel(XN)[0])


class TestBoundHierarchies(unittest.TestCase):
    def setUp(self):
        self.L = np.array([2.0, 1.0, 3.0])
        self.N = 5
        self.Z = 1.0
        self.Xex = _exact_X(self.L, self.N, self.Z)

    def _bracket(self, lo, hi):
        self.assertLessEqual(lo, self.Xex + 1e-9)
        self.assertLessEqual(self.Xex, hi + 1e-9)

    def test_pbh_brackets_and_converges(self):
        for lvl in (1, 2, 4):
            lo, hi, _, _ = pfqn_pbh(self.L, self.N, self.Z, lvl)
            self._bracket(lo, hi)
        lo, hi, _, _ = pfqn_pbh(self.L, self.N, self.Z, self.N)
        self.assertAlmostEqual(lo, self.Xex, places=6)
        self.assertAlmostEqual(hi, self.Xex, places=6)

    def test_pbk_bjbk_bracket(self):
        self._bracket(*pfqn_pbk(self.L, self.N, self.Z, 2))
        self._bracket(*pfqn_bjbk(self.L, self.N, self.Z, 1))

    def test_cbh_brackets_and_exact_at_M(self):
        self._bracket(*pfqn_cbh(self.L, self.N, self.Z, 2))
        lo, hi = pfqn_cbh(self.L, self.N, self.Z, self.L.size)
        self.assertAlmostEqual(lo, self.Xex, places=6)
        self.assertAlmostEqual(hi, self.Xex, places=6)

    def test_mcub_multiclass_brackets(self):
        L = np.array([[1.0, 0.8], [1.2, 1.5]])
        N = np.array([2, 2])
        Z = np.array([1.0, 1.0])
        Xub, Xlb = pfqn_mcub(L, N, Z)
        Xex = np.ravel(pfqn_mva(L, N, Z)[0])
        self.assertTrue(np.all(Xlb <= Xex + 1e-9))
        self.assertTrue(np.all(Xex <= Xub + 1e-9))

    def test_ssd_returns_valid_interval(self):
        lo, hi = pfqn_ssd(np.array([2.0, 1.0]), 4, 0.0, np.array([2.0, 1.0]))
        self.assertLessEqual(lo, hi + 1e-9)
        self.assertGreater(lo, 0.0)

    def test_ldbcmp_qhat_anchor(self):
        Xlo, Rhi, Qhat = pfqn_ldbcmp(np.array([1.0, 0.8, 0.6]), 20, 0.0,
                                     np.array([0.0, 2.0, 0.0]))
        self.assertAlmostEqual(Qhat, 13.5, places=4)
        self.assertLessEqual(Xlo, _exact_X([1.0, 0.8, 0.6], 20, 0.0) + 1e-6)

    def test_sib_z0_brackets_and_rejects_delay(self):
        lo, hi, _, _ = pfqn_sib(self.L, self.N, 0.0, 3)
        Xex0 = _exact_X(self.L, self.N, 0.0)
        self.assertLessEqual(lo, Xex0 + 1e-9)
        self.assertLessEqual(Xex0, hi + 1e-9)
        with self.assertRaises(ValueError):
            pfqn_sib(self.L, self.N, 1.0, 3)


class TestMarieSingleClass(unittest.TestCase):
    def test_exponential_reduces_to_exact(self):
        L = np.array([1.0, 1.5])
        N, Z = 4, 2.0
        X, Q, U, C, it, mu = pfqn_marie(L, N, Z, np.array([1.0, 1.0]))
        XN, CN, QN = pfqn_mva(L.reshape(-1, 1), np.array([N]), np.array([Z]))[:3]
        self.assertLess(np.max(np.abs(Q - np.ravel(QN))), 1e-9)
        self.assertAlmostEqual(X, float(np.ravel(XN)[0]), places=9)

    def test_coxian_converges_finite(self):
        X, Q, U, C, it, mu = pfqn_marie(np.array([1.0, 1.5]), 4, 2.0,
                                        np.array([0.5, 2.0]))
        self.assertTrue(np.all(np.isfinite(Q)))
        self.assertGreater(X, 0.0)

    def test_multiclass_exponential_reduces_to_exact(self):
        # Class-independent exponential FCFS -> genuine BCMP -> exact MVA.
        L = np.array([[1.0, 1.0], [1.5, 1.5]])
        N, Z = np.array([2, 2]), np.array([2.0, 2.0])
        X, Q, U, C, it, mu = pfqn_marie(L, N, Z, np.ones((2, 2)))
        XN, CN, QN = pfqn_mva(L, N, Z)[:3]
        self.assertEqual(it, 0)
        self.assertLess(np.max(np.abs(Q - QN)), 1e-9)
        self.assertLess(np.max(np.abs(np.ravel(X) - np.ravel(XN))), 1e-9)

    def test_multiclass_coxian_converges_finite(self):
        L = np.array([[0.5, 0.7], [1.0, 0.8]])
        scv = np.array([[0.5, 2.0], [0.25, 3.0]])
        X, Q, U, C, it, mu = pfqn_marie(L, np.array([2, 2]),
                                        np.array([1.0, 1.0]), scv)
        self.assertTrue(np.all(np.isfinite(Q)))
        self.assertTrue(np.all(X > 0.0))
        self.assertEqual(Q.shape, (2, 2))


if __name__ == '__main__':
    unittest.main()
