"""Tests for the hierarchical / multiserver / LD bound functions
(pfqn_pbh/pbk/bjbk/cbh/mcub/ssd/sib/ldbcmp) and single-class Marie, ported from
the MATLAB SolverBA bound suite. Bounds must bracket the exact solution.

The scb family is checked against the published tables of Dowdy et al. (1992)
rather than against a bracket on the given model, because that is not what it
bounds: it bounds the multiclass system the single-class model aggregates.
"""

import unittest
import numpy as np

from line_solver.api.pfqn.bound_hierarchies import (
    pfqn_pbh, pfqn_pbk, pfqn_bjbk, pfqn_cbh, pfqn_mcub, pfqn_ssd,
    pfqn_sib, pfqn_ldbcmp, pfqn_scb, pfqn_scbgap, pfqn_usumbound,
    pfqn_minclasses,
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


class TestScb(unittest.TestCase):
    """Dowdy, Carlson, Krantz, Tripathi (1992), JACM 39(1):188-213."""

    # Table I, p.200: maximum percent error of approximating a multiclass
    # system by its single-class counterpart, rows N=1..5, columns K=1..5.
    TABLE_I = [[0, 0, 0, 0, 0], [0, 33, 33, 33, 33], [0, 25, 40, 40, 40],
               [0, 20, 33, 43, 43], [0, 17, 29, 38, 44]]

    def test_scbgap_reproduces_table_i(self):
        for N in range(1, 6):
            for K in range(1, 6):
                self.assertEqual(round(100 * pfqn_scbgap(N, K)),
                                 self.TABLE_I[N - 1][K - 1], msg='N=%d K=%d' % (N, K))
        # the 50% ceiling of Theorem 3, approached as N = K -> inf
        self.assertLess(pfqn_scbgap(1000, 1000), 0.5)
        self.assertGreater(pfqn_scbgap(1000, 1000), 0.499)
        # merging a single class changes nothing
        self.assertEqual(pfqn_scbgap(6, 4, 1), 0.0)

    def test_scbgap_undominated_domain(self):
        # strictly tighter below r = K, equal at r = K, refused past it
        for r in range(2, 5):
            self.assertLess(pfqn_scbgap(8, 5, r, True), pfqn_scbgap(8, 5, r, False))
        self.assertAlmostEqual(pfqn_scbgap(8, 5, 5, True), pfqn_scbgap(8, 5, 5, False), places=12)
        with self.assertRaises(ValueError):
            pfqn_scbgap(8, 5, 6, True)

    def test_usumbound_and_minclasses_section_4_7(self):
        # the paper's worked example: K=2 devices, N=3 customers
        self.assertAlmostEqual(pfqn_usumbound(1, 2, 3), 1.5, places=12)
        self.assertAlmostEqual(pfqn_usumbound(2, 2, 3), 2.0, places=12)
        self.assertEqual(pfqn_minclasses(1.6, 2, 3), 2)
        self.assertEqual(pfqn_minclasses(1.4, 2, 3), 1)
        self.assertTrue(np.isnan(pfqn_minclasses(2.5, 2, 3)))
        # nondecreasing in R, which is what makes the inversion well posed
        for R in range(1, 5):
            self.assertLessEqual(pfqn_usumbound(R, 4, 5), pfqn_usumbound(R + 1, 4, 5))

    def test_scb_brackets_the_section_2_example(self):
        # single-class demands (0.114, 0.040, 0.062) at N = 4; the paper's
        # multiclass counterpart runs at 8.761 and the aggregate at 8.152
        L = np.array([0.114, 0.040, 0.062])
        Xlo, Xhi, Ulo, Uhi = pfqn_scb(L, 4)
        self.assertAlmostEqual(Xlo, 8.151894490678735, places=9)
        self.assertLessEqual(Xlo, 8.7615)
        self.assertGreaterEqual(Xhi, 8.7615)
        # the lower side IS the exact single-class solution, not an approximation
        self.assertAlmostEqual(Xlo, float(np.ravel(pfqn_mva(L.reshape(-1, 1), 4)[0])[0]), places=9)
        # Corollary 1: the utilization ratio is uniform across devices
        self.assertLess(np.max(np.abs(Ulo - Xlo * L)), 1e-12)
        self.assertLess(np.max(np.abs(Uhi / Ulo - Xhi / Xlo)), 1e-9)
        # the single-server cap binds, so the busiest device sits exactly at 1
        self.assertAlmostEqual(Uhi[0], 1.0, places=12)
        # N = 1 admits no aggregation error, so the bracket collapses
        b1 = pfqn_scb(L, 1)
        self.assertAlmostEqual(b1[1], b1[0], places=12)

    def test_scb_brackets_every_table_ii_row(self):
        # Each row is a multiclass demand matrix with one customer per class.
        # The exactly-weighted single-class aggregate must bracket the exact
        # multiclass throughput, and the realized error must respect the gap.
        rows = [
            (2, [[20, 0], [0, 8]]),
            (2, [[15, 5], [1, 7]]),
            (2, [[20, 0], [0, 20], [3.032, 18.069], [5.976, 12.356],
                 [17.007, 19.704], [19.678, 11.736]]),
            (2, [[4.491, 12.272], [14.381, 4.381], [2.589, 2.521],
                 [0.669, 11.404], [14.383, 12.737], [8.771, 11.055]]),
            (3, [[150, 0, 0], [0, 150, 0], [0, 0, 150],
                 [91.79, 78.68, 85.26], [17.43, 29.00, 118.42]]),
            (3, [[150, 0, 0], [80.24, 57.68, 126.27], [141.74, 108.17, 148.98],
                 [135.23, 45.73, 40.54], [84.12, 38.56, 44.66]]),
        ]
        for K, D in rows:
            D = np.asarray(D, dtype=float)
            R = D.shape[0]
            Xr = np.ravel(pfqn_mva(D.T, np.ones(R, dtype=int))[0])
            XR = float(np.sum(Xr))
            Dsc = (Xr @ D) / XR              # the paper's exactly-weighted aggregate
            Xlo, Xhi, Ulo, Uhi = pfqn_scb(Dsc, R)
            self.assertLessEqual(Xlo, XR + 1e-9, msg='K=%d R=%d' % (K, R))
            self.assertGreaterEqual(Xhi, XR - 1e-9, msg='K=%d R=%d' % (K, R))
            self.assertLessEqual((XR - Xlo) / XR, pfqn_scbgap(R, K) + 1e-9)


if __name__ == '__main__':
    unittest.main()
