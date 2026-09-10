"""Tests for pfqn_stdf, the exact sojourn time law at FCFS stations.

Every expected value here is a CLOSED FORM computed in the test, not a number
transcribed from MATLAB. That matters because pfqn_stdf had no caller in the
python tree and so was never exercised: a transcribed golden could have been
copied from a run of the broken code, whereas an analytic oracle cannot be.

Model throughout: Delay Exp(1) -> FCFS Exp(2), closed, single class, N jobs,
folded for pfqn_stdf as one queueing station with L=0.5, Z=1, rate=2.

By the arrival theorem the tagged job sees the N-1 network, so
P(n others at the queue) is proportional to D^n Z^(N-1-n)/(N-1-n)! and the
sojourn is Erlang(n+1, mu). The CDF is that Erlang mixture.
"""

import unittest
from math import factorial

import numpy as np

from line_solver.api.pfqn.stdf import pfqn_stdf
from line_solver.api.pfqn.mva import pfqn_mva

D, Z, MU = 0.5, 1.0, 2.0


def _stdf_cdf(N, tset, S=1, L=None, rates=None, Zv=None, k=0):
    L = np.array([[D]]) if L is None else L
    rates = np.array([[MU]]) if rates is None else rates
    Zv = np.array([Z]) if Zv is None else Zv
    S = np.atleast_1d(np.asarray(S, dtype=int))
    RD = pfqn_stdf(L, np.array([N]), Zv, S, np.array([k]), rates, tset)
    return RD[(k, 0)][:, 0]


def _oracle_cdf(N, t):
    """Exact Erlang(n+1, MU) mixture seen by the tagged job."""
    w = np.array([D ** n * Z ** (N - 1 - n) / factorial(N - 1 - n) for n in range(N)])
    w /= w.sum()
    t = np.atleast_1d(np.asarray(t, dtype=float))
    F = np.zeros_like(t)
    for n, wn in enumerate(w):
        s = np.zeros_like(t)
        for j in range(n + 1):
            s += np.exp(-MU * t) * (MU * t) ** j / factorial(j)
        F += wn * (1.0 - s)
    return F


class TestPfqnStdf(unittest.TestCase):

    def test_single_job_is_exactly_the_service_time(self):
        """N=1: a lone job never queues, so its sojourn IS its service time.

        F(t) must be exactly 1 - exp(-mu t). This is the sharpest localiser in
        the suite: N=1 exercises hkc, the empty-station-set constant lGk, the
        lGr normalization and the Hkrt sum on the smallest possible state
        space, so a failure here cannot be blamed on an oracle.
        """
        t = np.linspace(0.05, 6.0, 40)
        got = _stdf_cdf(1, t)
        want = 1.0 - np.exp(-MU * t)
        self.assertLess(float(np.max(np.abs(got - want))), 1e-10)

    def test_matches_erlang_mixture_oracle(self):
        """Pointwise against the exact law, N = 1..6."""
        t = np.linspace(0.05, 8.0, 32)
        for N in range(1, 7):
            got = _stdf_cdf(N, t)
            want = _oracle_cdf(N, t)
            self.assertLess(float(np.max(np.abs(got - want))), 1e-9,
                            "N=%d sojourn CDF deviates from the exact law" % N)

    def test_is_a_proper_distribution(self):
        """Monotone non-decreasing, in [0,1], and -> 1.

        The original defect returned F(t) = 1.0 at EVERY t, which is why a
        "reaches 1" check alone is not enough and monotonicity is asserted
        together with a value strictly below 1 at small t.
        """
        t = np.linspace(0.01, 40.0, 200)
        for N in (1, 3, 5):
            F = _stdf_cdf(N, t)
            self.assertTrue(np.all(np.diff(F) >= -1e-12), "N=%d CDF decreases" % N)
            self.assertTrue(np.all(F >= -1e-12) and np.all(F <= 1.0 + 1e-12))
            self.assertLess(F[0], 0.5, "N=%d CDF is already ~1 at t=0.01" % N)
            self.assertGreater(F[-1], 1.0 - 1e-6, "N=%d CDF does not reach 1" % N)

    def test_mean_equals_exact_mva_residence_time(self):
        """E[T] = int_0^inf (1-F) dt must equal the exact MVA residence time.

        Gauss-Legendre, because pfqn_stdf solves a normalizing constant at
        every time point and a fine grid is prohibitively slow.
        """
        tmax, ng = 60.0, 96
        x, w = np.polynomial.legendre.leggauss(ng)
        t = 0.5 * tmax * (x + 1.0)
        wt = 0.5 * tmax * w
        for N in (2, 3, 5):
            F = _stdf_cdf(N, t)
            mean = float(np.dot(wt, 1.0 - F))
            out = pfqn_mva(np.array([[D]]), np.array([N]), np.array([Z]))
            # R_k = Q_k / X by Little's law: pfqn_mva's output ORDER differs
            # between MATLAB and python, so indexing it by position is a trap.
            X = float(np.ravel(np.asarray(out[0]))[0])
            Rk = float(np.ravel(np.asarray(out[2]))[0]) / X
            self.assertAlmostEqual(mean, Rk, places=6,
                                   msg="N=%d mean sojourn != MVA residence time" % N)

    def test_multiserver_mean_matches_ld_mva(self):
        """The multiserver path (S>1) goes through pfqn_mvald, not comomrm_ld."""
        from line_solver.api.pfqn.mvald import pfqn_mvald
        tmax, ng = 60.0, 96
        x, w = np.polynomial.legendre.leggauss(ng)
        t = 0.5 * tmax * (x + 1.0)
        wt = 0.5 * tmax * w
        for N, S in ((4, 2), (6, 3)):
            F = _stdf_cdf(N, t, S=S)
            mean = float(np.dot(wt, 1.0 - F))
            mu = np.zeros((1, N))
            for n in range(N):
                mu[0, n] = min(S, n + 1)
            out = pfqn_mvald(np.array([[D]]), np.array([N]), np.array([Z]), mu)
            X = float(np.ravel(np.asarray(out[0]))[0])
            Rk = float(np.ravel(np.asarray(out[1]))[0]) / X
            self.assertAlmostEqual(mean, Rk, places=6,
                                   msg="N=%d S=%d multiserver mean mismatch" % (N, S))


if __name__ == '__main__':
    unittest.main()
