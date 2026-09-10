"""
Closed-form oracles for the api/wf phase-type algebra.

Python is the reference for api/wf (MATLAB has no counterpart), so these
oracles pin the reference itself rather than a port. Every convolution has a
closed-form mean its result must match, and none of those means depends on the
representation the formula happens to build:

  - sequence: the means add;
  - parallel: the maximum of two independent exponentials has mean
    1/l1 + 1/l2 - 1/(l1+l2);
  - loop: a geometric number of repetitions has mean m/(1-p);
  - branch: the means mix with the branch probabilities.

The mean of an (alpha, T) law is alpha (-T)^-1 e, computed here by a solve, so
the test never reuses the assembly it is checking. The same oracles are pinned
in jar/src/test/java/jline/api/wf/WfPatternUpdaterTest.java and in
cpp/tests/test_wf_pattern_updater.cpp.
"""

import os
import sys
import unittest

import numpy as np

# Add the python root for direct execution
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from line_solver.api.wf import (
    ServiceParameters,
    update_patterns,
    get_update_stats,
    validate_sequence,
)
from line_solver.api.wf.pattern_updater import (
    _convolve_sequence,
    _convolve_parallel,
    _convolve_loop,
    _convolve_branches,
    _remove_matrix_rows,
)

TOL = 1e-10


def expo(rate):
    """An exponential of the given rate as an (alpha, T) pair."""
    return ServiceParameters(np.array([1.0]), np.array([[-rate]]))


def erlang(rate, k):
    """Erlang-k of total mean k/rate."""
    alpha = np.zeros(k)
    alpha[0] = 1.0
    T = np.zeros((k, k))
    for i in range(k):
        T[i, i] = -rate
        if i + 1 < k:
            T[i, i + 1] = rate
    return ServiceParameters(alpha, T)


def ph_mean(p):
    """alpha (-T)^-1 e, by a solve rather than by the assembly under test."""
    alpha = np.asarray(p.alpha).flatten()
    x = np.linalg.solve(-np.asarray(p.T), np.ones(len(alpha)))
    return float(np.dot(alpha, x))


class TestWfPatternUpdater(unittest.TestCase):

    def test_convolve_sequence_adds_the_means(self):
        c = _convolve_sequence([expo(2.0), erlang(4.0, 3), expo(5.0)])
        self.assertEqual(len(np.asarray(c.alpha).flatten()), 5)
        self.assertAlmostEqual(ph_mean(c), 0.5 + 0.75 + 0.2, delta=TOL)

        self.assertEqual(len(np.asarray(_convolve_sequence([]).alpha).flatten()), 1)
        self.assertAlmostEqual(ph_mean(_convolve_sequence([expo(3.0)])), 1.0 / 3.0, delta=TOL)

    def test_convolve_parallel_gives_the_maximum_of_two_exponentials(self):
        l1, l2 = 2.0, 3.0
        c = _convolve_parallel([expo(l1), expo(l2)])
        self.assertEqual(len(np.asarray(c.alpha).flatten()), 3)
        self.assertAlmostEqual(ph_mean(c), 1.0 / l1 + 1.0 / l2 - 1.0 / (l1 + l2), delta=TOL)

    def test_convolve_parallel_dominates_each_branch(self):
        m1 = ph_mean(_convolve_parallel([erlang(4.0, 2), expo(1.0)]))
        m2 = ph_mean(_convolve_parallel([expo(1.0), erlang(4.0, 2)]))
        self.assertAlmostEqual(m1, m2, delta=TOL)
        self.assertGreater(m1, 1.0)  # the maximum exceeds the slower branch's own mean

    def test_convolve_loop_inflates_the_mean_by_the_geometric_factor(self):
        self.assertAlmostEqual(ph_mean(_convolve_loop(erlang(4.0, 2), 0.4)),
                               0.5 / (1.0 - 0.4), delta=TOL)
        # outside (0,1) the law is returned unchanged
        self.assertAlmostEqual(ph_mean(_convolve_loop(erlang(4.0, 2), 0.0)), 0.5, delta=TOL)
        self.assertAlmostEqual(ph_mean(_convolve_loop(erlang(4.0, 2), 1.0)), 0.5, delta=TOL)

    def test_convolve_branches_mixes_the_means(self):
        ps = [expo(2.0), erlang(4.0, 4)]
        c = _convolve_branches(ps, [0.25, 0.75])
        self.assertEqual(len(np.asarray(c.alpha).flatten()), 5)
        self.assertAlmostEqual(ph_mean(c), 0.25 * 0.5 + 0.75 * 1.0, delta=TOL)

        # unnormalized probabilities are renormalized, so the answer is the same
        self.assertAlmostEqual(ph_mean(_convolve_branches(ps, [1.0, 3.0])), ph_mean(c), delta=TOL)
        # a zero total falls back to a uniform choice
        self.assertAlmostEqual(ph_mean(_convolve_branches(ps, [0.0, 0.0])), 0.75, delta=TOL)

    def test_remove_matrix_rows_drops_every_listed_row(self):
        m = np.arange(5).reshape(5, 1).astype(float)
        r = _remove_matrix_rows(m, [1, 3])
        self.assertEqual(r.shape[0], 3)
        self.assertTrue(np.array_equal(r.flatten(), np.array([0.0, 2.0, 4.0])))

    def test_update_patterns_collapses_a_chain(self):
        # 1 -> 2 -> 3 -> 4 with 2,3,4 service nodes: the sequence 2-3-4 collapses
        # onto node 2, whose law becomes the convolution of the three.
        link = np.array([[1, 2, 1.0], [2, 3, 1.0], [3, 4, 1.0]])
        params = {2: expo(2.0), 3: expo(4.0), 4: expo(5.0)}

        w = update_patterns(link, [2, 3, 4], [], [], [], params)
        self.assertLess(w.link_matrix.shape[0], link.shape[0])
        self.assertIn(2, w.service_parameters)
        self.assertNotIn(3, w.service_parameters)
        self.assertNotIn(4, w.service_parameters)
        self.assertAlmostEqual(ph_mean(w.service_parameters[2]), 0.95, delta=TOL)

        stats = get_update_stats(link, w)
        self.assertEqual(stats['originalLinks'], 3)

    def test_update_patterns_is_the_identity_when_nothing_is_detected(self):
        link = np.array([[1, 2, 1.0]])
        params = {2: expo(1.0)}
        w = update_patterns(link, [2], [], [], [], params)
        self.assertEqual(w.link_matrix.shape[0], 1)
        self.assertEqual(len(w.service_parameters), 1)
        self.assertAlmostEqual(ph_mean(w.service_parameters[2]), 1.0, delta=TOL)

    def test_validate_sequence(self):
        link = np.array([[0, 1, 1.0], [1, 2, 1.0]])
        self.assertTrue(validate_sequence([0, 1, 2], link))
        self.assertFalse(validate_sequence([0, 2], link))


if __name__ == '__main__':
    unittest.main()
