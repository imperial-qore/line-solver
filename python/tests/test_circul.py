"""Regression tests for circul(), locked to the MATLAB and JAR reference.

circul(n) with a scalar n denotes the cyclic routing 1 -> 2 -> ... -> n -> 1.
Python previously returned the TRANSPOSE for n >= 3 (the cycle traversed
backwards), because scipy.linalg.circulant reads its argument as the first
COLUMN, not the first row. The error was invisible to every existing test:
circul(2) is symmetric, and in a closed cycle every station has visit ratio 1
regardless of order, so all product-form results are invariant under reversal.
It only shows up on non-product-form models with >= 3 stations.

Golden values below are the MATLAB circul.m output (matlab/util/circul.m),
which jline.util.Maths.circul reproduces.
"""

import numpy as np

from line_solver import circul


def test_circul_scalar_1():
    # MATLAB: circul(1) == 1 (not 0).
    np.testing.assert_array_equal(np.array(circul(1)), np.array([[1.0]]))


def test_circul_scalar_2():
    np.testing.assert_array_equal(np.array(circul(2)), np.array([[0., 1.],
                                                                 [1., 0.]]))


def test_circul_scalar_3():
    # 1 -> 2, 2 -> 3, 3 -> 1. The transpose would be the reversed cycle.
    np.testing.assert_array_equal(np.array(circul(3)), np.array([[0., 1., 0.],
                                                                 [0., 0., 1.],
                                                                 [1., 0., 0.]]))


def test_circul_scalar_4():
    np.testing.assert_array_equal(np.array(circul(4)), np.array([[0., 1., 0., 0.],
                                                                 [0., 0., 1., 0.],
                                                                 [0., 0., 0., 1.],
                                                                 [1., 0., 0., 0.]]))


def test_circul_scalar_is_forward_cycle():
    """Row i must route to i+1 (mod n) for every n, i.e. superdiagonal + wrap."""
    for n in range(2, 8):
        C = np.array(circul(n))
        for i in range(n):
            expected = np.zeros(n)
            expected[(i + 1) % n] = 1.0
            np.testing.assert_array_equal(C[i, :], expected,
                                          err_msg='circul(%d) row %d' % (n, i))


def test_circul_scalar_not_symmetric_for_n3():
    """Guard the exact failure mode: circul(3) must NOT equal its transpose."""
    C = np.array(circul(3))
    assert not np.array_equal(C, C.T)


def test_circul_vector():
    # MATLAB: circul([1,2,3]) == [1 3 2; 2 1 3; 3 2 1]. This branch was correct.
    np.testing.assert_array_equal(np.array(circul([1, 2, 3])),
                                  np.array([[1., 3., 2.],
                                            [2., 1., 3.],
                                            [3., 2., 1.]]))
