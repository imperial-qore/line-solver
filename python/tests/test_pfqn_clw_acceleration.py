"""
Regression tests for the two speed-ups of Choudhury, Leung and Whitt,
"Calculating normalization constants of closed queueing networks by numerically
inverting their generating functions", J. ACM 42(5):935-970, 1995, as applied by
pfqn_clw: Euler summation of the inner sums (Section 2.4) and dimension
reduction by decomposition (Section 3).

The expected values are the paper's own Tables I and III, which pin the
algorithm rather than the port; they agree across MATLAB, JAR, Python native and
C++. Table III is out of reach without the reduction: p = 11 chains cost
prod_j 2 l_j K_j contour points, and the reduction takes the inversion to
dimension 2 by inverting the hub chain and then each leaf separately.

Mirrors line-test.git/test/testsAPI/test_pfqn_clw_acceleration.m.
"""
import math
import time

import numpy as np
import pytest

from line_solver.api.pfqn.nc import pfqn_ca, pfqn_clw


def example81():
    """Example 8.1: p = 1, q' = 10 distinct queues of multiplicity 5."""
    return (0.1 * np.arange(1, 11)).reshape(-1, 1), np.array([5.0]), 5 * np.ones(10)


def star():
    """
    Examples 8.3 and 8.4: chain 1 is a hub visiting all ten queues, chain g+1
    has queue g to itself. The interdependence graph is a star, so D = {1}
    leaves ten single-variable components (the paper's Figure 1).
    """
    L = np.zeros((10, 11))
    Z = np.zeros(11)
    Z[0] = 50.0
    for g in range(1, 11):
        L[g - 1, 0] = 1 + 0.1 * g
        L[g - 1, g] = 0.1 * g
        Z[g] = 5.0 * (g + 1) - 10.0
    return L, Z


LOG10 = math.log(10.0)


def test_euler_reproduces_table_I():
    """Table I: every population from 2 to 2e7. Without Euler summation the
    last row alone would cost 4e7 contour points."""
    L, Z, m = example81()
    K = [2, 20, 200, 2000, 20000, 200000, 2000000, 20000000]
    ref = [(5.377500, 2), (1.906584, 13), (1.381312, 26), (1.284918, 31),
           (1.541538, 35), (1.569301, 39), (1.572100, 43), (1.572380, 47)]
    for k, (mant, expo) in zip(K, ref):
        _, lG = pfqn_clw(L, [k], Z, m)
        assert lG / LOG10 == pytest.approx(math.log10(mant) + expo, abs=1e-6), \
            "Table I at K1 = %d" % k


def test_euler_agrees_with_the_exact_sum():
    """The acceleration is adaptive: the Euler order is doubled until the
    paper's own estimate |E(m,n) - E(m,n+1)| settles, so switching it off must
    not move the answer. K1 = 40 is past the n+m = 31 threshold."""
    L, Z, m = example81()
    for k in [40, 100, 250]:
        _, a = pfqn_clw(L, [k], Z, m)
        _, b = pfqn_clw(L, [k], Z, m, euler=False)
        assert a == pytest.approx(b, rel=1e-9), \
            "Euler summation moved the answer at K1 = %d" % k


def test_euler_is_cheaper_than_the_exact_sum():
    """Cost is prod_j min(n+m+1, K_j) rather than prod_j K_j (eq. 2.26), so a
    population two orders larger must not cost two orders more time."""
    L, Z, m = example81()
    pfqn_clw(L, [1000], Z, m)                       # warm up
    t0 = time.time()
    pfqn_clw(L, [100000], Z, m)
    small = time.time() - t0
    t0 = time.time()
    pfqn_clw(L, [10000000], Z, m)
    big = time.time() - t0
    assert big < max(10 * small, 1.0)


def test_dimred_reproduces_table_III():
    """Table III rows 1-4: eleven chains, which only the reduction reaches."""
    L, Z = star()
    m5 = 5 * np.ones(10)
    # mantissa and exponent kept apart: 1.937826e683 is not a double
    ref = [(2, 1.235628, 25), (20, 7.503087, 45), (200, 5.970503, 129),
           (2000, 1.937826, 683)]
    for k1, mant, expo in ref:
        _, lG = pfqn_clw(L, [k1] + [2] * 10, Z, m5)
        assert lG / LOG10 == pytest.approx(math.log10(mant) + expo, abs=1e-5), \
            "Table III at K1 = %d" % k1


def test_dimred_rows_5_to_8_under_the_paper_scale_tuning():
    """Table III rows 5-8 need the manual tuning of page 956 on the hub chain,
    beta in [0.8, 1.2], which the paper prescribes for its largest examples."""
    L, Z = star()
    m5 = 5 * np.ones(10)
    Kg = list(5 * np.arange(1, 11))
    ref = [(2, 0.8, 3.004462, 107), (20, 0.8, 1.677866, 133),
           (200, 0.8, 8.032122, 260), (2000, 0.95, 1.617153, 926)]
    for k1, b1, mant, expo in ref:
        _, lG = pfqn_clw(L, [k1] + Kg, Z, m5, beta=[b1] + [1.0] * 10)
        assert lG / LOG10 == pytest.approx(math.log10(mant) + expo, abs=1e-3), \
            "Table III row for K1 = %d" % k1


def test_dimred_agrees_with_the_full_inversion():
    """Reduction reorders the inversion and splits the factors; it must not
    move the answer. Three leaves keep the undecomposed inversion affordable,
    and exact convolution adjudicates both."""
    L = np.zeros((3, 4))
    Z = np.zeros(4)
    Z[0] = 5.0
    for g in range(1, 4):
        L[g - 1, 0] = 1 + 0.1 * g
        L[g - 1, g] = 0.1 * g
        Z[g] = 5.0 * g - 5.0
    N = [4, 3, 3, 3]
    lca = pfqn_ca(L, N, Z)[1]
    lred = pfqn_clw(L, N, Z)[1]
    lfull = pfqn_clw(L, N, Z, dimred=False)[1]
    assert lred == pytest.approx(lca, rel=1e-7)
    assert lfull == pytest.approx(lca, rel=1e-6)


def test_dimred_is_inert_on_a_coupled_model():
    """Every chain visits every queue, so the interdependence graph is
    complete, no subset D reduces the dimension, and the classical path must be
    taken unchanged -- to the last bit, not merely to a tolerance."""
    L = np.array([[0.1, 0.2, 0.15], [0.3, 0.05, 0.1]])
    N = [2, 2, 1]
    Z = [1.0, 0.5, 0.2]
    a = pfqn_clw(L, N, Z)[1]
    b = pfqn_clw(L, N, Z, dimred=False)[1]
    assert a == b, "dimension reduction changed a model it cannot reduce"
    assert a == pytest.approx(-1.441982188151596, rel=1e-8)


def test_disconnected_model_needs_no_committed_variable():
    """Two chains that share no queue: the graph is already disconnected, so
    the reduction is exact with D empty and the constant is the product of the
    two single-chain constants."""
    L = np.array([[0.4, 0.0], [0.2, 0.0], [0.0, 0.5], [0.0, 0.3]])
    N = [3, 4]
    Z = [1.0, 2.0]
    both = pfqn_clw(L, N, Z)[1]
    one = pfqn_clw(L[0:2, 0:1], [3], [1.0])[1]
    two = pfqn_clw(L[2:4, 1:2], [4], [2.0])[1]
    assert both == pytest.approx(one + two, rel=1e-9)
