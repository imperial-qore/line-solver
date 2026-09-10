"""
lossn_rec: the exact normalising constant of a loss network by MDD-rec, and the
SolverNC route that selects it.

WHY THIS METHOD EXISTS. The admissible set of a Kelly loss network is
{n >= 0 : A n <= C} and the stationary law is independent Poisson counts
truncated to it, so the normalising constant is a sum of a product form over a
set that is finite and bounded per coordinate -- exactly what a decision diagram
holds and mdd_rec walks. The route it replaces is the Manjunath-Sikdar residue
transform, which is equally exact but whose residue argument COUNTS WHOLE UNITS
and so needs an integral A and C. On a fractional region the analyzer used to
fall back to the Erlang fixed point, an approximation; these tests pin the size
of the error that fallback was making and show that MDD-rec removes it.

THREE ORACLES:
 1. BRUTE FORCE. G, the carried load and the blocking probabilities summed
    state by state over the admissible set, sharing no code path with the walk.
 2. lossn_manjunath. On an INTEGRAL region the two exact methods must agree.
 3. THE OTHER CODEBASES. Pinned at 12 decimals.
"""

from itertools import product
from math import factorial, log

import numpy as np
import pytest

from line_solver import Delay, Exp, Network, OpenClass, Sink, Source, SolverNC
from line_solver.api.lossn import lossn_erlangfp, lossn_manjunath, lossn_rec

# 2 links, 3 routes: routes 0 and 1 take one link each, route 2 takes both.
A_INT = np.array([[1.0, 0.0, 1.0], [0.0, 1.0, 1.0]])
C_INT = np.array([6.0, 5.0])
NU_INT = np.array([2.5, 1.8, 1.2])

# One link with FRACTIONAL class sizes and capacity: the residue transform
# cannot count these, and erlangfp is what the default used to fall back to.
A_FRAC = np.array([[1.5, 0.75, 2.25]])
C_FRAC = np.array([7.5])
NU_FRAC = np.array([2.0, 3.0, 1.0])

# MATLAB, the JAR and the C++ port, at %.12f.
INT_QLEN = [2.282838019229, 1.619031944010, 0.994469884949]
INT_LOSS = [0.086864792309, 0.100537808883, 0.171275095875]
INT_LG = 5.342866909910
FRAC_QLEN = [1.383844183191, 2.545477546770, 0.537948865597]
FRAC_LOSS = [0.308077908404, 0.151507484410, 0.462051134403]
FRAC_LG = 5.451139581448


def _brute(nu, A, C, cap=14):
    """G, the carried load and the blocking, summed state by state."""
    K = nu.size
    tot = 0.0
    num = np.zeros(K)
    pts = []
    for n in product(*[range(cap)] * K):
        n = np.asarray(n, dtype=float)
        if np.all(A @ n <= C + 1e-12):
            w = float(np.prod([nu[r] ** n[r] / factorial(int(n[r])) for r in range(K)]))
            tot += w
            num += w * n
            pts.append((n, w))
    acc = np.array([sum(w for n, w in pts if np.all(A @ (n + np.eye(K)[r]) <= C + 1e-12)) / tot
                    for r in range(K)])
    return num / tot, 1.0 - acc, log(tot)


def test_matches_the_brute_force_sum_on_an_integral_region():
    q, l, lg, it = lossn_rec(NU_INT, A_INT, C_INT)
    bq, bl, blg = _brute(NU_INT, A_INT, C_INT)
    np.testing.assert_allclose(q, bq, atol=1e-12)
    np.testing.assert_allclose(l, bl, atol=1e-12)
    assert lg == pytest.approx(blg, rel=1e-12)
    assert it == NU_INT.size + 1          # one walk for G, one per class for blocking


def test_agrees_with_the_residue_transform_where_both_apply():
    q, l, lg, _ = lossn_rec(NU_INT, A_INT, C_INT)
    mq, ml, mlg, _ = lossn_manjunath(NU_INT, A_INT, C_INT)
    np.testing.assert_allclose(q, mq, atol=1e-11)
    np.testing.assert_allclose(l, ml, atol=1e-11)
    assert lg == pytest.approx(mlg, rel=1e-11)


def test_is_exact_where_the_residue_transform_cannot_count():
    q, l, lg, _ = lossn_rec(NU_FRAC, A_FRAC, C_FRAC)
    bq, bl, blg = _brute(NU_FRAC, A_FRAC, C_FRAC)
    np.testing.assert_allclose(q, bq, atol=1e-12)
    np.testing.assert_allclose(l, bl, atol=1e-12)
    assert lg == pytest.approx(blg, rel=1e-12)


def test_the_erlang_fallback_it_replaces_was_not_exact():
    # The point of the method: this is the error the old default was making on a
    # fractional region. Assert it is REAL, so that a future change silently
    # routing back to erlangfp cannot pass unnoticed.
    _, bl, _ = _brute(NU_FRAC, A_FRAC, C_FRAC)
    _, el, _, _ = lossn_erlangfp(NU_FRAC, A_FRAC, C_FRAC)
    assert np.max(np.abs(np.asarray(el) - bl)) > 1e-3


def test_values_are_pinned_across_the_codebases():
    q, l, lg, _ = lossn_rec(NU_INT, A_INT, C_INT)
    np.testing.assert_allclose(q, INT_QLEN, rtol=1e-11)
    np.testing.assert_allclose(l, INT_LOSS, rtol=1e-11)
    assert lg == pytest.approx(INT_LG, rel=1e-11)
    q, l, lg, _ = lossn_rec(NU_FRAC, A_FRAC, C_FRAC)
    np.testing.assert_allclose(q, FRAC_QLEN, rtol=1e-11)
    np.testing.assert_allclose(l, FRAC_LOSS, rtol=1e-11)
    assert lg == pytest.approx(FRAC_LG, rel=1e-11)


def test_a_class_consuming_nothing_is_refused_by_name():
    with pytest.raises(Exception, match='consumes no resource'):
        lossn_rec(np.array([1.0, 1.0]), np.array([[1.0, 0.0]]), np.array([3.0]))


def _fcr_model(memsize=None):
    model = Network('FCR Loss Network')
    source, delay, sink = Source(model, 'Source'), Delay(model, 'Delay'), Sink(model, 'Sink')
    c1, c2 = OpenClass(model, 'Class1', 0), OpenClass(model, 'Class2', 1)
    source.set_arrival(c1, Exp.fit_rate(0.3))
    source.set_arrival(c2, Exp.fit_rate(0.2))
    delay.set_service(c1, Exp.fit_rate(1.0))
    delay.set_service(c2, Exp.fit_rate(0.8))
    P = model.init_routing_matrix()
    P.set(c1, c1, source, delay, 1.0)
    P.set(c1, c1, delay, sink, 1.0)
    P.set(c2, c2, source, delay, 1.0)
    P.set(c2, c2, delay, sink, 1.0)
    model.link(P)
    fcr = model.add_region(delay)
    fcr.set_global_max_jobs(5)
    fcr.set_class_max_jobs(c1, 3)
    fcr.set_class_max_jobs(c2, 3)
    if memsize is not None:
        fcr.set_global_max_memory(memsize[0])
        fcr.set_class_size(c1, memsize[1])
        fcr.set_class_size(c2, memsize[2])
    fcr.set_drop_rule(c1, True)
    fcr.set_drop_rule(c2, True)
    return model


def test_solver_nc_takes_rec_on_a_fractional_region_and_exact_on_an_integral_one():
    integral = SolverNC(_fcr_model())
    integral.get_avg_qlen()
    assert integral._result.method == 'lossn.exact'

    fractional = SolverNC(_fcr_model(memsize=(7.5, 1.5, 0.75)))
    fractional.get_avg_qlen()
    assert fractional._result.method == 'lossn.rec'

    # and rec must agree with the transform wherever the transform can go
    asked = SolverNC(_fcr_model(), method='rec')
    np.testing.assert_allclose(np.ravel(asked.get_avg_qlen()),
                               np.ravel(integral.get_avg_qlen()), atol=1e-11)
