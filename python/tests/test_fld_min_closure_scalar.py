"""The scalar kernels of the fluid closures answer what their array paths answer.

`closures.min_closure` and `closures.share_closure` are evaluated once per
station per right-hand side by the fluid ODE, always on ONE point: a scalar
population for the first, one coordinate per class at a single station for the
second. Both therefore carry a scalar kernel that skips the array machinery,
and this file is what keeps the two implementations one function.

WHY THE ASSERTIONS ARE NOT BITWISE. Neither kernel can be, and for a reason that
is numpy's rather than this code's:

  * `np.exp` over an ARRAY uses a SIMD kernel and libm evaluates the scalar, and
    the two disagree in the last ulp on ~5% of arguments. The array path is not
    self-consistent either -- its answer depends on the length of the array it
    is handed -- so there is no bitwise reference to be identical to.
  * `add.reduce` and matmul accumulate in an unrolled/pairwise ORDER chosen from
    the length, where a scalar kernel sums left to right.

So each kernel is held to the bound that was measured for it, several orders
below the 1e-4 the fluid ODE is integrated to, plus the properties the callers
actually READ: the branch a value selects, and the sign of a share.
"""
import numpy as np
import pytest

from line_solver.solvers.solver_fld.utils import closures as CL


def _share_array_path(x, wv, C):
    """`share_closure`'s array path, verbatim, for the no-Jacobian case."""
    x = np.asarray(x, dtype=float).ravel()
    wv = np.asarray(wv, dtype=float).ravel()
    n = len(x)
    u = wv * x
    v = float(np.sum(u))
    if v <= 0:
        return np.zeros(n)
    s = u / v
    if C is None or not np.any(C):
        return s
    C = np.asarray(C, dtype=float)
    cuv = wv * (C @ wv)
    cvv = float(wv @ C @ wv)
    # the correction is admitted only where the series converges, exactly as the
    # array path does it; outside that region the share stays first-order
    tau, _ = CL._expansion_weight(cvv / v ** 2)
    if tau <= 0:
        return s
    s = s + tau * (-cuv / v ** 2 + (u * cvv) / v ** 3)
    if np.all(s >= -CL.GlobalConstants.Zero):
        return np.where(s < 0, 0.0, s)
    act = s > 0
    if not np.any(act):
        return u / v
    T = float(np.sum(s[act]))
    snew = np.zeros(n)
    snew[act] = s[act] / T
    return snew


def _min_cases(rng, count):
    for _ in range(count):
        n = float(rng.choice([rng.uniform(0, 200), rng.uniform(0, 1e-9), 0.0, 1.0]))
        c = float(rng.choice([rng.uniform(0.5, 10), 1.0, np.inf, n, n + 1e-15]))
        s2 = float(rng.choice([0.0, rng.uniform(0, 50), 1e-18]))
        vc = float(rng.choice([0.0, rng.uniform(0, 5)]))
        cov = float(rng.choice([0.0, rng.uniform(-5, 5)]))
        yield n, c, s2, vc, cov


def test_min_closure_scalar_matches_the_vector_path():
    rng = np.random.default_rng(101)
    seen_degenerate = seen_smooth = 0
    for n, c, s2, vc, cov in _min_cases(rng, 4000):
        hs, ds, _ = CL._min_closure_scalar(n, c, s2, vc, cov)
        hv, dv = CL.min_closure(np.array([n]), np.array([c]), np.array([s2]),
                                np.array([vc]), np.array([cov]))
        hv, dv = float(hv[0]), float(dv[0])

        th2 = max(s2 - 2 * cov + vc, 0.0)
        degenerate = (not (th2 > 0)) or np.isinf(c)
        if degenerate:
            seen_degenerate += 1
            # no transcendental on this branch, so it IS exact
            assert hs == hv
        else:
            seen_smooth += 1
            assert abs(hs - hv) <= 1e-12 * max(1.0, abs(hv))
        # the derivative is 1 - Phi and never cancels, on either branch
        assert ds == dv
    assert seen_degenerate > 100 and seen_smooth > 100


def test_min_closure_scalar_is_reached_through_the_public_entry():
    """Scalar arguments must take the kernel; anything array-shaped must not."""
    h, dh = CL.min_closure(3.0, 2.0, 1.5)
    assert h.shape == (1,) and dh.shape == (1,)
    assert (h[0], dh[0]) == CL._min_closure_scalar(3.0, 2.0, 1.5, 0.0, 0.0)[:2]
    # a vector argument still broadcasts against the scalars, as before
    hv, dhv = CL.min_closure(np.array([1.0, 3.0, 9.0]), 2.0, 0.0)
    assert hv.shape == (3,)
    assert list(hv) == [1.0, 2.0, 2.0]


def _share_cases(rng, count):
    for _ in range(count):
        n = int(rng.integers(1, CL._SHARE_SCALAR_MAX + 1))
        x = rng.uniform(0, 100, n)
        if rng.random() < 0.15:
            x = np.zeros(n)
        elif rng.random() < 0.2:
            x[rng.integers(0, n)] = 0.0
        wv = rng.uniform(0, 3, n) if rng.random() < 0.7 else np.ones(n)
        if rng.random() < 0.1:
            wv = np.zeros(n)
        if rng.random() < 0.15:
            C = None
        else:
            A = rng.normal(0, rng.choice([0.3, 30.0, 3000.0]), (n, n))
            C = A @ A.T
            if rng.random() < 0.1:
                C = np.zeros((n, n))
        yield x, wv, C


def test_share_closure_scalar_matches_the_array_path():
    rng = np.random.default_rng(102)
    worst = 0.0
    for x, wv, C in _share_cases(rng, 4000):
        a = np.asarray(_share_array_path(x, wv, C), dtype=float)
        b = np.asarray(CL.share_closure(x, wv, C), dtype=float)
        assert a.shape == b.shape
        worst = max(worst, float(np.max(np.abs(a - b))))
        # the branches below the correction read only the SIGN of a share, so a
        # last-ulp disagreement must never move one across zero
        assert ((a > 0) == (b > 0)).all()
    # measured 2.8e-14 over a 60k sweep; the bound is that with headroom
    assert worst <= 1e-12, worst


def test_share_closure_scalar_conserves_capacity():
    """The invariant the closure exists to preserve: the shares sum to one.

    Stated on `share_closure` itself -- ``sum_j Cov(u_j,v) = Var(v)`` makes the
    two correction terms cancel in the sum -- and it is what a work-conserving
    discipline requires, so it must survive the scalar kernel.

    The cancellation is not exact in floating point: with a large covariance the
    correction is large and the two terms that cancel are each far bigger than
    their difference, which costs digits. The array path loses the SAME digits on
    the SAME inputs (both reach 2.9e-11 on this sweep), so what is asserted is
    that the kernel does not lose MORE -- an absolute bound alone would be a
    statement about the closure rather than about this port.
    """
    rng = np.random.default_rng(103)
    for x, wv, C in _share_cases(rng, 1500):
        s = np.asarray(CL.share_closure(x, wv, C), dtype=float)
        if float(np.sum(np.asarray(wv) * np.asarray(x))) <= 0:
            assert not s.any()
            continue
        ref = np.asarray(_share_array_path(x, wv, C), dtype=float)
        err = abs(float(s.sum()) - 1.0)
        assert err <= max(1e-14, 2.0 * abs(float(ref.sum()) - 1.0))


def test_share_closure_keeps_the_array_path_for_a_jacobian():
    """`want_jac` has no scalar kernel and must still return both outputs."""
    x = np.array([3.0, 5.0])
    wv = np.ones(2)
    C = np.array([[2.0, 0.5], [0.5, 1.0]])
    s, ds = CL.share_closure(x, wv, C, want_jac=True)
    assert s.shape == (2,) and ds.shape == (2, 2)
    assert np.allclose(s, np.asarray(CL.share_closure(x, wv, C), dtype=float))


@pytest.mark.parametrize('n', [1, 2, CL._SHARE_SCALAR_MAX])
def test_share_closure_without_covariance_is_the_plug_in_ratio(n):
    x = np.arange(1.0, n + 1.0)
    wv = np.ones(n)
    s = np.asarray(CL.share_closure(x, wv, None), dtype=float)
    assert np.array_equal(s, x / x.sum())
