"""Regression tests for pfqn_qlen_joint_moments and the tail vertex it uses.

The routine turns normalizing constants into the joint moments of the
queue-length vector of a closed product-form network. Both of its routes are
refereed by the same independent oracle: the exact product-form state
distribution, enumerated over the whole state space and summed directly. That
oracle knows nothing about normalizing constants, survival identities or the
house of moments, so an agreement is a genuine check of the chain

    G ratios -> survival array -> binomial -> factorial -> raw -> central

rather than a round trip. The identities behind the two routes are proved
symbolically in sage/proofs/qlen_tail_moments.py.
"""

from itertools import product
from math import comb, factorial

import numpy as np
import pytest

from line_solver.api.moment import (moment_binomial_from_tail,
                                    moment_joint_binomial_from_tail,
                                    moment_joint_central_from_tail,
                                    moment_joint_tail_from_binomial,
                                    moment_tail_from_binomial)
from line_solver.api.pfqn import pfqn_nc, pfqn_qlen_joint_moments

TOL = 1e-9


def _product_form_states(L, N, Z):
    """
    Exact state distribution of a closed product-form network, by enumeration.

    Args:
        L: Demand matrix (M x R) of the queueing stations.
        N: Population vector (R,).
        Z: Think time vector (R,).

    Returns:
        Tuple (states, probs) where states[s] is an (M+1) x R array holding the
        per-class population at every queueing station and, in its last row, at
        the delay.
    """
    M, R = L.shape
    per = []
    for r in range(R):
        per.append([c for c in product(range(N[r] + 1), repeat=M + 1) if sum(c) == N[r]])
    states, weights = [], []
    for combo in product(*per):
        s = np.array([[combo[r][j] for r in range(R)] for j in range(M + 1)])
        w = 1.0
        for i in range(M):
            w *= factorial(int(s[i].sum()))
            for r in range(R):
                w *= L[i, r] ** s[i, r] / factorial(int(s[i, r]))
        for r in range(R):
            w *= Z[r] ** s[M, r] / factorial(int(s[M, r]))
        states.append(s)
        weights.append(w)
    weights = np.asarray(weights)
    return states, weights / weights.sum()


def _oracle_moments(L, N, Z, pairs):
    """
    Mean vector, covariance matrix and third central moment of the first
    coordinate, by direct summation over the exact state distribution.

    Args:
        L: Demand matrix.
        N: Population vector.
        Z: Think time vector.
        pairs: Sequence of 0-based (station, class) pairs.

    Returns:
        Tuple (mean, cov, m3).
    """
    states, probs = _product_form_states(L, N, Z)
    d = len(pairs)
    mean = np.zeros(d)
    for j, (i, r) in enumerate(pairs):
        mean[j] = sum(p * s[i, r] for p, s in zip(probs, states))
    cov = np.zeros((d, d))
    for j, (ij, rj) in enumerate(pairs):
        for l, (il, rl) in enumerate(pairs):
            cov[j, l] = sum(p * (s[ij, rj] - mean[j]) * (s[il, rl] - mean[l])
                            for p, s in zip(probs, states))
    i0, r0 = pairs[0]
    m3 = sum(p * (s[i0, r0] - mean[0]) ** 3 for p, s in zip(probs, states))
    return mean, cov, m3


def test_tail_edge_on_a_deterministic_variable():
    """N = 3 with probability one has t = (1,1,1,1,0) and b_j = C(3,j)."""
    t = [1.0, 1.0, 1.0, 1.0, 0.0]
    b = moment_binomial_from_tail(t)
    np.testing.assert_allclose(b, [1, 3, 3, 1, 0], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(moment_tail_from_binomial(b), t, rtol=1e-12, atol=1e-12)


def test_tail_edge_against_the_pmf():
    """b_j = sum_m C(m,j) p_m for an arbitrary law on a bounded support."""
    rng = np.random.default_rng(7)
    p = rng.random(6)
    p /= p.sum()
    t = np.array([p[m:].sum() for m in range(6)])
    b = moment_binomial_from_tail(t)
    for j in range(6):
        exp = sum(float(comb(m, j)) * p[m] for m in range(6))
        assert abs(b[j] - exp) < TOL, 'order %d' % j
    np.testing.assert_allclose(moment_tail_from_binomial(b), t, rtol=1e-9, atol=1e-12)


def test_joint_tail_edge_and_central_moments():
    """The joint tail edge reproduces the joint binomial moments, and the
    composed path gives the covariance of a dependent bivariate law."""
    rng = np.random.default_rng(11)
    P = rng.random((4, 4))
    P /= P.sum()
    t = np.zeros((4, 4))
    for a in np.ndindex(4, 4):
        t[a] = P[a[0]:, a[1]:].sum()
    b = moment_joint_binomial_from_tail(t)
    for a in np.ndindex(4, 4):
        exp = sum(float(comb(u, a[0])) * float(comb(v, a[1])) * P[u, v]
                  for u in range(4) for v in range(4))
        assert abs(b[a] - exp) < TOL, 'multi-order %s' % (a,)
    np.testing.assert_allclose(moment_joint_tail_from_binomial(b), t, rtol=1e-9, atol=1e-12)
    mc = moment_joint_central_from_tail(t)
    mu = [sum(u * P[u, :].sum() for u in range(4)), sum(v * P[:, v].sum() for v in range(4))]
    cov = sum(P[u, v] * (u - mu[0]) * (v - mu[1]) for u in range(4) for v in range(4))
    assert abs(mc[1, 1] - cov) < TOL


def test_single_class_tail_route_against_the_state_space():
    """Two queues, one class: the survival identity route must reproduce the
    exact mean, covariance and third central moment."""
    L = np.array([[2.0], [1.0]])
    N = np.array([5])
    Z = np.array([0.7])
    pairs = [(0, 0), (1, 0)]
    out = pfqn_qlen_joint_moments(L, N, Z, pairs=pairs)
    mean, cov, m3 = _oracle_moments(L, N, Z, pairs)
    assert out['info']['route'] == 'tail'
    np.testing.assert_allclose(out['mean'], mean, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(out['cov'], cov, rtol=1e-9, atol=1e-9)
    assert abs(out['central'][3, 0] - m3) < 1e-8
    # the survival array itself is the queue-length tail
    states, probs = _product_form_states(L, N, Z)
    exp = sum(p for p, s in zip(probs, states) if s[0, 0] >= 1 and s[1, 0] >= 2)
    assert abs(out['tail'][1, 2] - exp) < TOL


def test_multiclass_pmf_route_against_the_state_space():
    """One station, two classes: the shape a class-oriented method of moments
    produces. The naive survival formula fails here, so the complementary
    network route is the one under test."""
    L = np.array([[2.0, 1.0]])
    N = np.array([4, 3])
    Z = np.array([0.5, 0.8])
    pairs = [(0, 0), (0, 1)]
    out = pfqn_qlen_joint_moments(L, N, Z, pairs=pairs)
    mean, cov, m3 = _oracle_moments(L, N, Z, pairs)
    assert out['info']['route'] == 'pmf'
    np.testing.assert_allclose(out['mean'], mean, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(out['cov'], cov, rtol=1e-9, atol=1e-9)
    assert abs(out['central'][3, 0] - m3) < 1e-8


def test_cross_station_and_cross_class_pairs():
    """Coordinates at different stations and different classes, together."""
    L = np.array([[2.0, 1.0], [1.0, 3.0]])
    N = np.array([3, 2])
    Z = np.array([0.5, 0.0])
    for pairs in ([(0, 0), (1, 1)], [(0, 0), (0, 1), (1, 0), (1, 1)]):
        out = pfqn_qlen_joint_moments(L, N, Z, pairs=pairs)
        mean, cov, _ = _oracle_moments(L, N, Z, pairs)
        np.testing.assert_allclose(out['mean'], mean, rtol=1e-9, atol=1e-9)
        np.testing.assert_allclose(out['cov'], cov, rtol=1e-9, atol=1e-9)


def test_pmf_route_agrees_with_the_tail_route_when_both_apply():
    """With one class both routes are valid and must give the same arrays,
    although they touch different networks."""
    L = np.array([[2.0], [1.0], [0.5]])
    N = np.array([4])
    Z = np.array([0.3])
    pairs = [(0, 0), (2, 0)]
    a = pfqn_qlen_joint_moments(L, N, Z, pairs=pairs, route='tail')
    b = pfqn_qlen_joint_moments(L, N, Z, pairs=pairs, route='pmf')
    np.testing.assert_allclose(a['tail'], b['tail'], rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(a['central'], b['central'], rtol=1e-9, atol=1e-8)


def test_injected_source_is_used_and_counted():
    """A batched oracle must be called once and spare every pfqn_nc call it
    serves; the rows it declines fall back."""
    L = np.array([[2.0], [1.0]])
    N = np.array([4])
    Z = np.array([0.6])
    calls = []

    def full(Lsub, pops):
        calls.append(pops.shape[0])
        return np.array([pfqn_nc(Lsub, pops[p, :].astype(float), Z, method='ca')[1]
                         for p in range(pops.shape[0])])

    out = pfqn_qlen_joint_moments(L, N, Z, pairs=[(0, 0), (1, 0)], lg_source=full)
    assert len(calls) == 1 and calls[0] == out['info']['points']
    assert out['info']['served'] == out['info']['points']
    assert out['info']['evals'] == 0

    def partial(Lsub, pops):
        lg = np.full(pops.shape[0], np.nan)
        for p in range(pops.shape[0]):
            if pops[p, 0] % 2 == 0:
                lg[p] = pfqn_nc(Lsub, pops[p, :].astype(float), Z, method='ca')[1]
        return lg

    ref = pfqn_qlen_joint_moments(L, N, Z, pairs=[(0, 0), (1, 0)])
    got = pfqn_qlen_joint_moments(L, N, Z, pairs=[(0, 0), (1, 0)], lg_source=partial)
    assert 0 < got['info']['served'] < got['info']['points']
    assert got['info']['evals'] == got['info']['points'] - got['info']['served']
    np.testing.assert_allclose(got['central'], ref['central'], rtol=1e-9, atol=1e-9)


def test_error_paths():
    L = np.array([[2.0, 1.0]])
    N = np.array([3, 2])
    with pytest.raises(ValueError):
        pfqn_qlen_joint_moments(L, N, [0.1, 0.1], route='tail')
    with pytest.raises(ValueError):
        pfqn_qlen_joint_moments(L, N, [0.1, 0.1], pairs=[(0, 0), (0, 0)])
    with pytest.raises(ValueError):
        pfqn_qlen_joint_moments(L, N, [0.1, 0.1], pairs=[(1, 0)])
    with pytest.raises(ValueError):
        pfqn_qlen_joint_moments(L, [3], [0.1, 0.1])
    with pytest.raises(ValueError):
        pfqn_qlen_joint_moments(L, N, [0.1, 0.1], route='nosuch')
