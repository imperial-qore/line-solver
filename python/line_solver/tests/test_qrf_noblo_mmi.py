"""Regression tests for the population-free QRF no-blocking polytope (MMI/MEM).

Ground truth is the AMPL model shipped with the QRF paper
(qrboundsbas_skel.mod, specialised to no blocking) solved with glpsol 5.0, plus
an exact oracle that does not need an LP at all.

The exact oracle. With M == 2 and N1 + N2 == N the pairwise joint p2 is fully
determined by the marginal P(N1 = n1), so the stationary distribution of the
underlying CTMC gives a point the polytope MUST contain. That is far stronger
than an LP range: an LP range only shows that SOME constraint mentioning the
service rates is present, whereas the exact point checks the sign and every
index of THM30 and THM3, which is where these families go wrong silently.

Two defect classes are pinned:

  * Missing marginal-balance families. THM30 and THM3 are the only blocks
    besides THM1 that mention q, and THM1 alone balances phases, so it is
    identically zero at K == 1. Without them no constraint refers to the
    service rates and the U bound is exactly [0, 1] for every instance, mu and
    all. This is what the JAR's buildConstraints did before the port.
  * An infeasible NLP start. From x0 = 0 (which violates ONE by a unit and COR1
    by N^2) SLSQP exits with mode 4 and returns the start point untouched, so
    UN and QN read out as identically zero with nothing in the result to say
    so. See qrf_noblo_common.feasible_start.

Comparing ROW COUNTS across implementations is misleading here: the Python
equality block is QR-reduced to an independent set (reduce_equalities) while
MATLAB and the JAR keep the raw block. Compare rank, or compare the polytope.
"""

import numpy as np
import pytest
from scipy.optimize import linprog

from line_solver.api.mapqn.qrf_noblo_common import (
    affine_constraint_matrices, build_q_from_mu_v_rt, compute_num_vars,
    extract_mu_v_from_maps, feasible_start, sub_qrfcon_noblo)
from line_solver.api.mapqn.qrf_noblo_mem import qrf_noblo_mem
from line_solver.api.mapqn.qrf_noblo_mmi import qrf_noblo_mmi

_MR = 1
_CYC2 = np.array([[0.0, 1.0], [1.0, 0.0]])

# A genuine (correlated) 2-phase MAP, from the QRF paper's instance set.
_MAP_D1 = np.array([[1.016186e+00, 2.585708e-05],
                    [1.569888e-03, 1.413298e-02]])
_MAP_D0 = -np.diag(_MAP_D1.sum(axis=1))


def _expmap(rate):
    return [np.array([[-rate]]), np.array([[rate]])]


def _erlang2(rate):
    """Erlang-2 of the given mean rate, in the layout extract_mu_v_from_maps reads.

    That helper takes mu[i,h,k] = D1[h,k] but v[i,k,h] = D0[h,k], i.e. D0
    transposed, so the phase-advance rate must sit in D0's LOWER triangle for
    the model to see a genuine Erlang-2 rather than its reverse.
    """
    D0 = np.array([[-2 * rate, 0.0], [2 * rate, -2 * rate]])
    D1 = np.array([[0.0, 0.0], [2 * rate, 0.0]])
    return [D0, D1]


def _polytope(MAPs, N):
    """(Aeq, beq, Aub, bub, n, K) for the population-free no-blocking model."""
    M = len(MAPs)
    K = np.array([MAPs[i][0].shape[0] for i in range(M)], dtype=int)
    F = np.full(M, N, dtype=int)
    BB = np.zeros((_MR, M))
    mu, v = extract_mu_v_from_maps(MAPs, M, K)
    q = build_q_from_mu_v_rt(M, K, mu, v, _CYC2)
    n = compute_num_vars(M, N, K, _MR)
    Aeq, beq = affine_constraint_matrices(
        lambda z: sub_qrfcon_noblo(z, q, M, _MR, BB, F, N, K)[1], n)
    Aub, bub = affine_constraint_matrices(
        lambda z: sub_qrfcon_noblo(z, q, M, _MR, BB, F, N, K)[0], n)
    return Aeq, beq, Aub, bub, n, K


def _compact_pos(M, N, K):
    """Column index of each p2 entry, in the order sub_qrfvar consumes x.

    The layout is COMPACT: the phase index runs over K[i], not over max(K), so
    assuming a padded layout points at the wrong columns whenever the K differ.
    """
    pos, ctr = {}, 0
    for j in range(M):
        for nj in range(N + 1):
            for k in range(K[j]):
                for i in range(M):
                    for ni in range(N + 1):
                        for h in range(K[i]):
                            for m in range(_MR):
                                pos[(j, nj, k, i, ni, h, m)] = ctr
                                ctr += 1
    return pos, ctr


def _u1_range(MAPs, N):
    Aeq, beq, Aub, bub, n, K = _polytope(MAPs, N)
    pos, _ = _compact_pos(len(MAPs), N, K)
    obj = np.zeros(n)
    for k in range(K[0]):
        for n1 in range(1, N + 1):
            obj[pos[(0, n1, k, 0, n1, k, 0)]] += 1.0
    out = []
    for s in (1.0, -1.0):
        r = linprog(s * obj, A_ub=Aub if Aub.shape[0] else None,
                    b_ub=bub if Aub.shape[0] else None, A_eq=Aeq, b_eq=beq,
                    bounds=[(0.0, 1.0)] * n, method='highs')
        assert r.success, 'LP failed: %s' % r.message
        out.append(s * r.fun)
    return out


# glpsol --math on the no-blocking specialisation of qrboundsbas_skel.mod.
GLPSOL_CASES = [
    ('exp+exp N=2', [_expmap(1.0), _expmap(1.0)], 2, 0.666666667, 0.666666667),
    ('mu1=4 N=3', [_expmap(4.0), _expmap(1.0)], 3, 0.247058824, 0.247058824),
    ('erl2+exp N=2', [_erlang2(1.0), _expmap(1.0)], 2, 0.600000000, 0.750000000),
    ('erl2+exp N=3', [_erlang2(1.0), _expmap(1.0)], 3, 0.636363636, 0.875000000),
    ('erl2+erl2 N=3', [_erlang2(1.0), _erlang2(1.0)], 3, 0.600000000, 1.000000000),
    ('MAP+MAP N=2', [[_MAP_D0, _MAP_D1], [_MAP_D0, _MAP_D1]], 2,
     0.663140224, 0.670230816),
    ('MAP+exp N=3', [[_MAP_D0, _MAP_D1], _expmap(1.0)], 3,
     0.746988870, 0.753019329),
]


@pytest.mark.parametrize('name,MAPs,N,lo_ref,hi_ref', GLPSOL_CASES)
def test_polytope_u1_range_matches_glpsol(name, MAPs, N, lo_ref, hi_ref):
    """The emitted polytope reproduces the AMPL model's U1 bound exactly.

    Heterogeneous K (the erl2+exp and MAP+exp rows) is covered deliberately:
    the decision vector is laid out compactly over K[i], so a padded index
    lands on the wrong variables only when the K differ.
    """
    lo, hi = _u1_range(MAPs, N)
    assert abs(lo - lo_ref) < 1e-8, '%s min %.10f' % (name, lo)
    assert abs(hi - hi_ref) < 1e-8, '%s max %.10f' % (name, hi)


EXACT_CASES = [(1.0, 1.0, 2), (4.0, 1.0, 3), (0.5, 1.0, 3), (1.0, 1.0, 4)]


def _exact_pairwise_joint(mu1, mu2, N):
    """The exact pairwise joint of the 2-station cyclic M/M/1//N network."""
    M, K = 2, np.array([1, 1])
    rho = mu2 / mu1
    w = np.array([rho ** n1 for n1 in range(N + 1)])
    p = w / w.sum()

    n = compute_num_vars(M, N, K, _MR)
    x = np.zeros(n)
    pos, p2_end = _compact_pos(M, N, K)
    for n1 in range(N + 1):
        n2 = N - n1
        x[pos[(0, n1, 0, 0, n1, 0, 0)]] = p[n1]
        x[pos[(1, n2, 0, 1, n2, 0, 0)]] = p[n1]
        x[pos[(0, n1, 0, 1, n2, 0, 0)]] = p[n1]
        x[pos[(1, n2, 0, 0, n1, 0, 0)]] = p[n1]
    x[p2_end] = 1.0 - p[0]        # e[0,0] = P(station 1 busy)
    x[p2_end + 1] = 1.0 - p[N]    # e[1,0] = P(station 2 busy)
    return x


@pytest.mark.parametrize('mu1,mu2,N', EXACT_CASES)
def test_exact_joint_is_feasible(mu1, mu2, N):
    """The exact CTMC pairwise joint satisfies every constraint of the polytope.

    A missing THM30/THM3 would pass this trivially, but a WRONG one would not:
    this is the check that pins their signs, their q lookups and the population
    index on their right-hand sides (p2[..., i, ni+1, ...], not ni).
    """
    MAPs = [_expmap(mu1), _expmap(mu2)]
    Aeq, beq, Aub, bub, n, _ = _polytope(MAPs, N)
    x = _exact_pairwise_joint(mu1, mu2, N)
    assert np.max(np.abs(Aeq @ x - beq)) < 1e-12
    if Aub.shape[0]:
        assert np.max(Aub @ x - bub) < 1e-12


def test_polytope_depends_on_the_service_rates():
    """Pin the pre-port symptom directly: the bound was mu-independent."""
    lo1, hi1 = _u1_range([_expmap(1.0), _expmap(1.0)], 3)
    lo4, hi4 = _u1_range([_expmap(4.0), _expmap(1.0)], 3)
    assert abs(lo1 - lo4) > 0.4 and abs(hi1 - hi4) > 0.4
    # ... and in particular is not the vacuous [0, 1].
    for lo, hi in ((lo1, hi1), (lo4, hi4)):
        assert lo > 1e-6 and hi < 1.0 - 1e-6


@pytest.mark.parametrize('mu1,expected_u', [(0.5, 0.933333333),
                                            (1.0, 0.750000000),
                                            (4.0, 0.247058824)])
def test_mmi_and_mem_end_to_end_exponential(mu1, expected_u):
    """Both objectives reach the exact utilization where the bound is tight."""
    M, N = 2, 3
    K = np.array([1, 1])
    mu = np.zeros((M, 1, 1))
    mu[0, 0, 0] = mu1
    mu[1, 0, 0] = 1.0
    v = np.zeros((M, 1, 1))
    MAPs = [_expmap(mu1), _expmap(1.0)]

    for label, UN, QN in (('mmi',) + qrf_noblo_mmi(M, _MR, K, N, mu, v, _CYC2),
                          ('mem',) + qrf_noblo_mem(MAPs, N, _CYC2)):
        assert abs(UN[0] - expected_u) < 1e-6, '%s: %s' % (label, UN)
        assert abs(sum(QN) - N) < 1e-6, '%s: population %s' % (label, QN)


def test_mmi_multiphase_returns_a_feasible_point():
    """Regression for the all-zero return on multiphase instances.

    With the old x0 = 0 start SLSQP exited with mode 4 on this instance and
    handed the start point back, so UN and QN were identically zero while the
    call reported no error at all.
    """
    MAPs = [[_MAP_D0, _MAP_D1], [_MAP_D0, _MAP_D1]]
    N = 2
    M = len(MAPs)
    K = np.array([2, 2])
    mu, v = extract_mu_v_from_maps(MAPs, M, K)
    UN, QN = qrf_noblo_mmi(M, _MR, K, N, mu, v, _CYC2)

    assert np.all(UN > 1e-6), 'utilization collapsed to zero: %s' % UN
    assert abs(sum(QN) - N) < 1e-6, 'population %s not conserved' % QN
    lo, hi = _u1_range(MAPs, N)
    assert lo - 1e-6 <= UN[0] <= hi + 1e-6, \
        'UN[0]=%.9f outside the polytope bound [%.9f, %.9f]' % (UN[0], lo, hi)


def test_feasible_start_lands_on_the_polytope():
    """The phase-1 start point satisfies the constraints it was built from."""
    MAPs = [_erlang2(1.0), _expmap(1.0)]
    Aeq, beq, Aub, bub, n, _ = _polytope(MAPs, 2)
    x0 = feasible_start(Aeq, beq, Aub, bub, n)
    assert np.max(np.abs(Aeq @ x0 - beq)) < 1e-8
    assert np.max(Aub @ x0 - bub) < 1e-8


def test_feasible_start_rejects_an_empty_polytope():
    """An empty polytope is a modelling error and must not yield a point."""
    n = 3
    Aeq = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    beq = np.array([0.0, 1.0])
    with pytest.raises(ValueError, match='infeasible'):
        feasible_start(Aeq, beq, np.zeros((0, n)), np.zeros(0), n)
