"""Regression tests for the load-dependent QRF no-blocking polytope.

Ground truth is glpsol 5.0 on GMPL instances of the paper's
qrboundsrsrd_skel.mod. On a 2-station cyclic exponential network the QRF
polytope is tight: min and max of the station-1 utilization coincide and equal
the exact product-form value, so a single number checks both sides of the
bound and any slackening shows up immediately.

    N=3, mu2=1, cyclic routing
      mu1=0.5              -> 14/15   = 0.933333333
      mu1=1                ->  3/4    = 0.750000000
      mu1=4                -> 21/85   = 0.247058824
      mu1=1, alpha1=[1,2,3]->          0.625000000
      mu1=1, alpha1=[1,3,9]->          0.578125000

Three defects are pinned, each of which independently collapsed the polytope
to the vacuous [0, 1] for every case above:

  * q arity. qrboundsrsrd_skel.mod:11 declares q with FIVE indices, the fifth
    being the population of the EMITTING station, and folds alpha[i,n] into
    BOTH the background term v and the completion term. A population-free q
    cannot express a load-dependent rate at all.
  * THM30's right-hand side is p2[j,nj,h, i,0+1,k, m] with q[i,j,k,u,1], i.e.
    population 1. Writing the population-0 entry instead is an index-base slip
    that forces the boundary term to zero.
  * THM30 and THM3 must be emitted at all. They are the marginal-balance
    families (THM3b and THM3a of the skeleton); without them nothing ties
    consecutive populations together.

The polytope, not an SLSQP iterate, is the primary oracle here: an empty or
vacuous polytope still yields a plausible-looking iterate, so the LP is what
decides. The end-to-end entry points are checked separately, against the same
numbers, once the polytope itself is known to be right.
"""

import importlib

import numpy as np
import pytest
from scipy.optimize import linprog

from line_solver.api.mapqn.qrf_noblo_common import (
    affine_constraint_matrices, build_q_ld, compute_num_vars,
    extract_mu_v_from_maps, reduce_equalities, sub_qrfcon_noblo)
from line_solver.api.mapqn.qrf_noblo_common import build_q_from_mu_v_rt
from line_solver.api.mapqn.qrf_noblo_mem import qrf_noblo_mem
from line_solver.api.mapqn.qrf_noblo_mmi import qrf_noblo_mmi
from line_solver.api.mapqn.qrf_noblo_mmi_ld import qrf_noblo_mmi_ld
from line_solver.api.mapqn.qrf_noblo_mmi_linear import qrf_noblo_mmi_linear

_LIN = importlib.import_module('line_solver.api.mapqn.qrf_noblo_mmi_linear')

# glpsol --math on qrboundsrsrd_skel.mod; min and max coincide in every case.
CASES = [
    ('mu1=0.5', 0.5, None, 0.933333333),
    ('mu1=1', 1.0, None, 0.750000000),
    ('mu1=4', 4.0, None, 0.247058824),
    ('alpha=[1,2,3]', 1.0, [[1., 2., 3.], [1., 1., 1.]], 0.625000000),
    ('alpha=[1,3,9]', 1.0, [[1., 3., 9.], [1., 1., 1.]], 0.578125000),
]

_N = 3
_MR = 1


def _instance(mu1, alpha):
    """2-station cyclic exponential network, station 2 at rate 1."""
    M = 2
    MAPs = [[np.array([[-mu1]]), np.array([[mu1]])],
            [np.array([[-1.0]]), np.array([[1.0]])]]
    rt = np.array([[0.0, 1.0], [1.0, 0.0]])
    K = np.array([1, 1])
    F = np.full(M, _N, dtype=int)
    mu, v = extract_mu_v_from_maps(MAPs, M, K)
    q = build_q_ld(M, K, mu, v, rt, _N, alpha)
    return M, K, F, q


def _u1_range(A_eq, b_eq, A_ub, b_ub, obj, n):
    """Minimize and maximize obj over the polytope."""
    out = []
    for sense in ('min', 'max'):
        c = obj if sense == 'min' else -obj
        res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                      bounds=[(0.0, 1.0)] * n, method='highs')
        assert res.success, '%s: %s' % (sense, res.message)
        out.append(res.fun if sense == 'min' else -res.fun)
    return out


@pytest.mark.parametrize('name,mu1,alpha,expected', CASES)
def test_linear_builder_polytope_matches_glpsol(name, mu1, alpha, expected):
    """_build_linear_constraints pins U1 to the glpsol optimum."""
    M, K, F, q = _instance(mu1, alpha)
    n = _LIN._compact_num_vars(K, _N, _MR)
    Aeq, beq, Aub, bub = _LIN._build_linear_constraints(
        q, M, _MR, np.zeros((_MR, M)), F, _N, K, n)

    obj = np.zeros(n)
    for k in range(K[0]):
        for n1 in range(1, F[0] + 1):
            for m in range(_MR):
                obj[_LIN._deltap2(0, n1, k, 0, n1, k, m, M, _N, K, _MR)] += 1.0

    lo, hi = _u1_range(Aeq, np.asarray(beq, float),
                       Aub if Aub.shape[0] else None,
                       np.asarray(bub, float) if Aub.shape[0] else None,
                       obj, n)
    assert abs(lo - expected) < 1e-8, '%s min: %.10f' % (name, lo)
    assert abs(hi - expected) < 1e-8, '%s max: %.10f' % (name, hi)


@pytest.mark.parametrize('name,mu1,alpha,expected', CASES)
def test_callback_polytope_matches_glpsol(name, mu1, alpha, expected):
    """sub_qrfcon_noblo with a 5D q pins U1 to the same optimum.

    The residuals are affine in x, so the constraint matrices are recovered
    column by column and the same LP is solved.
    """
    M, K, F, q = _instance(mu1, alpha)
    BB = np.zeros((_MR, M))
    n = compute_num_vars(M, _N, K, _MR)

    c0, e0 = sub_qrfcon_noblo(np.zeros(n), q, M, _MR, BB, F, _N, K)
    A_eq = np.zeros((len(e0), n))
    A_ub = np.zeros((len(c0), n))
    basis = np.zeros(n)
    for col in range(n):
        basis[col] = 1.0
        ci, ei = sub_qrfcon_noblo(basis, q, M, _MR, BB, F, _N, K)
        A_eq[:, col] = ei - e0
        A_ub[:, col] = ci - c0
        basis[col] = 0.0

    Kmax = int(max(K))
    obj = np.zeros(n)
    for k in range(K[0]):
        for n1 in range(1, F[0] + 1):
            for m in range(_MR):
                idx = ((((((0 * (_N + 1) + n1) * Kmax + k) * M + 0)
                         * (_N + 1) + n1) * Kmax + k) * _MR + m)
                obj[idx] += 1.0

    lo, hi = _u1_range(A_eq, -e0, A_ub if len(c0) else None,
                       -c0 if len(c0) else None, obj, n)
    assert abs(lo - expected) < 1e-8, '%s min: %.10f' % (name, lo)
    assert abs(hi - expected) < 1e-8, '%s max: %.10f' % (name, hi)


def test_build_q_ld_arity_zero_population_and_alpha_on_v():
    """q is 5D, vanishes at population 0, and scales BOTH terms by alpha."""
    M, N = 1, 3
    K = np.array([2])
    mu = np.zeros((1, 2, 2)); mu[0, 0, 1] = 5.0
    v = np.zeros((1, 2, 2)); v[0, 1, 0] = 7.0
    rt = np.array([[1.0]])
    alpha = np.array([[1.0, 2.0, 4.0]])
    q = build_q_ld(M, K, mu, v, rt, N, alpha)

    assert q.shape == (1, 1, 2, 2, N + 1)
    assert np.all(q[..., 0] == 0.0)
    # self-loop entry: v*alpha + r*mu*alpha, so alpha multiplies v as well
    for n in (1, 2, 3):
        assert q[0, 0, 0, 1, n] == pytest.approx(5.0 * alpha[0, n - 1])
        assert q[0, 0, 1, 0, n] == pytest.approx(7.0 * alpha[0, n - 1])
    # default alpha is all ones
    q1 = build_q_ld(M, K, mu, v, rt, N)
    assert q1[0, 0, 1, 0, 2] == pytest.approx(7.0)


def test_equality_block_is_over_determined_and_reduces_to_full_rank():
    """The raw equality block has more rows than variables; reduction fixes it.

    This is the precondition that keeps SLSQP alive. scipy's SLSQP cannot
    accept meq > n (its exit mode 2 is "More equality constraints than
    independent variables"), and scipy 1.16.3 does not refuse such a problem
    cleanly: its Fortran workspace is sized by a closed form valid only for
    meq <= n, so it under-allocates and the process dies with "double free or
    corruption". Asserting the reduced block satisfies meq <= n fails loudly if
    the reduction is ever dropped, instead of aborting the test process.
    """
    M, K, F, q = _instance(1.0, None)
    n = _LIN._compact_num_vars(K, _N, _MR)
    Aeq, beq, _, _ = _LIN._build_linear_constraints(
        q, M, _MR, np.zeros((_MR, M)), F, _N, K, n)
    A = Aeq.toarray()
    b = np.asarray(beq, float)
    assert A.shape[0] > n, 'fixture no longer over-determined'

    Ar, br, keep = reduce_equalities(A, b)
    assert Ar.shape[0] <= n
    assert Ar.shape[0] == np.linalg.matrix_rank(A)
    # dropped rows are implied, so the feasible set is unchanged
    x = np.linalg.lstsq(Ar, br, rcond=None)[0]
    assert np.max(np.abs(A @ x - b)) < 1e-9


def test_reduce_equalities_rejects_inconsistent_system():
    """An inconsistent system is a modelling error, not something to drop."""
    A = np.array([[1.0, 0.0], [2.0, 0.0]])
    with pytest.raises(ValueError, match='inconsistent'):
        reduce_equalities(A, np.array([1.0, 3.0]))
    # the consistent version reduces cleanly
    Ar, br, _ = reduce_equalities(A, np.array([1.0, 2.0]))
    assert Ar.shape[0] == 1


def test_affine_constraint_matrices_rejects_nonlinear_callback():
    """The matrix recovery is exact only for affine residuals; verify it."""
    A, b = affine_constraint_matrices(lambda x: np.array([x[0] + 2 * x[1] - 3]), 2)
    assert np.allclose(A, [[1.0, 2.0]]) and np.allclose(b, [3.0])
    with pytest.raises(ValueError, match='not affine'):
        affine_constraint_matrices(lambda x: np.array([float(x[0] * x[1])]), 2)


@pytest.mark.parametrize('name,mu1,alpha,expected', CASES)
def test_nlp_entry_points_run_to_completion(name, mu1, alpha, expected):
    """qrf_noblo_mmi_linear and qrf_noblo_mmi_ld reach the oracle end to end.

    Before the equality reduction these aborted the interpreter inside SLSQP.
    The polytope is tight here, so every feasible point carries the same U1 and
    the nonlinear objective cannot move the answer.
    """
    M = 2
    MAPs = [[np.array([[-mu1]]), np.array([[mu1]])],
            [np.array([[-1.0]]), np.array([[1.0]])]]
    rt = np.array([[0.0, 1.0], [1.0, 0.0]])
    a = np.asarray(alpha, dtype=float) if alpha is not None else None

    for label, fn in (('linear', qrf_noblo_mmi_linear),
                      ('ld', qrf_noblo_mmi_ld)):
        # Three values since the arms gained the alpha-weighted BN readout:
        # the mean number IN SERVICE, which this case does not exercise.
        UN, QN, _BN = fn(MAPs, _N, rt, a)
        assert np.all(np.isfinite(UN)) and np.all(np.isfinite(QN))
        assert np.all(UN <= 1.0 + 1e-8) and np.all(UN >= -1e-8)
        assert abs(UN[0] - expected) < 1e-6, '%s/%s: U1=%.9f' % (name, label, UN[0])
        assert abs(sum(QN) - _N) < 1e-6, '%s/%s: population not conserved' % (
            name, label)


# ---------------------------------------------------------------------------
# The population-free (4D q) path used by qrf_noblo_mmi / qrf_noblo_mem.
# ---------------------------------------------------------------------------

# glpsol on a GMPL transcription of qrboundsbas_skel.mod specialized to the
# no-blocking case (MR=1, no finite-capacity station f, so ZERO4/ZERO5/ZERO7/
# ZERO8/THM3f/THM3I/THM3L are vacuous and THM1old is a weaker duplicate of
# THM1). The K=2 entry is a MAP, where the bound is genuinely not tight.
_MAP_D1 = np.array([[1.016186e+00, 2.585708e-05],
                    [1.569888e-03, 1.413298e-02]])
_MAP_D0 = -np.diag(_MAP_D1.sum(axis=1))

POP_FREE_CASES = [
    ('mu1=0.5', 0.5, 0.933333333, 0.933333333),
    ('mu1=1', 1.0, 0.750000000, 0.750000000),
    ('mu1=4', 4.0, 0.247058824, 0.247058824),
]


def _pop_free_polytope(MAPs, N):
    """LP range of U1 over the 4D-q no-blocking polytope."""
    M = len(MAPs)
    K = np.array([MAPs[i][0].shape[0] for i in range(M)], dtype=int)
    MR = 1
    F = np.full(M, N, dtype=int)
    BB = np.zeros((MR, M))
    rt = np.array([[0.0, 1.0], [1.0, 0.0]])
    mu, v = extract_mu_v_from_maps(MAPs, M, K)
    q = build_q_from_mu_v_rt(M, K, mu, v, rt)
    n = compute_num_vars(M, N, K, MR)
    A_eq, b_eq = affine_constraint_matrices(
        lambda z: sub_qrfcon_noblo(z, q, M, MR, BB, F, N, K)[1], n)
    A_ub, b_ub = affine_constraint_matrices(
        lambda z: sub_qrfcon_noblo(z, q, M, MR, BB, F, N, K)[0], n)

    # sub_qrfvar consumes x in COMPACT order (phase ranges over K[i], not
    # max(K)); assuming a padded layout points the objective at the wrong
    # columns whenever the K differ.
    pos = {}
    ctr = 0
    for j in range(M):
        for nj in range(N + 1):
            for k in range(K[j]):
                for i in range(M):
                    for ni in range(N + 1):
                        for h in range(K[i]):
                            for m in range(MR):
                                pos[(j, nj, k, i, ni, h, m)] = ctr
                                ctr += 1
    obj = np.zeros(n)
    for k in range(K[0]):
        for n1 in range(1, F[0] + 1):
            obj[pos[(0, n1, k, 0, n1, k, 0)]] += 1.0

    return _u1_range(A_eq, b_eq, A_ub if A_ub.shape[0] else None,
                     b_ub if A_ub.shape[0] else None, obj, n)


@pytest.mark.parametrize('name,mu1,lo_ref,hi_ref', POP_FREE_CASES)
def test_population_free_polytope_matches_glpsol(name, mu1, lo_ref, hi_ref):
    """qrf_noblo_mmi/mem's polytope reproduces glpsol on the no-blocking model.

    Regression for the missing THM30/THM3 marginal-balance families. They and
    THM1 are the only blocks that mention q, and THM1 alone balances phases, so
    it is identically zero at K == 1: without THM30/THM3 no constraint referred
    to the service rates at all and the bound was exactly [0, 1] for every
    instance.
    """
    MAPs = [[np.array([[-mu1]]), np.array([[mu1]])],
            [np.array([[-1.0]]), np.array([[1.0]])]]
    lo, hi = _pop_free_polytope(MAPs, _N)
    assert abs(lo - lo_ref) < 1e-8, '%s min: %.10f' % (name, lo)
    assert abs(hi - hi_ref) < 1e-8, '%s max: %.10f' % (name, hi)


def test_population_free_polytope_multiphase_matches_glpsol():
    """A K=2 MAP at both stations; here the bound is genuinely not tight."""
    MAPs = [[_MAP_D0, _MAP_D1], [_MAP_D0, _MAP_D1]]
    lo, hi = _pop_free_polytope(MAPs, _N)
    assert abs(lo - 0.744064826) < 1e-8, 'min: %.10f' % lo
    assert abs(hi - 0.756030622) < 1e-8, 'max: %.10f' % hi


@pytest.mark.parametrize('mu1,expected', [
    (0.5, [0.933333333, 0.466666667]),
    (1.0, [0.750000000, 0.750000000]),
    (4.0, [0.247058824, 0.988235294]),
])
def test_noblo_mmi_and_mem_end_to_end(mu1, expected):
    """qrf_noblo_mmi and qrf_noblo_mem reach the exact utilizations.

    Before THM30/THM3 were ported these returned output that did not depend on
    the service rates at all, so the assertion that the three parameter values
    give three DIFFERENT answers is as important as the values themselves.
    """
    M, MR = 2, 1
    K = np.array([1, 1])
    rt = np.array([[0.0, 1.0], [1.0, 0.0]])
    mu = np.zeros((M, 1, 1)); mu[0, 0, 0] = mu1; mu[1, 0, 0] = 1.0
    v = np.zeros((M, 1, 1))
    MAPs = [[np.array([[-mu1]]), np.array([[mu1]])],
            [np.array([[-1.0]]), np.array([[1.0]])]]

    for label, UN, QN in (('mmi',) + qrf_noblo_mmi(M, MR, K, _N, mu, v, rt),
                          ('mem',) + qrf_noblo_mem(MAPs, _N, rt)):
        assert np.allclose(UN, expected, atol=1e-6), '%s: %s' % (label, UN)
        assert abs(sum(QN) - _N) < 1e-6, '%s: population not conserved' % label


def test_noblo_mmi_output_depends_on_service_rate():
    """Pin the exact pre-port symptom: output identical across service rates."""
    M, MR = 2, 1
    K = np.array([1, 1])
    rt = np.array([[0.0, 1.0], [1.0, 0.0]])

    def run(mu1):
        mu = np.zeros((M, 1, 1)); mu[0, 0, 0] = mu1; mu[1, 0, 0] = 1.0
        return qrf_noblo_mmi(M, MR, K, _N, mu, np.zeros((M, 1, 1)), rt)[0]

    assert not np.allclose(run(1.0), run(4.0), atol=1e-6)


def test_linear_builder_rejects_population_free_q():
    """A 4D q is the wrong model for this builder and must be refused."""
    M, K, F, _ = _instance(1.0, None)
    n = _LIN._compact_num_vars(K, _N, _MR)
    q4 = np.zeros((M, M, int(max(K)), int(max(K))))
    with pytest.raises(ValueError, match='5D'):
        _LIN._build_linear_constraints(q4, M, _MR, np.zeros((_MR, M)),
                                       F, _N, K, n)
