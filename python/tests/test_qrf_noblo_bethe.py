"""Regression tests for `qrf.bethe`, the tree-reweighted (Bethe) arm of the QRF
no-blocking reduction.

WHAT IS ASSERTED, AND WHY IT IS NOT A TABLE OF DIGITS
----------------------------------------------------
`qrf.bethe` minimises

    lambda * sum_{i!=j} I(n_i;n_j) - sum_i H(n_i),      lambda = 1/M

over the polytope `qrf.mmi` already uses -- the negative of a tree-reweighted
entropy with uniform edge weight rho_ij = 2*lambda on the complete station
graph. H_rho is a convex combination of tree entropies, hence concave on the
local marginal polytope, exactly when rho lies in the spanning tree polytope of
K_M, whose uniform point is rho_ij = 2/M. So lambda = 1/M is the largest
uniform weight at which minimising it is a CONVEX program, and the property
that buys is not accuracy but well-posedness: every local optimum is global, so
the answer is a property of the model rather than of where the solver started.

That is what these tests assert -- start-point independence, feasibility, and
containment in the LP range of the same polytope -- rather than frozen digits.
The reason to avoid frozen digits here is specific. Restoring the n = 0 cells
(the range the AMPL source states and the coded `mmi` does not use) brings the
structurally zero entries inside the sum, where they contribute 0*log(tau) = 0
to the VALUE but log(tau) ~= -13.8 to the GRADIENT: the objective is
insensitive to LOGTOL while the descent direction is not.
`test_optimum_is_insensitive_to_logtol` measures how much that actually moves
the optimum (under 1e-6 over four orders of magnitude in tau on these
fixtures), and every other tolerance here is set above it.

Existing `qrf.mmi` / `qrf.mem` / `qrf.mmi.linear` values are untouched: the new
token has its own objective and shares only the polytope.
"""
import numpy as np
import pytest
from scipy.linalg import null_space
from scipy.optimize import linprog

import line_solver.api.mapqn.qrf_noblo_common as C
from line_solver.api.mapqn.qrf_noblo_common import (
    affine_constraint_matrices, bethe_gradient, bethe_objective,
    build_q_from_mu_v_rt, compute_num_vars, extract_mu_v_from_maps,
    feasible_start, mem_objective, mmi_objective, qrf_index_map,
    reduce_equalities, solve_qrf_nlp, sub_qrfvar)
from line_solver.api.mapqn.qrf_noblo_bethe import qrf_noblo_bethe

_MR = 1


def _expmap(rate):
    return [np.array([[-rate]]), np.array([[rate]])]


def _erlang2(rate):
    """Erlang-2 of the given mean rate, phase advance in D0's upper triangle."""
    D0 = np.array([[-2 * rate, 2 * rate], [0.0, -2 * rate]])
    D1 = np.array([[0.0, 0.0], [2 * rate, 0.0]])
    return [D0, D1]


def _cyclic(M):
    rt = np.zeros((M, M))
    for i in range(M):
        rt[i, (i + 1) % M] = 1.0
    return rt


# name -> (MAPs, N). Kept small on purpose: the reduced free dimension is 0, 66,
# 3 and 9 here, so the four cover both the pinned polytope and a loose one, at
# K == 1 and K == 2, in a few seconds.
FIXTURES = {
    'm2k1_N3': ([_expmap(1.0), _expmap(1.5)], 3),
    'm2k2_N2': ([_expmap(1.0), _erlang2(1.5)], 2),
    'm3k1_N3': ([_expmap(1.0), _expmap(1.5), _expmap(2.0)], 3),
    'm3k1_N4': ([_expmap(1.0), _expmap(1.5), _expmap(2.0)], 4),
}


def _setup(name):
    """(M, N, K, F, Aeq, beq, Aub, bub, n) for one fixture."""
    MAPs, N = FIXTURES[name]
    M = len(MAPs)
    K = np.array([MAPs[i][0].shape[0] for i in range(M)], dtype=int)
    F = np.full(M, N, dtype=int)
    BB = np.zeros((_MR, M))
    mu, v = extract_mu_v_from_maps(MAPs, M, K)
    q = build_q_from_mu_v_rt(M, K, mu, v, _cyclic(M))
    n = compute_num_vars(M, N, K, _MR)
    Aeq, beq = affine_constraint_matrices(
        lambda z: C.sub_qrfcon_noblo(z, q, M, _MR, BB, F, N, K)[1], n)
    Aub, bub = affine_constraint_matrices(
        lambda z: C.sub_qrfcon_noblo(z, q, M, _MR, BB, F, N, K)[0], n)
    Aeq, beq, _ = reduce_equalities(Aeq, beq)
    return M, N, K, F, Aeq, beq, Aub, bub, n


def _solve(name, x0=None):
    M, N, K, F, Aeq, beq, Aub, bub, n = _setup(name)
    if x0 is None:
        x0 = feasible_start(Aeq, beq, Aub, bub, n)
    idx = qrf_index_map(M, N, K, _MR)
    return solve_qrf_nlp(
        lambda x: bethe_objective(x, M, N, K, F, _MR),
        lambda x: bethe_gradient(x, M, N, K, F, _MR, idx),
        x0, Aeq, beq, Aub, bub, 'qrf_noblo_bethe')


def _random_vertex(Aeq, beq, Aub, bub, n, seed):
    """A feasible point of the polytope that is NOT the phase-1 point.

    An LP with a random cost lands on a vertex, i.e. as far from the phase-1
    interior point as the polytope allows -- which is what makes it a real test
    of start-point independence rather than a perturbation of it.
    """
    rng = np.random.default_rng(seed)
    r = linprog(rng.normal(size=n),
                A_ub=Aub if len(Aub) else None, b_ub=bub if len(Aub) else None,
                A_eq=Aeq, b_eq=beq, bounds=[(0.0, 1.0)] * n, method='highs')
    assert r.success, 'random-vertex LP failed: %s' % r.message
    return r.x


def _u1_range(name):
    """glpsol-equivalent [min, max] of U at station 1 over the same polytope."""
    M, N, K, F, Aeq, beq, Aub, bub, n = _setup(name)
    idx = qrf_index_map(M, N, K, _MR)
    obj = np.zeros(n)
    for k in range(K[0]):
        for n1 in range(1, N + 1):
            obj[idx[0, n1, k, 0, n1, k, 0]] += 1.0
    out = []
    for s in (1.0, -1.0):
        r = linprog(s * obj, A_ub=Aub if len(Aub) else None,
                    b_ub=bub if len(Aub) else None, A_eq=Aeq, b_eq=beq,
                    bounds=[(0.0, 1.0)] * n, method='highs')
        assert r.success, 'LP failed: %s' % r.message
        out.append(s * r.fun)
    return out


@pytest.mark.parametrize('name', sorted(FIXTURES))
def test_gradient_matches_finite_differences(name):
    """bethe_gradient is the gradient of bethe_objective.

    Differenced in the REDUCED space (x0 + Z t with Z spanning null(Aeq)),
    which is where the solver actually works, and at a step small enough that
    no coordinate crosses zero -- a central difference through the boundary
    would measure the clip, not the derivative.

    WHERE THE EQUALITIES PIN A SINGLE POINT (m2k1_N3) the reduced space offers
    no direction at all. That used to skip, which asserted nothing: the
    identity under test is a property of the FUNCTION, and bethe_objective is
    defined on all of R^n rather than only on the polytope. So difference
    along the ambient coordinates instead, restricted to the ones strictly
    positive at the point -- 18 of 66 there, the smallest 0.123 -- which is
    what keeps the step clear of the boundary the reduced route avoids by
    staying inside.
    """
    M, N, K, F, Aeq, beq, Aub, bub, n = _setup(name)
    idx = qrf_index_map(M, N, K, _MR)
    x0 = feasible_start(Aeq, beq, Aub, bub, n)
    Z = null_space(Aeq)
    if Z.shape[1] == 0:
        Z = np.eye(n)[:, np.flatnonzero(x0 > 1e-6)]
        assert Z.shape[1] > 0, '%s: no interior coordinate to difference' % name
        x = x0
    else:
        rng = np.random.default_rng(0)
        x = x0 + Z @ (rng.normal(size=Z.shape[1]) * 1e-6)
    gz = Z.T @ bethe_gradient(x, M, N, K, F, _MR, idx)
    h = 1e-9
    for d in range(min(Z.shape[1], 12)):
        e = np.zeros(Z.shape[1])
        e[d] = h
        num = (bethe_objective(x + Z @ e, M, N, K, F, _MR)
               - bethe_objective(x - Z @ e, M, N, K, F, _MR)) / (2 * h)
        assert abs(num - gz[d]) < 1e-4 * max(1.0, abs(num)), (
            '%s direction %d: analytic %.9f vs finite difference %.9f'
            % (name, d, gz[d], num))


@pytest.mark.parametrize('name', ['m2k2_N2', 'm3k1_N3'])
def test_optimum_is_independent_of_the_start_point(name):
    """The property convexity buys, asserted directly.

    f_{1/M} is convex on this polytope, so every local minimum is global and
    the reported point must not depend on where the descent started. This is
    exactly what `qrf.mmi` cannot promise: its objective is concave in some
    directions of these same feasible sets, which is why its second digit moves
    with the restart schedule.
    """
    M, N, K, F, Aeq, beq, Aub, bub, n = _setup(name)
    base = _solve(name)
    fbase = bethe_objective(base, M, N, K, F, _MR)
    UNbase, _ = C.extract_results(sub_qrfvar(base, M, N, K, _MR)[0], M, K, F, _MR)
    for seed in (1, 2):
        x0 = _random_vertex(Aeq, beq, Aub, bub, n, seed)
        xo = _solve(name, x0)
        UN, _ = C.extract_results(sub_qrfvar(xo, M, N, K, _MR)[0], M, K, F, _MR)
        assert abs(bethe_objective(xo, M, N, K, F, _MR) - fbase) < 1e-8, (
            '%s seed %d: objective %.12f vs %.12f from the phase-1 point'
            % (name, seed, bethe_objective(xo, M, N, K, F, _MR), fbase))
        np.testing.assert_allclose(UN, UNbase, atol=1e-5, rtol=0,
                                   err_msg='%s seed %d' % (name, seed))


@pytest.mark.parametrize('name', sorted(FIXTURES))
def test_optimum_is_feasible(name):
    """The returned point satisfies the constraints it was optimised over."""
    M, N, K, F, Aeq, beq, Aub, bub, n = _setup(name)
    x = _solve(name)
    assert np.max(np.abs(Aeq @ x - beq)) < 1e-7
    if len(Aub):
        assert np.max(Aub @ x - bub) < 1e-7
    assert np.min(x) >= -1e-12 and np.max(x) <= 1.0 + 1e-12


def test_optimum_is_insensitive_to_logtol():
    """The caveat that keeps digits out of these tests, measured.

    The restored n = 0 cells are structurally zero, so they add 0*log(tau) = 0
    to the objective VALUE but log(tau) to the GRADIENT. The value is therefore
    insensitive to LOGTOL while the descent direction is not, and the optimum
    has to be checked rather than assumed stable. It moves by well under 1e-6
    over four orders of magnitude here, which is why the tolerances above are
    1e-5 and not 1e-9.
    """
    name = 'm3k1_N3'
    M, N, K, F, _Aeq, _beq, _Aub, _bub, _n = _setup(name)
    saved = C.LOGTOL
    seen = []
    try:
        for tau in (1e-4, 1e-6, 1e-8):
            C.LOGTOL = tau
            x = _solve(name)
            UN, _ = C.extract_results(sub_qrfvar(x, M, N, K, _MR)[0], M, K, F, _MR)
            seen.append(np.asarray(UN, dtype=float))
    finally:
        C.LOGTOL = saved
    spread = max(float(np.max(np.abs(u - seen[0]))) for u in seen)
    assert spread < 1e-6, 'optimum moved by %.3e across LOGTOL' % spread


def test_the_idle_cells_are_inside_the_sum():
    """`qrf.bethe` decomposes exactly as lambda*MMI + MEM, n = 0 cells included.

    The three bodies do NOT all share a lower bound, and the reason is the AMPL
    source each transcribes. The MI line reads `sum {ni, nj in 0..F}`, so
    `mmi_objective` starting at n = 1 was a transcription defect (D1) and was
    corrected on 2026-09-02. The MEM line reads `sum {ni in 1..F[i]}`, so
    `mem_objective` starting at n = 1 MATCHES its own spec and was left alone.
    `bethe_objective` sums both terms from n = 0 by construction.

    Hence bethe == lambda*MMI + MEM + (the MEM n = 0 block), and that residual
    is a real number rather than a structural zero: p_ii at n = 0 is
    P(n_i = 0) = 1 - U_i, far from zero. If a future edit reintroduced the n = 1
    bound in the MI body, `coded` would fall short by the MI n = 0 block as well
    and the final assertion would fail.
    """
    name = 'm3k1_N3'
    M, N, K, F, Aeq, beq, Aub, bub, n = _setup(name)
    x = feasible_start(Aeq, beq, Aub, bub, n)
    p2, _ = sub_qrfvar(x, M, N, K, _MR)
    lam = 1.0 / M

    # The MEM n = 0 block: the ONLY piece of bethe that lambda*MMI + MEM lacks.
    idle = 0.0
    for i in range(M):
        for k in range(K[i]):
            pv = p2[i, 0, k, i, 0, k, 0]
            idle += pv * np.log(C.LOGTOL + pv)

    coded = (lam * mmi_objective(x, M, N, K, F, _MR)
             + mem_objective(x, M, N, K, F, _MR))
    assert abs(idle) > 1e-3, 'the n = 0 block is numerically empty on %s' % name
    assert abs(bethe_objective(x, M, N, K, F, _MR) - coded - idle) < 1e-10


@pytest.mark.parametrize('name', sorted(FIXTURES))
def test_within_the_lp_range_of_its_own_polytope(name):
    """Containment, the one invariant any point of this polytope must satisfy.

    Where the LP range is degenerate the polytope pins the answer and this is
    an equality test against the exact value; where it is not, containment is
    all that is claimed -- the objective selects a point, it does not bound.
    """
    MAPs, N = FIXTURES[name]
    M = len(MAPs)
    K = np.array([MAPs[i][0].shape[0] for i in range(M)], dtype=int)
    mu, v = extract_mu_v_from_maps(MAPs, M, K)
    UN, QN = qrf_noblo_bethe(M, _MR, K, N, mu, v, _cyclic(M))
    lo, hi = _u1_range(name)
    assert lo - 1e-6 <= UN[0] <= hi + 1e-6, (
        '%s: U1 = %.9f outside the LP range [%.9f, %.9f]' % (name, UN[0], lo, hi))
    assert abs(sum(QN) - N) < 1e-6, '%s: population %s' % (name, QN)
    assert np.all(np.asarray(UN) > 1e-6) and np.all(np.asarray(UN) <= 1.0 + 1e-9)


def test_token_is_advertised_and_runs_end_to_end():
    """SolverBA -> solver_ctmc_qrf_analyzer -> qrf_noblo_bethe."""
    from line_solver import (Network, Queue, ClosedClass, SchedStrategy, Exp)
    from line_solver.solvers.solver_ba.solver_ba import BA_METHODS
    from line_solver.api.solvers.ctmc.solver_ctmc_qrf_analyzer import (
        solver_ctmc_qrf_analyzer)

    assert 'qrf.bethe' in BA_METHODS

    model = Network('bethe_chain')
    queues = [Queue(model, 'Q%d' % (i + 1), SchedStrategy.FCFS) for i in range(3)]
    cls = ClosedClass(model, 'C1', 3, queues[0])
    for q, rate in zip(queues, (1.0, 1.5, 2.0)):
        q.setService(cls, Exp(rate))
    model.link(model.serialRouting(*queues))

    class _Opts(object):
        method = 'qrf.bethe'
        config = {}

    QN, UN, _RN, _TN, _CN, XN, _rt = solver_ctmc_qrf_analyzer(
        model.getStruct(), _Opts())
    UN = np.asarray(UN, dtype=float).flatten()
    QN = np.asarray(QN, dtype=float).flatten()
    assert np.all(UN > 1e-6) and np.all(UN <= 1.0 + 1e-9)
    assert abs(QN.sum() - 3.0) < 1e-6
    assert float(np.asarray(XN, dtype=float).flatten()[0]) > 0.0
