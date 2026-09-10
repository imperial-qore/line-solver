"""Regression test for `solve_qrf_nlp` refusing a start it can escape from.

WHAT WENT WRONG
---------------
`solve_qrf_nlp` guards against SLSQP silently returning the phase-1 LP vertex
as though it were the QRF optimum: if a feasible descent direction exists at
the start and the solver did not move, the reported metrics would be those of
the vertex. The guard was asked BEFORE the vertex-probe escape and raised
`RuntimeError` on the spot -- but "a descent direction exists and SLSQP did not
move" is not by itself a solver failure. The phase-1 LP returns a VERTEX, SLSQP
declines a degenerate one with status 4 ("Inequality constraints incompatible")
without taking a step, and the finite step the escape takes from that same
point descends perfectly well. The guard now records the condition and raises
only if the escape cannot improve either, which keeps what it is for without
discarding an instance that is merely awkward for the QP subproblem.

WHY THIS FIXTURE AND NOT A SMALLER ONE
--------------------------------------
The condition is rare, which is how it survived: a search over the 3- and
4-station cyclic models at N in {2,3,4}, all three no-blocking objectives and
60 random-cost phase-1 vertices each (720 combinations) triggers it nowhere,
and at 5 stations only `mem` at N = 3 from the seed-22 vertex does, out of 120.
So the fixture is that instance, and it costs the ~15 s the constraint assembly
takes at 405 variables.

WHAT IS ASSERTED
----------------
That the refused vertex now reaches the SAME optimum as the ordinary phase-1
start, not merely that no exception escapes. Five other vertices of this
polytope always converged to it; the sixth is the one that used to raise.
"""
import numpy as np
import pytest
from scipy.optimize import linprog

import line_solver.api.mapqn.qrf_noblo_common as C
from line_solver.api.mapqn.qrf_noblo_common import (
    affine_constraint_matrices, build_q_from_mu_v_rt, compute_num_vars,
    extract_results, feasible_start, mem_gradient, mem_objective,
    qrf_index_map, reduce_equalities, solve_qrf_nlp, sub_qrfvar)

_MR = 1
_RATES = [1.0, 2.0, 4.0, 8.0, 16.0]
_N = 3
_STUCK_SEED = 22          # the one vertex of this polytope SLSQP will not leave


def _setup():
    """(M, N, K, F, Aeq, beq, Aub, bub, n) for the 5-station cyclic fixture."""
    M = len(_RATES)
    K = np.ones(M, dtype=int)
    F = np.full(M, _N, dtype=int)
    BB = np.zeros((_MR, M))
    mu = np.zeros((M, 1, 1))
    v = np.zeros((M, 1, 1))
    rt = np.zeros((M, M))
    for i in range(M):
        mu[i, 0, 0] = _RATES[i]
        rt[i, (i + 1) % M] = 1.0
    q = build_q_from_mu_v_rt(M, K, mu, v, rt)
    n = compute_num_vars(M, _N, K, _MR)
    Aeq, beq = affine_constraint_matrices(
        lambda z: C.sub_qrfcon_noblo(z, q, M, _MR, BB, F, _N, K)[1], n)
    Aub, bub = affine_constraint_matrices(
        lambda z: C.sub_qrfcon_noblo(z, q, M, _MR, BB, F, _N, K)[0], n)
    Aeq, beq, _ = reduce_equalities(Aeq, beq)
    return M, _N, K, F, Aeq, beq, Aub, bub, n


def _random_vertex(Aeq, beq, Aub, bub, n, seed):
    """A vertex of the polytope other than the phase-1 one."""
    rng = np.random.default_rng(seed)
    r = linprog(rng.standard_normal(n),
                A_ub=Aub if len(Aub) else None, b_ub=bub if len(Aub) else None,
                A_eq=Aeq, b_eq=beq, bounds=[(0.0, 1.0)] * n, method='highs')
    assert r.success, 'random-vertex LP failed: %s' % r.message
    return r.x


def _solve_mem(setup, x0):
    M, N, K, F, Aeq, beq, Aub, bub, n = setup
    idx = qrf_index_map(M, N, K, _MR)
    xopt = solve_qrf_nlp(lambda x: mem_objective(x, M, N, K, F, _MR),
                         lambda x: mem_gradient(x, M, N, K, F, _MR, idx),
                         x0, Aeq, beq, Aub, bub, 'qrf_noblo_mem')
    p2opt, _ = sub_qrfvar(xopt, M, N, K, _MR)
    UN, QN = extract_results(p2opt, M, K, F, _MR)
    return UN, QN, float(mem_objective(np.clip(xopt, 0.0, 1.0), M, N, K, F, _MR))


@pytest.fixture(scope='module')
def setup():
    return _setup()


def test_a_vertex_slsqp_will_not_leave_still_reaches_the_optimum(setup):
    """The seed-22 vertex used to raise; it must now agree with the default."""
    M, N, K, F, Aeq, beq, Aub, bub, n = setup
    U_ref, Q_ref, f_ref = _solve_mem(setup, feasible_start(Aeq, beq, Aub, bub, n))
    x0 = _random_vertex(Aeq, beq, Aub, bub, n, _STUCK_SEED)
    U, Q, f = _solve_mem(setup, x0)

    assert f == pytest.approx(f_ref, abs=1e-6)
    assert np.allclose(U, U_ref, atol=1e-6)
    assert np.allclose(Q, Q_ref, atol=1e-6)
    # It really is a different start, not a perturbation of the phase-1 point.
    assert np.linalg.norm(x0 - feasible_start(Aeq, beq, Aub, bub, n)) > 0.1


def test_the_escape_moves_off_the_start(setup):
    """The answer must not be the start point dressed up as an optimum.

    This is what the guard exists to prevent, and it is the half of the guard
    the fix keeps: the returned point has to differ from the vertex it began
    at, and its objective has to be strictly lower.
    """
    M, N, K, F, Aeq, beq, Aub, bub, n = setup
    x0 = _random_vertex(Aeq, beq, Aub, bub, n, _STUCK_SEED)
    idx = qrf_index_map(M, N, K, _MR)
    f0 = float(mem_objective(np.clip(x0, 0.0, 1.0), M, N, K, F, _MR))
    xopt = solve_qrf_nlp(lambda x: mem_objective(x, M, N, K, F, _MR),
                         lambda x: mem_gradient(x, M, N, K, F, _MR, idx),
                         x0, Aeq, beq, Aub, bub, 'qrf_noblo_mem')
    f1 = float(mem_objective(np.clip(xopt, 0.0, 1.0), M, N, K, F, _MR))

    assert np.linalg.norm(xopt - x0) > 1e-8
    assert f1 < f0 - 1e-8


def test_the_optimum_is_feasible_from_the_refused_vertex(setup):
    """Escaping must not step outside the polytope it was escaping within."""
    M, N, K, F, Aeq, beq, Aub, bub, n = setup
    x0 = _random_vertex(Aeq, beq, Aub, bub, n, _STUCK_SEED)
    idx = qrf_index_map(M, N, K, _MR)
    xopt = solve_qrf_nlp(lambda x: mem_objective(x, M, N, K, F, _MR),
                         lambda x: mem_gradient(x, M, N, K, F, _MR, idx),
                         x0, Aeq, beq, Aub, bub, 'qrf_noblo_mem')

    assert np.max(np.abs(np.asarray(Aeq) @ xopt - np.asarray(beq))) < 1e-6
    if len(Aub):
        assert np.max(np.asarray(Aub) @ xopt - np.asarray(bub)) < 1e-6
    assert xopt.min() > -1e-9 and xopt.max() < 1.0 + 1e-9
