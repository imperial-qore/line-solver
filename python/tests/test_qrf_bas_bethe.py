"""Regression tests for `qrf.bas.bethe`, the Bethe arm of the BAS-blocking QRF.

WHAT THE METHOD IS
------------------
The objective of `qrf.bethe` evaluated over the BAS decision vector:

    lambda * sum_{i!=j} I(n_i;n_j) - sum_i H(n_i),      lambda = 1/M

over the polytope `qrf.bas` and `qrf.bas.mmi` already use, so the blocking
configurations and the per-station capacities enter through the ranges alone.
Both population sums run from n = 0, as they do in `qrf_noblo_bethe` and unlike
`qrf.bas.mem`, whose entropy starts at 1: the idle cell is the strongest
correlation in a closed chain, and an entropy taken over a different range than
the mutual information it is combined with is not a free entropy of anything.

WHAT IS ASSERTED, AND WHY CONVEXITY IS MEASURED HERE RATHER THAN ASSUMED
-----------------------------------------------------------------------
lambda = 1/M is the largest uniform edge weight at which the tree-reweighted
entropy is concave on the LOCAL MARGINAL polytope, which is what makes the
no-blocking arm's optimum a property of the model rather than of the start
point. The BAS polytope adds the blocking families on top of the marginal ones,
so it is a convex subset and the objective's convexity carries -- but the
per-configuration marginal consistency the argument rests on is not proved
under ZERO5, which pins a blocked station's idle cell to zero. So
`test_optimum_is_independent_of_the_start_point` MEASURES start-point
independence on these fixtures instead of taking it on trust.

The accuracy assertions are containment and ordering, not frozen digits: the
value is insensitive to LOGTOL while the descent direction is not, the same
reason test_qrf_noblo_bethe.py gives.
"""
import numpy as np
import pytest
from scipy.optimize import linprog

from line_solver import (ClosedClass, DropStrategy, Exp, Network, Queue,
                         SchedStrategy, SolverBA, SolverCTMC)
from line_solver.api.mapqn.qrf_bas_nlp import (
    _bas_polytope, _bethe_pair, _mem_terms, _metrics, _mmi_terms, _p2_columns,
    qrf_bas_bethe)
from line_solver.api.mapqn.parameters import QRBoundsBasParameters
from line_solver.api.mapqn.qrf_noblo_common import solve_qrf_nlp
from line_solver.api.sn.qrf_blocking import sn_to_qrf_blocking


def _cqn_bas_blocking():
    """Queue1 -BAS-> Queue2(cap 1), N = 2. The in-tree example."""
    m = Network('cqn_bas_blocking')
    q1 = Queue(m, 'Queue1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Queue2', SchedStrategy.FCFS)
    c = ClosedClass(m, 'Class1', 2, q1, 0)
    q1.set_service(c, Exp(1.0))
    q2.set_service(c, Exp(0.8))
    q2.setCap(1)
    q1.setDropRule(c, DropStrategy.BAS)
    m.link(Network.serialRouting(q1, q2))
    return m


def _two_feeders(N=3, cap3=1):
    """Q1 -> {Q2, Q3(cap 1)}, Q2 -> Q3, Q3 -> Q1: two queues can block behind Q3."""
    m = Network('twoFeeders')
    qs = [Queue(m, 'Q%d' % (i + 1), SchedStrategy.FCFS) for i in range(3)]
    c = ClosedClass(m, 'C', N, qs[0], 0)
    for q, r in zip(qs, (1.0, 0.9, 1.2)):
        q.set_service(c, Exp(r))
    qs[2].setCap(cap3)
    for q in qs[:2]:
        q.setDropRule(c, DropStrategy.BAS)
    P = m.initRoutingMatrix()
    P.set(c, c, qs[0], qs[1], 0.5)
    P.set(c, c, qs[0], qs[2], 0.5)
    P.set(c, c, qs[1], qs[2], 1.0)
    P.set(c, c, qs[2], qs[0], 1.0)
    m.link(P)
    return m


FIXTURES = {'cqn_bas_blocking': _cqn_bas_blocking, 'two_feeders': _two_feeders}


def _params(make):
    """The QRBoundsBasParameters the token derives for one fixture."""
    model = make()
    sn = model.getStruct()
    qp, msg = sn_to_qrf_blocking(sn)
    assert msg == '', 'fixture is not a QRF BAS model: %s' % msg
    M = int(sn.nstations)
    N = int(np.sum(np.asarray(sn.njobs).ravel()))
    K = np.ones(M, dtype=int)
    mu = [np.array([[float(sn.rates[i, 0])]]) for i in range(M)]
    v = [np.zeros((1, 1)) for _ in range(M)]
    rt = np.asarray(sn.rt, dtype=float)[:M, :M]
    return QRBoundsBasParameters(
        _M=M, _N=N, MR=qp['MR'], f=qp['f'], K=K, F=np.asarray(qp['F']),
        MM=np.asarray(qp['MM']), MM1=np.asarray(qp['MM1']),
        ZZ=np.asarray(qp['ZZ']), BB=np.asarray(qp['BB']),
        mu=mu, v=v, r=rt)


def _random_vertex(Aeq, beq, Aub, bub, n, seed):
    """A vertex of the BAS polytope that is not the phase-1 point."""
    rng = np.random.default_rng(seed)
    r = linprog(rng.standard_normal(n),
                A_ub=Aub if len(Aub) else None, b_ub=bub if len(Aub) else None,
                A_eq=Aeq if len(Aeq) else None, b_eq=beq if len(Aeq) else None,
                bounds=[(0.0, 1.0)] * n, method='highs')
    assert r.success, 'random-vertex LP failed: %s' % r.message
    return r.x


def _solve_from(params, x0=None):
    """(UN, QN, f) of the Bethe arm, optionally from a supplied start point."""
    Aeq, beq, Aub, bub, model = _bas_polytope(params)
    cols = _p2_columns(model, params)
    ij, ii, jj = _mmi_terms(cols, params, n_from=0)
    diag = _mem_terms(cols, params, n_from=0)
    objective, gradient = _bethe_pair(ij, ii, jj, diag, params.M)
    if x0 is None:
        n = Aeq.shape[1] if Aeq.size else Aub.shape[1]
        r = linprog(np.zeros(n), A_ub=Aub if len(Aub) else None,
                    b_ub=bub if len(Aub) else None,
                    A_eq=Aeq if len(Aeq) else None,
                    b_eq=beq if len(Aeq) else None,
                    bounds=[(0.0, 1.0)] * n, method='highs')
        assert r.success
        x0 = r.x
    x = solve_qrf_nlp(objective, gradient, x0, Aeq, beq, Aub, bub, 'qrf_bas_bethe')
    UN, QN = _metrics(x, cols, model, params)
    return UN, QN, float(objective(np.clip(x, 0.0, 1.0)))


@pytest.mark.parametrize('name', sorted(FIXTURES))
def test_token_runs_end_to_end(name):
    """SolverBA must serve the method name by name on a blocked model."""
    util = np.asarray(SolverBA(FIXTURES[name](), 'qrf.bas.bethe').getAvgUtil()).ravel()
    assert util.size == 3 if name == 'two_feeders' else util.size == 2
    assert np.all(np.isfinite(util))
    assert np.all(util >= -1e-9) and np.all(util <= 1.0 + 1e-9)


def test_the_token_is_advertised_on_a_blocked_model():
    """It is a blocking-aware bound, so the blocked-model list must keep it."""
    valid = SolverBA(_cqn_bas_blocking()).list_valid_methods()
    assert 'qrf.bas.bethe' in valid


def test_pinned_polytope_reproduces_the_exact_chain():
    """On cqn_bas_blocking the polytope pins U, so every arm must be exact.

    States (n1,n2,blocked) give pi = (0.262295, 0.327869, 0.409836); Queue1 is
    serving in the first two and BLOCKED in the third, so U1 = 0.590164, and
    Queue2 is busy in the last two, U2 = 0.737705. The LP token already returns
    these, and an objective chosen over a polytope with one feasible utilization
    cannot return anything else.
    """
    util = np.asarray(SolverBA(_cqn_bas_blocking(), 'qrf.bas.bethe').getAvgUtil()).ravel()
    assert util[0] == pytest.approx(0.590164, abs=1e-5)
    assert util[1] == pytest.approx(0.737705, abs=1e-5)


@pytest.mark.parametrize('name', sorted(FIXTURES))
def test_optimum_is_independent_of_the_start_point(name):
    """Convexity, measured. See the module docstring for why it is not assumed."""
    params = _params(FIXTURES[name])
    Aeq, beq, Aub, bub, model = _bas_polytope(params)
    n = Aeq.shape[1] if Aeq.size else Aub.shape[1]

    U_ref, Q_ref, f_ref = _solve_from(params)
    for seed in (7, 19, 31):
        U, Q, f = _solve_from(params, _random_vertex(Aeq, beq, Aub, bub, n, seed))
        assert f == pytest.approx(f_ref, abs=1e-6), 'seed %d moved the optimum' % seed
        assert np.allclose(U, U_ref, atol=1e-5), 'seed %d moved U' % seed


@pytest.mark.parametrize('name', sorted(FIXTURES))
def test_optimum_is_feasible(name):
    """The reported point must lie in the polytope it was optimised over."""
    params = _params(FIXTURES[name])
    Aeq, beq, Aub, bub, model = _bas_polytope(params)
    cols = _p2_columns(model, params)
    ij, ii, jj = _mmi_terms(cols, params, n_from=0)
    diag = _mem_terms(cols, params, n_from=0)
    objective, gradient = _bethe_pair(ij, ii, jj, diag, params.M)
    n = Aeq.shape[1] if Aeq.size else Aub.shape[1]
    r = linprog(np.zeros(n), A_ub=Aub if len(Aub) else None,
                b_ub=bub if len(Aub) else None,
                A_eq=Aeq if len(Aeq) else None, b_eq=beq if len(Aeq) else None,
                bounds=[(0.0, 1.0)] * n, method='highs')
    x = solve_qrf_nlp(objective, gradient, r.x, Aeq, beq, Aub, bub, 'qrf_bas_bethe')

    if Aeq.size:
        assert np.max(np.abs(Aeq @ x - beq)) < 1e-6
    if len(Aub):
        assert np.max(Aub @ x - bub) < 1e-6
    assert x.min() > -1e-9


def test_the_idle_cells_are_inside_the_sums():
    """Both blocks must span n = 0, which is what separates this from mem+mmi.

    Asserted on the term sets rather than on a value: n_from=0 has to reach BOTH
    the mutual-information triples and the entropy diagonal, and a regression
    that restored only one of them would still produce a plausible number.
    """
    params = _params(_cqn_bas_blocking)
    _, _, _, _, model = _bas_polytope(params)
    cols = _p2_columns(model, params)

    ij0, _, _ = _mmi_terms(cols, params, n_from=0)
    ij1, _, _ = _mmi_terms(cols, params, n_from=1)
    d0 = _mem_terms(cols, params, n_from=0)
    d1 = _mem_terms(cols, params, n_from=1)
    assert ij0.size > ij1.size
    assert d0.size > d1.size


def test_bethe_is_the_lambda_weighted_sum_of_the_other_two_bodies():
    """The objective must BE lam*MI + entropy, evaluated on a real point."""
    params = _params(_two_feeders)
    Aeq, beq, Aub, bub, model = _bas_polytope(params)
    cols = _p2_columns(model, params)
    ij, ii, jj = _mmi_terms(cols, params, n_from=0)
    diag = _mem_terms(cols, params, n_from=0)
    obj, _ = _bethe_pair(ij, ii, jj, diag, params.M)

    n = Aeq.shape[1] if Aeq.size else Aub.shape[1]
    x = _random_vertex(Aeq, beq, Aub, bub, n, 5)

    from line_solver.api.mapqn.qrf_bas_nlp import _mem_pair, _mmi_pair
    mi_obj, _ = _mmi_pair(ij, ii, jj)
    en_obj, _ = _mem_pair(diag)
    assert obj(x) == pytest.approx(mi_obj(x) / params.M + en_obj(x), rel=1e-12)


def test_it_is_not_looser_than_the_lp_on_the_feeder_model():
    """The point selection must beat the LP vertex it starts from.

    The LP token optimises one station's utilization and lands on a vertex; an
    information objective over the same polytope picks an interior point that is
    closer to the exact chain. This is the reason to have the arm at all, so a
    regression that quietly reverted to the vertex has to fail here.
    """
    model = _two_feeders()
    exact = np.asarray(SolverCTMC(_two_feeders()).getAvgUtil()).ravel()
    lp = np.asarray(SolverBA(_two_feeders(), 'qrf.bas').getAvgUtil()).ravel()
    bethe = np.asarray(SolverBA(_two_feeders(), 'qrf.bas.bethe').getAvgUtil()).ravel()
    assert np.max(np.abs(bethe - exact)) < np.max(np.abs(lp - exact))
