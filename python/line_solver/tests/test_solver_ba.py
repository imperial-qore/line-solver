"""Tests for SolverBA (bound-analysis solver) and the SolverMVA 'marie' method."""
import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, ClosedClass, Exp, Coxian, Erlang,
                         SchedStrategy, SolverBA, SolverMVA, SolverCTMC)


def _cqn(N=6, m1=1.0, m2=1.5, Z=2.0, sched=SchedStrategy.PS,
         scv1=1.0, scv2=1.0):
    m = Network('cqn')
    d = Delay(m, 'Z'); q1 = Queue(m, 'Q1', sched); q2 = Queue(m, 'Q2', sched)
    c = ClosedClass(m, 'C', N, d, 0)
    d.setService(c, Exp(1.0 / Z))
    q1.setService(c, Exp(1.0 / m1) if scv1 == 1.0 else Coxian.fitMeanAndSCV(m1, scv1))
    q2.setService(c, Exp(1.0 / m2) if scv2 == 1.0 else Coxian.fitMeanAndSCV(m2, scv2))
    m.link(Network.serialRouting(d, q1, q2))
    return m, c


def _chainX(solver):
    return float(np.asarray(solver.get_bounds()['Tupper'])[1, 0]), \
        float(np.asarray(solver.get_bounds()['Tlower'])[1, 0])


@pytest.mark.parametrize('fam', ['aba', 'bjb', 'pb', 'gb', 'pbh', 'pbk', 'bjbk', 'cbh'])
def test_single_class_bounds_bracket_exact(fam):
    m, _ = _cqn()
    exX = float(SolverMVA(m, 'exact').getAvgTable().Tput.iloc[0])
    b = SolverBA(m, method=fam + '.upper', level=6).get_bounds()
    xl = float(np.asarray(b['Tlower'])[1, 0])
    xu = float(np.asarray(b['Tupper'])[1, 0])
    assert xl - 1e-9 <= exX <= xu + 1e-9


def test_hierarchy_converges_to_exact():
    m, _ = _cqn()
    exX = float(SolverMVA(m, 'exact').getAvgTable().Tput.iloc[0])
    b = SolverBA(m, method='pbh.upper', level=6).get_bounds()
    assert abs(float(np.asarray(b['Tlower'])[1, 0]) - exX) < 1e-6
    assert abs(float(np.asarray(b['Tupper'])[1, 0]) - exX) < 1e-6


def test_multiclass_cub_mbjb_bracket():
    m = Network('mc')
    d = Delay(m, 'Z'); q = Queue(m, 'Q', SchedStrategy.PS)
    c1 = ClosedClass(m, 'A', 2, d, 0); c2 = ClosedClass(m, 'B', 2, d, 0)
    d.setService(c1, Exp(1)); d.setService(c2, Exp(1))
    q.setService(c1, Exp(2)); q.setService(c2, Exp(1.5))
    P = m.initRoutingMatrix()
    P[(c1, c1)] = Network.serialRouting(d, q)
    P[(c2, c2)] = Network.serialRouting(d, q)
    m.link(P)
    Tex = np.asarray(SolverMVA(m, 'exact').getAvg()[3])[1, :]
    ub = SolverBA(m, method='cub.upper'); ub.runAnalyzer()
    lb = SolverBA(m, method='mbjb.lower'); lb.runAnalyzer()
    Tub = np.asarray(ub._result['TN'])[1, :]
    Tlb = np.asarray(lb._result['TN'])[1, :]
    assert np.all(Tlb <= Tex + 1e-9) and np.all(Tex <= Tub + 1e-9)


def test_cub_get_bounds_one_sided_nan():
    m = Network('mc')
    d = Delay(m, 'Z'); q = Queue(m, 'Q', SchedStrategy.PS)
    c1 = ClosedClass(m, 'A', 2, d, 0); c2 = ClosedClass(m, 'B', 2, d, 0)
    d.setService(c1, Exp(1)); d.setService(c2, Exp(1))
    q.setService(c1, Exp(2)); q.setService(c2, Exp(1.5))
    P = m.initRoutingMatrix()
    P[(c1, c1)] = Network.serialRouting(d, q)
    P[(c2, c2)] = Network.serialRouting(d, q)
    m.link(P)
    b = SolverBA(m, method='cub.upper').get_bounds()
    assert np.all(np.isnan(np.asarray(b['Tlower'])))
    assert np.all(np.isfinite(np.asarray(b['Tupper'])))


def test_marie_exponential_matches_exact_mva():
    m, _ = _cqn(sched=SchedStrategy.FCFS)
    Qex = np.asarray(SolverMVA(m, 'exact').getAvg()[0])
    Qma = np.asarray(SolverMVA(m, method='marie').getAvg()[0])
    assert np.max(np.abs(Qex - Qma)) < 1e-9


def test_marie_coxian_close_to_ctmc():
    m, _ = _cqn(N=4, sched=SchedStrategy.FCFS, scv1=0.5, scv2=2.0)
    Qct = np.asarray(SolverCTMC(m).getAvg()[0])
    Qma = np.asarray(SolverMVA(m, method='marie').getAvg()[0])
    assert np.max(np.abs(Qct - Qma)) / np.max(Qct) < 0.05


@pytest.mark.parametrize('meth', ['bjb.upper', 'pbh.upper', 'cub.upper', 'sib.lower'])
def test_mva_redirects_bound_methods_to_ba(meth):
    m, _ = _cqn()
    with pytest.raises(ValueError, match='SolverBA'):
        SolverMVA(m, method=meth).runAnalyzer()


from line_solver import Erlang, HyperExp   # noqa: E402


def _cqn2(N=2, Z=1.0):
    """Two-class closed FCFS network with class-independent PH service per
    station (Erlang at Q1, HyperExp at Q2) so SolverCTMC is an exact reference
    for the non-product-form Marie approximation."""
    m = Network('cqn2')
    d = Delay(m, 'Z'); q1 = Queue(m, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Q2', SchedStrategy.FCFS)
    c1 = ClosedClass(m, 'C1', N, d, 0); c2 = ClosedClass(m, 'C2', N, d, 0)
    d.setService(c1, Exp(1.0 / Z)); d.setService(c2, Exp(1.0 / Z))
    q1.setService(c1, Erlang.fitMeanAndSCV(1.0, 0.5))
    q1.setService(c2, Erlang.fitMeanAndSCV(1.0, 0.5))
    q2.setService(c1, HyperExp.fitMeanAndSCV(1.0, 3.0))
    q2.setService(c2, HyperExp.fitMeanAndSCV(1.0, 3.0))
    P = m.initRoutingMatrix()
    P[c1] = Network.serialRouting(d, q1, q2)
    P[c2] = Network.serialRouting(d, q1, q2)
    m.link(P)
    return m


def test_marie_multiclass_exponential_matches_exact_mva():
    # Class-independent exponential FCFS multiclass -> genuine BCMP -> exact MVA.
    m = Network('cqn2e')
    d = Delay(m, 'Z'); q1 = Queue(m, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Q2', SchedStrategy.FCFS)
    c1 = ClosedClass(m, 'C1', 2, d, 0); c2 = ClosedClass(m, 'C2', 2, d, 0)
    for c in (c1, c2):
        d.setService(c, Exp(1.0)); q1.setService(c, Exp(1.0))
        q2.setService(c, Exp(1.0 / 1.5))
    P = m.initRoutingMatrix()
    P[c1] = Network.serialRouting(d, q1, q2)
    P[c2] = Network.serialRouting(d, q1, q2)
    m.link(P)
    Qex = np.asarray(SolverMVA(m, 'exact').getAvg()[0])
    Qma = np.asarray(SolverMVA(m, method='marie').getAvg()[0])
    assert np.max(np.abs(Qex - Qma)) < 1e-9


def test_marie_multiclass_close_to_ctmc():
    m = _cqn2()
    Qct = np.asarray(SolverCTMC(m).getAvg()[0])
    Qma = np.asarray(SolverMVA(m, method='marie').getAvg()[0])
    # Non-product-form multiclass FCFS: AMVA-class heuristic, within ~10%.
    assert np.max(np.abs(Qct - Qma)) / np.max(Qct) < 0.10


@pytest.mark.parametrize('meth', ['qr', 'qrf.mmi', 'qrf.mmi.linear'])
def test_qrf_rejects_delay_stations(meth):
    # The NLP-backed QRF methods are supported natively since 2026-07-20, but
    # qrf_noblo_* models every station as a single server, so an infinite-server
    # station would silently be a different model. _cqn carries a Delay, hence
    # the rejection. End-to-end numerical coverage lives in
    # test_solver_ba_qrf_path.py.
    m, _ = _cqn()
    with pytest.raises(ValueError, match='infinite-server'):
        SolverBA(m, method=meth).runAnalyzer()


@pytest.mark.parametrize("rates,N", [([1.0, 1.0], 3), ([1.0, 2.0, 0.5], 3), ([1.0, 2.0], 5)])
def test_lr_lp_bounds_bracket_exact(rates, N):
    """lr is a pure-LP linear-reduction bound: lr.lower <= exact <= lr.upper."""
    m = Network('lr')
    qs = [Queue(m, 'Q%d' % (i + 1), SchedStrategy.PS) for i in range(len(rates))]
    for q in qs:
        q.setNumberOfServers(1)
    c = ClosedClass(m, 'C', N, qs[0], 0)
    for q, r in zip(qs, rates):
        q.setService(c, Exp(r))
    m.link(Network.serialRouting(*qs))

    def util(t):
        return np.asarray(t.Util if hasattr(t, 'Util') else t['Util'], dtype=float)

    exact = util(SolverMVA(m, 'exact').getAvgTable())
    lo = util(SolverBA(m, 'lr.lower').getAvgTable())
    hi = util(SolverBA(m, 'lr.upper').getAvgTable())
    for i in range(len(rates)):
        assert lo[i] - 1e-6 <= exact[i] <= hi[i] + 1e-6


def test_lr_registered_and_not_qrf_alias():
    """lr/lr.lower/lr.upper are advertised and no longer alias qrf.mmi.linear."""
    from line_solver.solvers.solver_ba.solver_ba import BA_METHODS
    for tok in ('lr', 'lr.lower', 'lr.upper'):
        assert tok in BA_METHODS


# ---- getBoundsTable ----

def test_bounds_table_agrees_with_get_bounds():
    """The table carries exactly the numbers get_bounds returns."""
    m, _ = _cqn(6)
    b = SolverBA(m, 'gb.upper').get_bounds()
    t = SolverBA(m, 'gb.upper').get_bounds_table()
    assert list(t.columns) == ['Station', 'JobClass', 'Qlower', 'Qupper',
                               'Tlower', 'Tupper']
    assert len(t) == 3
    for i in range(len(t)):
        assert abs(t['Qlower'][i] - np.asarray(b['Qlower'])[i, 0]) < 1e-12
        assert abs(t['Qupper'][i] - np.asarray(b['Qupper'])[i, 0]) < 1e-12
        assert abs(t['Tlower'][i] - np.asarray(b['Tlower'])[i, 0]) < 1e-12
        assert abs(t['Tupper'][i] - np.asarray(b['Tupper'])[i, 0]) < 1e-12


def test_bounds_table_brackets_are_ordered():
    m, _ = _cqn(6)
    t = SolverBA(m, 'gb.upper').get_bounds_table()
    assert np.all(np.asarray(t['Qlower']) <= np.asarray(t['Qupper']) + 1e-9)
    assert np.all(np.asarray(t['Tlower']) <= np.asarray(t['Tupper']) + 1e-9)


def test_bounds_table_one_sided_keeps_nan():
    """cub is upper-only: the lower side stays NaN, rows are not dropped."""
    m, _ = _cqn(6)
    t = SolverBA(m, 'cub.upper').get_bounds_table()
    assert len(t) > 0
    assert np.all(np.isnan(np.asarray(t['Qlower'], dtype=float)))
    assert np.all(np.isnan(np.asarray(t['Tlower'], dtype=float)))
    assert np.all(np.isfinite(np.asarray(t['Qupper'], dtype=float)))


def test_get_bounds_propagates_level():
    """Regression: get_bounds must re-run each side under the caller's options.

    Without propagation, options.level reverts to its default of 2 and a
    hierarchical family never tightens however high level is set. At level = N
    the bracket must collapse onto the exact solution.
    """
    N = 5
    m, _ = _cqn(N)
    exX = float(np.asarray(SolverMVA(m, 'exact').getAvgTable().Tput,
                           dtype=float)[1])

    b2 = SolverBA(m, 'pbh.upper', level=2).get_bounds()
    w2 = float(np.asarray(b2['Tupper'])[1, 0]) - float(np.asarray(b2['Tlower'])[1, 0])

    bN = SolverBA(m, 'pbh.upper', level=N).get_bounds()
    wN = float(np.asarray(bN['Tupper'])[1, 0]) - float(np.asarray(bN['Tlower'])[1, 0])

    assert w2 > 1e-6, 'level 2 bracket should be strictly loose'
    assert wN < w2 - 1e-9, 'level N bracket must be tighter than level 2'
    assert abs(float(np.asarray(bN['Tlower'])[1, 0]) - exX) < 1e-6
    assert abs(float(np.asarray(bN['Tupper'])[1, 0]) - exX) < 1e-6


def test_bounds_table_propagates_level():
    """get_bounds_table inherits the level fix by routing through get_bounds."""
    N = 5
    m, _ = _cqn(N)
    exX = float(np.asarray(SolverMVA(m, 'exact').getAvgTable().Tput,
                           dtype=float)[1])
    t = SolverBA(m, 'pbh.upper', level=N).get_bounds_table()
    assert abs(t['Tlower'][1] - exX) < 1e-6
    assert abs(t['Tupper'][1] - exX) < 1e-6


# ---------------------------------------------------------------------------
# LP-backed QRF reduction bounds (qrf.bas, qrf.rsrd)
# ---------------------------------------------------------------------------

# The paper's BAS blocking structure (M=5, MR=5) from tomacs_qrf/example_bas.mod,
# shrunk in N/F/K. glpsol on the matching .mod gives U1min=0.3, U1max=0.4.
# The full-size instance is exercised by test_qr_bounds_bas_full_instance below;
# this surrogate has K=1 and therefore cannot see phase-indexed defects.
_BAS_BB = np.array([[0, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 1, 0, 0, 1],
                    [0, 0, 0, 0, 1], [0, 1, 0, 0, 1]], dtype=float)
_BAS_R = np.array([[0.0, .5, 0.0, 0.0, .5], [.5, 0.0, .5, 0.0, 0.0],
                   [0.0, .5, 0.0, .5, 0.0], [0.0, 0.0, .5, 0.0, .5],
                   [.5, 0.0, 0.0, .5, 0.0]])


def _bas_params():
    from line_solver.api.mapqn.parameters import QRBoundsBasParameters
    M, MR = 5, 5
    K = [1] * M
    MM = np.zeros((MR, M)); MM[:, 0] = [0, 2, 2, 5, 5]; MM[:, 1] = [0, 0, 5, 0, 2]
    MM1 = np.zeros((MR, M)); MM1[1] = [-1, -1, -1, -1, 3]; MM1[3] = [-1, 5, -1, -1, -1]
    return QRBoundsBasParameters(
        _M=M, _N=2, MR=MR, f=1, K=np.array(K), F=np.array([1, 2, 2, 2, 2]),
        MM=MM, MM1=MM1, ZZ=np.array([0, 1, 2, 1, 2], dtype=int), BB=_BAS_BB,
        mu=[np.full((1, 1), 1.0) for _ in range(M)],
        v=[np.zeros((1, 1)) for _ in range(M)], r=_BAS_R)


@pytest.mark.parametrize('sense,expected', [('min', 0.3), ('max', 0.4)])
def test_qr_bounds_bas_matches_glpk(sense, expected):
    """qr_bounds_bas reproduces the AMPL/glpsol optimum on a BAS instance.

    Regression for the five constraints THM30/THM3/THM3f/THM3I/THM3L, absent
    from this port, and for the UEFF aggregation defect: AMPL indexes UEFF over
    j (one row per j), so summing over j asserted e = M*e and forced e == 0.
    """
    from line_solver.api.mapqn.qr_bounds_bas import mapqn_qr_bounds_bas
    got = mapqn_qr_bounds_bas(_bas_params(), 1, sense).objective_value
    assert abs(got - expected) < 1e-5, '%s: got %.10f want %.10f' % (sense, got, expected)


# The paper's reference instances, at full size. Ground truth is glpsol 5.0 on
# tomacs_qrf/example_bas.mod and example_rsrd.mod (10 printed digits), so the
# tolerance below is set by the reference printout, not by solver accuracy.
_QRF_MAP_MU = np.array([[1.016186e+00, 2.585708e-05],
                        [1.569888e-03, 1.413298e-02]])
_QRF_R5 = np.array([[0.0, .5, 0.0, 0.0, .5], [.5, 0.0, .5, 0.0, 0.0],
                    [0.0, .5, 0.0, .5, 0.0], [0.0, 0.0, .5, 0.0, .5],
                    [.5, 0.0, 0.0, .5, 0.0]])


def _bas_params_full():
    """example_bas.mod: M=5, N=10, f=1, F=[5,10,10,10,10], K=2, MR=5."""
    from line_solver.api.mapqn.parameters import QRBoundsBasParameters
    M, MR = 5, 5
    MM = np.array([[0, 0], [2, 0], [2, 5], [5, 0], [5, 2]], dtype=float)
    MM1 = np.array([[0, 0, 0, 0, 0], [-1, -1, -1, -1, 3], [0, 0, 0, 0, 0],
                    [-1, 5, -1, -1, -1], [0, 0, 0, 0, 0]], dtype=float)
    return QRBoundsBasParameters(
        _M=M, _N=10, MR=MR, f=1, K=np.array([2] * M),
        F=np.array([5, 10, 10, 10, 10]), MM=MM, MM1=MM1,
        ZZ=np.array([0, 1, 2, 1, 2], dtype=int), BB=_BAS_BB,
        mu=[_QRF_MAP_MU.copy() for _ in range(M)],
        v=[np.zeros((2, 2)) for _ in range(M)], r=_QRF_R5)


def _rsrd_params_full():
    """example_rsrd.mod: M=5, N=20, F=[5]*5, K=2, alpha=1."""
    from line_solver.api.mapqn.parameters import QRBoundsRsrdParameters
    M, N = 5, 20
    return QRBoundsRsrdParameters(
        _M=M, _N=N, F=np.array([5] * M), K=np.array([2] * M),
        mu=[_QRF_MAP_MU.copy() for _ in range(M)],
        v=[np.zeros((2, 2)) for _ in range(M)],
        alpha=np.ones((M, N)), r=_QRF_R5)


@pytest.mark.parametrize('sense,expected', [('min', 0.4776625414),
                                            ('max', 0.8077316296)])
def test_qr_bounds_bas_full_instance(sense, expected):
    """Full-size BAS instance reproduces glpsol on example_bas.mod.

    Two regressions are pinned here, neither of which the K=1 surrogate above
    can see:
      * MapqnLpModel must store constraints sparsely. This LP is 60510 columns
        by 60429 rows with 467250 nonzeros; one dense row per constraint is
        29 GB and dies with MemoryError.
      * THM30's blocked-handover term carries a coefficient
        sum_{p,w} q[f,w,y,p] that depends on the phase y, so it must be formed
        inside the y loop. Summing over y as well over-tightened the LP and
        raised U1min to 0.47946 against the true 0.47766.
    """
    from line_solver.api.mapqn.qr_bounds_bas import mapqn_qr_bounds_bas
    got = mapqn_qr_bounds_bas(_bas_params_full(), 1, sense).objective_value
    assert abs(got - expected) < 1e-8, '%s: got %.10f want %.10f' % (
        sense, got, expected)


@pytest.mark.parametrize('sense,expected', [('min', 0.8705798470),
                                            ('max', 1.0000000000)])
def test_qr_bounds_rsrd_full_instance(sense, expected):
    """Full-size RS-RD instance reproduces glpsol on example_rsrd.mod."""
    from line_solver.api.mapqn.qr_bounds_rsrd import mapqn_qr_bounds_rsrd
    got = mapqn_qr_bounds_rsrd(_rsrd_params_full(), 1, sense).objective_value
    assert abs(got - expected) < 1e-8, '%s: got %.10f want %.10f' % (
        sense, got, expected)


def test_lpmodel_is_sparse_and_sums_duplicate_terms():
    """MapqnLpModel assembles sparse matrices and accumulates repeated terms.

    add_term on the same variable must accumulate (the dense builder did
    coefficients[idx] += v); with triplet storage that becomes duplicate (i,j)
    entries, which COO -> CSR sums.
    """
    import scipy.sparse as sp
    from line_solver.api.mapqn.lpmodel import MapqnLpModel

    model = MapqnLpModel()
    model.add_variable('x', lb=0.0, ub=None)
    model.add_variable('y', lb=0.0, ub=None)
    b = model.constraint_builder()
    b.add_term('x', 1.0).add_term('x', 2.0).add_term('y', -1.0)
    model.add_constraint(b.eq(3.0))
    model.add_constraint(model.constraint_builder().add_term('y', 1.0).leq(5.0))

    A_ub, b_ub, A_eq, b_eq = model._assemble()
    assert sp.issparse(A_eq) and sp.issparse(A_ub)
    assert A_eq.shape == (1, 2) and A_ub.shape == (1, 2)
    assert np.allclose(A_eq.toarray(), [[3.0, -1.0]])
    assert np.allclose(b_eq, [3.0]) and np.allclose(b_ub, [5.0])
    # 3x - y = 3, y <= 5, minimize x  ->  x = 1 at y = 0
    sol = model.solve('x', minimize=True)
    assert abs(sol.objective_value - 1.0) < 1e-12


def test_lpmodel_negates_geq_rows():
    """'>=' rows are negated into '<=' form on both sides during assembly."""
    from line_solver.api.mapqn.lpmodel import MapqnLpModel

    model = MapqnLpModel()
    model.add_variable('x')
    model.add_constraint(model.constraint_builder().add_term('x', 2.0).geq(4.0))
    A_ub, b_ub, A_eq, b_eq = model._assemble()
    assert A_eq is None and b_eq is None
    assert np.allclose(A_ub.toarray(), [[-2.0]])
    assert np.allclose(b_ub, [-4.0])
    assert abs(model.solve('x', minimize=True).objective_value - 2.0) < 1e-12


def test_qrf_lp_methods_are_advertised():
    m, _ = _cqn(3)
    methods = SolverBA(m).listValidMethods()
    assert 'qrf.bas' in methods
    assert 'qrf.rsrd' in methods


# glpsol max of U[station] on the paper's AMPL no-blocking model
# (noblo_skel.mod), for the closed cyclic chains used below. The LP tokens call
# the backend with sense='max', so this is the value they must reproduce.
# SolverMVA('exact') is NOT a valid oracle for the phase-type chain: it is
# product-form and returns the all-exponential answer regardless of the Erlang.
_QRF_LP_GLPSOL_MAX = {
    'm2sym_N2': [0.666666667, 0.666666667],
    'm3asym_N3': [0.824385805, 0.549590537, 0.412192903],
    'm3pht_N3': [0.870629371, 0.580419580, 0.435314685],
}

_QRF_LP_CHAINS = {
    'm2sym_N2': (2, [(1.0, 1), (1.0, 1)]),
    'm3asym_N3': (3, [(1.0, 1), (1.5, 1), (2.0, 1)]),
    'm3pht_N3': (3, [(1.0, 1), (1.5, 2), (2.0, 1)]),
}


def _qrf_chain(name):
    N, spec = _QRF_LP_CHAINS[name]
    m = Network(name)
    qs = [Queue(m, 'Q%d' % (i + 1), SchedStrategy.FCFS) for i in range(len(spec))]
    c = ClosedClass(m, 'C', N, qs[0], 0)
    for q, (rate, phases) in zip(qs, spec):
        q.setService(c, Exp(rate) if phases == 1 else Erlang(rate * phases, phases))
    P = m.initRoutingMatrix()
    P[c, c] = Network.serialRouting(*qs)
    m.link(P)
    return m


@pytest.mark.parametrize('method', ['qrf.bas', 'qrf.rsrd'])
@pytest.mark.parametrize('name', sorted(_QRF_LP_CHAINS))
def test_qrf_lp_methods_run_end_to_end(method, name):
    """The LP-backed QRF tokens dispatch and return finite, valid bounds.

    Widened from a single symmetric 2-station chain to M=3, asymmetric rates and
    a phase-type station. Unlike the NLP tokens these optimize a face of the
    polytope, so the upper-bound direction IS guaranteed and is asserted.
    """
    U = np.asarray(SolverBA(_qrf_chain(name), method).getAvgTable().Util,
                   dtype=float).flatten()
    assert np.all(np.isfinite(U))
    assert np.all(U >= -1e-9) and np.all(U <= 1.0 + 1e-9)

    # Both tokens maximize over a relaxation that contains the true joint
    # distribution, so both are upper bounds on utilization. SolverCTMC is the
    # oracle, not SolverMVA('exact'): the latter is product-form and returns the
    # all-exponential answer regardless of the Erlang, so it cannot adjudicate
    # the phase-type chain.
    exU = np.asarray(SolverCTMC(_qrf_chain(name)).getAvgTable().Util,
                     dtype=float).flatten()
    assert np.all(U >= exU - 1e-6)

    # qrf.bas optimizes the same LP glpsol solves, so it must attain that
    # optimum. qrf.rsrd is a DIFFERENT relaxation (RS-RD, with the alpha
    # load-dependence) rather than a superset of the no-blocking polytope, so it
    # is not bounded below by the no-blocking maximum: measured 0.85968 against
    # 0.87063 on m3pht_N3. Only qrf.bas is checked against glpsol.
    if method == 'qrf.bas':
        np.testing.assert_allclose(U, np.array(_QRF_LP_GLPSOL_MAX[name]),
                                   rtol=1e-5)


@pytest.mark.parametrize('method', ['qr', 'qrf.mmi', 'qrf.mem', 'qrf.mmi.ld',
                                    'qrf.mmi.linear'])
def test_qrf_nlp_methods_are_advertised(method):
    """The no-blocking NLP QRF tokens are advertised natively since 2026-07-20.

    They were withheld while the SolverBA -> solver_ctmc_qrf_analyzer ->
    qrf_noblo_* path was uncovered; exercising it fixed a transposed v in the
    adapter and a throughput inversion valid only for a delay reference station.
    Numerical adjudication against glpsol and the CTMC is in
    test_solver_ba_qrf_path.py; this only pins the dispatch surface.
    """
    m, _ = _cqn(3)
    assert method in SolverBA(m).listValidMethods()


@pytest.mark.parametrize('method', ['qrf.bas.mmi', 'qrf.bas.mem'])
def test_qrf_bas_blocking_variants_unimplemented(method):
    """The BAS-blocking NLP variants have no native backend: api.mapqn ships the
    no-blocking qrf_noblo_* family only. Unlike the tokens above these were never
    merely withheld pending coverage."""
    m, _ = _cqn(3)
    assert method not in SolverBA(m).listValidMethods()
    with pytest.raises(NotImplementedError):
        SolverBA(m, method).getAvgTable()
