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


def _cqn_nodelay(N):
    """Three FCFS queues in a cycle, no delay: what harel requires."""
    m = Network('cqnB')
    q1 = Queue(m, 'Queue1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Queue2', SchedStrategy.FCFS)
    q3 = Queue(m, 'Queue3', SchedStrategy.FCFS)
    c = ClosedClass(m, 'C1', N, q1)
    q1.setService(c, Exp(2.0)); q2.setService(c, Exp(3.0)); q3.setService(c, Exp(5.0))
    m.link(Network.serialRouting(q1, q2, q3))
    return m


@pytest.mark.parametrize('N', [3, 5, 10, 20])
def test_harel_bounds_bracket_exact(N):
    m = _cqn_nodelay(N)
    exX = float(SolverMVA(m, 'exact').getAvgTable().Tput.iloc[0])
    xl = float(np.atleast_1d(SolverBA(m, 'harel.lower').getAvgSysTput()).ravel()[0])
    xu = float(np.atleast_1d(SolverBA(m, 'harel.upper').getAvgSysTput()).ravel()[0])
    assert xl <= exX * (1 + 1e-9)
    assert xu >= exX * (1 - 1e-9)


def test_harel_upper_is_tighter_than_sb_at_moderate_population():
    # the two families share the paper and the lower side; the harel upper
    # extrapolates from the EXACT normalizing constant, sb from power sums
    m = _cqn_nodelay(5)
    xu = float(np.atleast_1d(SolverBA(m, 'harel.upper').getAvgSysTput()).ravel()[0])
    su = float(np.atleast_1d(SolverBA(m, 'sb.upper').getAvgSysTput()).ravel()[0])
    exX = float(SolverMVA(m, 'exact').getAvgTable().Tput.iloc[0])
    assert xu <= su
    assert xu == pytest.approx(exX, rel=1e-9)   # at N = 5 <= 7 the point is N itself


def test_harel_refuses_a_delay_station():
    m, _ = _cqn()
    with pytest.raises(Exception):
        SolverBA(m, 'harel.upper').getAvgSysTput()


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


@pytest.mark.parametrize('meth', ['qr', 'qrf.mmi', 'qrf.mem', 'qrf.bethe'])
def test_qrf_alpha_free_arms_reject_delay_stations(meth):
    # The ALPHA-FREE arms build a population-free q, so they model every station
    # as a single server and an infinite server would silently be a different
    # model. _cqn carries a Delay, hence the rejection -- which now names the
    # two arms that DO serve it rather than the construct that rules it out.
    # 'qrf.mmi.ld' and 'qrf.mmi.linear' moved to the test below.
    m, _ = _cqn()
    with pytest.raises(ValueError, match="Use 'qrf.mmi.ld' or 'qrf.mmi.linear'"):
        SolverBA(m, method=meth).runAnalyzer()


@pytest.mark.parametrize('meth', ['qrf.mmi.ld', 'qrf.mmi.linear'])
def test_qrf_ld_arms_serve_delay_stations(meth):
    # alpha(i,n) = n IS the rate law of an infinite server, so the two
    # load-dependent arms answer the model's own chain rather than refusing it.
    # Small N here because this actually solves the NLP; the numerical oracle
    # against exact CTMC lives in test_solver_ba_qrf_ld.py.
    m = Network('cqnZ')
    d = Delay(m, 'Z')
    q1 = Queue(m, 'Q1', SchedStrategy.FCFS)
    c = ClosedClass(m, 'C', 2, d, 0)
    d.setService(c, Exp(1.0))
    q1.setService(c, Exp(1.5))
    m.link(Network.serialRouting(d, q1))
    QN = np.asarray(SolverBA(m, method=meth).getAvgQLen(), dtype=float).ravel()
    assert np.all(np.isfinite(QN)) and abs(QN.sum() - 2.0) < 1e-6


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
    # NO DELAY: the QR/LR/QRF reduction models every station as a single
    # server, so an INF station rules the whole family out of the advertised
    # list, in SolverBA.m and SolverBA.java as here. _cqn() carries a Delay.
    m = _cqn_nodelay(3)
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
    """The LP-backed QRF method names dispatch and return finite, valid bounds.

    Widened from a single symmetric 2-station chain to M=3, asymmetric rates and
    a phase-type station. Unlike the NLP tokens these optimize a face of the
    polytope, so the upper-bound direction IS guaranteed and is asserted.

    The blocking parameters are passed explicitly because both tokens now
    REFUSE a model without them (they used to fabricate exactly the values
    below). These are the no-blocking ones -- one configuration, no blocked
    station, capacity N -- which is the parameterisation the glpsol reference
    values were computed under.
    """
    M = len(_QRF_LP_CHAINS[name][1])
    N = _QRF_LP_CHAINS[name][0]
    qrf_params = {'f': 1, 'MR': 1, 'ZM': 0, 'ZZ': [0],
                  'MM': np.zeros((1, 2)), 'MM1': np.zeros((1, M)),
                  'BB': np.zeros((1, M)), 'F': [N] * M}
    solver = SolverBA(_qrf_chain(name), method,
                      config={'qrf_params': qrf_params})
    U = np.asarray(solver.getAvgTable().Util, dtype=float).flatten()
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
    # NO DELAY: an INF station rules the reduction family out; see the LP twin.
    m = _cqn_nodelay(3)
    assert method in SolverBA(m).listValidMethods()


@pytest.mark.parametrize('method', ['qrf.bas.mmi', 'qrf.bas.mem'])
def test_qrf_bas_blocking_variants_require_parameters(method):
    """The BAS-blocking NLP variants are native since 2026-07-29.

    They optimize MMI/MEM over the polytope qr_bounds_bas builds, the same one
    the LP token qrf.bas is validated on, so they carry the same contract: the
    blocking parameters are REQUIRED, there being no meaningful default for a
    bound on a blocking network. Numerical adjudication is in
    test_qrf_bas_nlp_matches_exact_two_station below.
    """
    m, _ = _cqn(3)
    with pytest.raises((ValueError, NotImplementedError)):
        SolverBA(m, method).getAvgTable()


def test_qrf_bas_nlp_matches_exact_two_station():
    """On a 2-station cycle at N=2 the pairwise joint is determined by the
    marginal, so the polytope contains the CTMC point and the bound is tight
    there: p(n1) = 2^n1/7 gives U = (6/7, 3/7) and QN = (10/7, 4/7)."""
    import numpy as np
    from line_solver.api.mapqn.parameters import QRBoundsBasParameters
    from line_solver.api.mapqn.qrf_bas_nlp import qrf_bas_mem, qrf_bas_mmi

    params = QRBoundsBasParameters(
        _M=2, _N=2, MR=1, f=1,
        K=np.array([1, 1]), F=np.array([2, 2]),
        MM=np.zeros((1, 2)), MM1=np.zeros((1, 2)),
        ZZ=np.array([0]), BB=np.zeros((1, 2)),
        mu=[np.array([[1.0]]), np.array([[2.0]])],
        v=[np.zeros((1, 1)), np.zeros((1, 1))],
        r=np.array([[0.0, 1.0], [1.0, 0.0]]),
    )
    for fn in (qrf_bas_mmi, qrf_bas_mem):
        UN, QN = fn(params)
        np.testing.assert_allclose(UN, [6 / 7, 3 / 7], atol=1e-6)
        np.testing.assert_allclose(QN, [10 / 7, 4 / 7], atol=1e-6)


def test_qrf_bas_mmi_spans_the_idle_cells():
    """The MI block runs from n = 0, the AMPL source's range and MATLAB's.

    This port used 1..F until 2026-09-03, which dropped the idle cells. Since
    U_i = 1 - p_i(0) that excluded the strongest correlation in a closed chain
    from the functional meant to measure coupling, and it was the last thing
    keeping 'qrf.bas.mmi' refused in MATLAB and C++: the token could not be
    enabled while the four codebases disagreed on what it computes. Asserted on
    the term set rather than on a value, because the fixtures whose polytope
    pins U would pass either way.
    """
    import numpy as np
    from line_solver.api.mapqn.parameters import QRBoundsBasParameters
    from line_solver.api.mapqn.qrf_bas_nlp import (
        _bas_polytope, _mmi_terms, _p2_columns)

    params = QRBoundsBasParameters(
        _M=2, _N=2, MR=1, f=1,
        K=np.array([1, 1]), F=np.array([2, 2]),
        MM=np.zeros((1, 2)), MM1=np.zeros((1, 2)),
        ZZ=np.array([0]), BB=np.zeros((1, 2)),
        mu=[np.array([[1.0]]), np.array([[2.0]])],
        v=[np.zeros((1, 1)), np.zeros((1, 1))],
        r=np.array([[0.0, 1.0], [1.0, 0.0]]),
    )
    _, _, _, _, model = _bas_polytope(params)
    cols = _p2_columns(model, params)

    default, _, _ = _mmi_terms(cols, params)
    from_zero, _, _ = _mmi_terms(cols, params, n_from=0)
    from_one, _, _ = _mmi_terms(cols, params, n_from=1)
    assert default.size == from_zero.size
    assert from_zero.size > from_one.size


def test_scb_brackets_the_multiclass_system_not_this_model():
    """scb is the only SolverBA family that does not bracket the given model.

    Its lower side IS the model's exact single-class throughput; what the pair
    brackets is the multiclass system the single-class model aggregates
    (Dowdy et al. 1992). Asserting the usual "bracket the exact answer"
    property would therefore be asserting the wrong thing.
    """
    from line_solver.api.pfqn import pfqn_scbgap
    m = _cqn_nodelay(3)                     # scb needs a delay-free model
    exX = float(SolverMVA(m, 'exact').getAvgTable().Tput.iloc[0])
    b = SolverBA(m, method='scb.upper').get_bounds()
    xl = float(np.asarray(b['Tlower'])[0, 0])
    xu = float(np.asarray(b['Tupper'])[0, 0])
    assert abs(xl - exX) < 1e-9             # exact, not a bound on this model
    assert xu >= xl
    # the gap never exceeds Theorem 3 at N = 3 queueing stations, K = 3
    gap = pfqn_scbgap(3, 3)
    assert xu / xl <= 1.0 + gap / (1.0 - gap) + 1e-9
    # Corollary 1: one utilization ratio for every device, all capped at 1
    lo = SolverBA(m, method='scb.lower'); lo.runAnalyzer()
    hi = SolverBA(m, method='scb.upper'); hi.runAnalyzer()
    Ulo = np.asarray(lo._result['UN'])[:, 0]
    Uhi = np.asarray(hi._result['UN'])[:, 0]
    assert np.max(np.abs(Uhi / Ulo - xu / xl)) < 1e-9
    assert np.all(Uhi <= 1.0 + 1e-9)


def test_scb_refuses_delay_multiserver_and_multiclass():
    m, _ = _cqn()                           # has a delay station
    with pytest.raises(Exception):
        SolverBA(m, method='scb.lower').runAnalyzer()
    # multiclass: the aggregate is the INPUT to scb, not its output
    mc = Network('mc')
    d = Delay(mc, 'Z'); q = Queue(mc, 'Q', SchedStrategy.PS)
    c1 = ClosedClass(mc, 'A', 2, d, 0); c2 = ClosedClass(mc, 'B', 2, d, 0)
    d.setService(c1, Exp(1)); d.setService(c2, Exp(1))
    q.setService(c1, Exp(2)); q.setService(c2, Exp(1.5))
    P = mc.initRoutingMatrix()
    P[(c1, c1)] = Network.serialRouting(d, q)
    P[(c2, c2)] = Network.serialRouting(d, q)
    mc.link(P)
    with pytest.raises(Exception):
        SolverBA(mc, method='scb.upper').runAnalyzer()


def test_scb_is_not_an_auto_candidate():
    """auto.lower must stay a bound on THIS model, so scb cannot enter it.

    scb.lower equals the exact throughput, so admitting it would collapse the
    composite lower bound onto the exact answer and report it as a bound.
    """
    m = _cqn_nodelay(3)
    exX = float(SolverMVA(m, 'exact').getAvgTable().Tput.iloc[0])
    b = SolverBA(m, method='auto.lower').get_bounds()
    assert float(np.asarray(b['Tlower'])[0, 0]) < exX


def _cqn_bas_blocking():
    """The cqn_bas_blocking example.

    Queue1 declares the BAS drop rule and Queue2 holds 1 job against a
    population of 2, so its buffer BINDS.
    """
    from line_solver import DropStrategy
    m = Network('cqn_bas_blocking')
    q1 = Queue(m, 'Queue1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Queue2', SchedStrategy.FCFS)
    c = ClosedClass(m, 'Class1', 2, q1, 0)
    q1.setService(c, Exp(1.0))
    q2.setService(c, Exp(0.8))
    q2.setCap(1)
    q1.setDropRule(c, DropStrategy.BAS)
    m.link(Network.serialRouting(q1, q2))
    return m


@pytest.mark.parametrize('method', ['gb.upper', 'gb.lower', 'aba.upper'])
def test_blocking_is_refused_not_bounded(method):
    """A blocking-blind family asked for BY NAME must be REFUSED on a model
    with a binding buffer.

    Before the gate, gb.upper on this model returned QLen 1.28 at Queue2 -- a
    station that can never hold more than one job -- because every family here
    is parameterized by demands and a population alone and never read the cap
    or the drop rule.

    'default' and 'auto.upper' are NOT in this list any more: since
    sn_to_qrf_blocking derives the QRF BAS tables from the model, they resolve
    to 'qrf.bas' on a blocked model instead of to the blind 'gb.upper'. That
    routing is asserted by test_blocking_routes_the_default_to_qrf_bas below.
    Asking for a blind family by name is still refused, which is what keeps the
    AUTO composite from probing its way past the gate: 'auto.lower' has no
    blocking counterpart, since the analyzer solves qrf.bas in the 'max'
    direction alone.
    """
    with pytest.raises(Exception):
        SolverBA(_cqn_bas_blocking(), method=method).getAvgTable()


@pytest.mark.parametrize('method', ['default', 'auto', 'auto.upper'])
def test_blocking_routes_the_default_to_qrf_bas(method):
    """On a blocked model the upper-side aliases MEAN the QRF BAS bound.

    They used to be refused with the blind families, which left
    cqn_bas_blocking with no usable SolverBA method at all: the advice was
    'qrf.bas', and that in turn demanded hand-built blocking tables. The tables
    are derived now, so the alias resolves rather than refusing -- the same
    routing SolverMVA performs to reach 'sqd'.
    """
    routed = np.asarray(SolverBA(_cqn_bas_blocking(), method=method).getAvgUtil()).ravel()
    named = np.asarray(SolverBA(_cqn_bas_blocking(), method='qrf.bas').getAvgUtil()).ravel()
    assert np.allclose(routed, named)


def test_blocking_still_refuses_the_lower_side_of_auto():
    """qrf.bas is UPPER-only, so 'auto.lower' has nothing to route to."""
    with pytest.raises(Exception):
        SolverBA(_cqn_bas_blocking(), method='auto.lower').getAvgTable()


def test_blocking_narrows_the_method_list_to_the_qrf_blocking_bounds():
    s = SolverBA(_cqn_bas_blocking())
    valid = s.list_valid_methods()
    assert valid, 'the QRF blocking bounds stay available on this model'
    # 'default' is offered back because it now MEANS 'qrf.bas' on this model;
    # it is the one entry whose RESOLVED form ('gb.upper') is blind, which is
    # exactly why runAnalyzer rewrites it rather than dispatching it.
    assert 'default' in valid
    assert all(m.startswith('qrf.bas') or m.startswith('qrf.rsrd')
               for m in valid if m != 'default')
    with pytest.raises(Exception):
        SolverBA(_cqn_bas_blocking(), method='gb.upper').get_bounds()


def test_blocking_get_bounds_refuses_the_routed_default_as_one_sided():
    """getBounds must not contradict the run that just succeeded.

    'default' resolves to qrf.bas here, which the analyzer solves in the 'max'
    direction alone, so there is no bracket -- but the reason is one-sidedness,
    not blindness to the buffer.
    """
    with pytest.raises(Exception) as e:
        SolverBA(_cqn_bas_blocking()).get_bounds()
    assert 'UPPER-only' in str(e.value)


def test_the_blocking_gate_leaves_an_unblocked_model_alone():
    m = _cqn_nodelay(3)
    assert SolverBA(m, method='gb.upper').getAvgTable() is not None
    assert len(SolverBA(m).list_valid_methods()) > 10
