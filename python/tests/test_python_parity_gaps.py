"""Regression tests for the python/MATLAB parity gaps closed on 2026-08-14.

Every expected value here was READ OFF THE MATLAB REFERENCE on the same fixture
(through the matlab MCP server), not off another port. Where a value has an
independent analytic form, that form is asserted too, so a shared defect in both
codebases could not make the test pass.

See ../PYTHON-PARITY-GAPS.md and git show 449847e7b:_kb/log.md [2026-08-14].
"""

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, Source, Sink, ClosedClass,
                         OpenClass, Exp, Erlang, HyperExp, MarkedMAP, Prior,
                         SchedStrategy, SolverCTMC, SolverCTMCOptions, SolverMAM,
                         SolverAUTO, SolverAUTOOptions, SolverMVA)


# ---------------------------------------------------------------------------
# 1. Exact MAP/G/1/K family (api/qsys/mapg1k.py)
# ---------------------------------------------------------------------------

def test_qsys_mapg1k_reproduces_exact_mm1k():
    """M/M/1/5 with rho = 2/3: the analytic loss is 0.04812030."""
    from line_solver.api.qsys import qsys_mapg1k

    r = qsys_mapg1k([[-2.0]], [[2.0]], {'type': 'gamma', 'alpha': 1, 'theta': 1 / 3.0}, 5)

    rho = 2.0 / 3.0
    p = np.array([rho ** i for i in range(6)])
    p = p / p.sum()
    assert r['pK'] == pytest.approx(p[-1], abs=1e-14)
    assert r['p0'] == pytest.approx(p[0], abs=1e-14)
    assert r['meanQueueLength'] == pytest.approx(float(np.arange(6) @ p), abs=1e-13)
    assert r['plevel'].sum() == pytest.approx(1.0, abs=1e-12)


def test_qsys_mmapg1k_resolves_loss_by_class():
    """MATLAB qsys_mmapg1k on an asymmetric two-class MMAP, gamma service, K=6.

    The point of the exact branch is that the loss ratio DIFFERS by class; an
    aggregate-only analysis makes them equal by construction.
    """
    from line_solver.api.qsys import qsys_mmapg1k

    D0 = [[-3.0, 0.5], [0.2, -1.0]]
    D1a = [[1.5, 0.2], [0.1, 0.3]]
    D1b = [[0.6, 0.2], [0.2, 0.2]]
    r = qsys_mmapg1k(D0, [D1a, D1b], {'type': 'gamma', 'alpha': 2, 'theta': 0.25}, 6)

    assert r['throughput'][0] == pytest.approx(0.816666219297204, rel=1e-12)
    assert r['throughput'][1] == pytest.approx(0.517204158696982, rel=1e-12)
    assert r['lossRatio'][0] == pytest.approx(0.0550969363503425, rel=1e-11)
    assert r['lossRatio'][1] == pytest.approx(0.0472554971371378, rel=1e-11)
    assert r['lossRatio'][0] != pytest.approx(r['lossRatio'][1], rel=1e-3)
    assert r['pK'] == pytest.approx(0.0397252261466794, rel=1e-12)


def test_qsys_mapg1k_perflow_matches_matlab():
    """MATLAB qsys_mapg1k_perflow on three heterogeneous MAPs, det service."""
    from line_solver.api.qsys import qsys_mapg1k_perflow

    maps = [([[-2, 0], [0, -3]], [[1.5, 0.5], [1, 2]]),
            ([[-1.2]], [[1.2]]),
            ([[-5, 1], [0.5, -2]], [[3, 1], [0.5, 1]])]
    r = qsys_mapg1k_perflow(maps, {'type': 'det', 'd': 0.12}, 8)

    want = [2.32718248083633, 1.19728022318715, 2.32343182343482]
    np.testing.assert_allclose(r['throughput'], want, rtol=1e-11)
    assert r['lambdaAggregate'] == pytest.approx(5.86666666666667, rel=1e-12)


def _finite_cap_two_class():
    model = Network('mamfc')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue1', SchedStrategy.FCFS)
    queue.setNumberOfServers(1)
    sink = Sink(model, 'Sink')
    oc1 = OpenClass(model, 'C1')
    oc2 = OpenClass(model, 'C2')
    source.setArrival(oc1, Erlang.fitMeanAndOrder(1.0, 3))
    source.setArrival(oc2, HyperExp.fitMeanAndSCV(1.0, 8.0))
    queue.setService(oc1, Exp(2.5))
    queue.setService(oc2, Exp(2.5))
    queue.setCapacity(5)
    P = model.initRoutingMatrix()
    P[oc1] = Network.serialRouting(source, queue, sink)
    P[oc2] = Network.serialRouting(source, queue, sink)
    model.link(P)
    return model


def test_mam_finite_buffer_takes_the_exact_branch():
    """SolverMAM must resolve the finite-buffer loss BY CLASS.

    Both classes arrive at rate 1 but differ in variability, so equal
    throughputs would prove the truncate-and-renormalize fallback ran.
    """
    t = SolverMAM(_finite_cap_two_class()).getAvgTable()
    qlen = np.asarray(t['QLen'], dtype=float)
    tput = np.asarray(t['Tput'], dtype=float)

    assert qlen[2] == pytest.approx(0.95216, abs=5e-6)
    assert qlen[3] == pytest.approx(0.90961, abs=5e-6)
    assert tput[2] == pytest.approx(0.93009, abs=5e-6)
    assert tput[3] == pytest.approx(0.88853, abs=5e-6)
    assert abs(tput[2] - tput[3]) > 1e-3


def test_mam_finite_buffer_single_class_matches_matlab():
    """Erlang-3 arrivals into a single-server FCFS queue of capacity 6."""
    model = Network('fc1')
    s = Source(model, 'S')
    q = Queue(model, 'Q', SchedStrategy.FCFS)
    k = Sink(model, 'K')
    c = OpenClass(model, 'C')
    s.setArrival(c, Erlang.fitMeanAndOrder(1.0, 3))
    q.setService(c, Exp(2.0))
    q.setCapacity(6)
    model.link(Network.serialRouting(s, q, k))

    t = SolverMAM(model).getAvgTable()
    assert float(np.asarray(t['QLen'])[1]) == pytest.approx(0.7422, abs=5e-5)
    assert float(np.asarray(t['Tput'])[1]) == pytest.approx(0.99929, abs=5e-6)


# ---------------------------------------------------------------------------
# 2. Marked-MMAP CTMC state expansion (api/sn/proc_form.py)
# ---------------------------------------------------------------------------

def _marked_source_model():
    """Asymmetric modulating chain: pi = (0.75, 0.25).

    A SYMMETRIC D0 gives pi = (0.5, 0.5), where the old defect looked like an
    exact factor of two and hid its own mechanism.
    """
    D0 = np.array([[-4.0, 0.5], [1.5, -2.5]])
    D11 = np.array([[3.0, 0.0], [0.0, 0.8]])
    D12 = np.array([[0.5, 0.0], [0.0, 0.2]])
    model = Network('mm')
    src = Source(model, 'S')
    q = Queue(model, 'Q', SchedStrategy.FCFS)
    snk = Sink(model, 'K')
    r1 = OpenClass(model, 'R1')
    r2 = OpenClass(model, 'R2')
    src.setMarkedArrival(MarkedMAP([D0, D11, D12]), [r1, r2])
    q.setService(r1, Exp(10.0))
    q.setService(r2, Exp(10.0))
    P = model.initRoutingMatrix()
    P[r1] = Network.serialRouting(src, q, snk)
    P[r2] = Network.serialRouting(src, q, snk)
    model.link(P)
    return model, D0, D11, D12


def test_proc_to_map_reads_the_marked_m3a_layout():
    """A four-block {D0,D1,D11,D12} entry is a MAP of order 2, not of order 1."""
    from line_solver.api.sn.proc_form import proc_to_map, proc_n_phases

    _model, D0, D11, D12 = _marked_source_model()
    entry = [D0, D11 + D12, D11, D12]
    d0, d1 = proc_to_map(entry)
    assert d0 is not None and d0.shape == (2, 2)
    np.testing.assert_allclose(d0, D0)
    assert proc_n_phases(entry) == 2


def test_marked_mmap_ctmc_carries_the_carrier_phase():
    """The CTMC must expand the carrier's modulating phase.

    MATLAB gives a 38 x 9 space at cutoff=2 with Source throughput
    2.2988 / 0.42355; the analytic rates are pi*D1k*e = 2.45 / 0.425 and the
    truncated values approach them as the cutoff rises.
    """
    model, D0, D11, D12 = _marked_source_model()
    solver = SolverCTMC(model, SolverCTMCOptions(cutoff=2))
    space, _local = solver.getStateSpace()
    space = np.asarray(space)

    assert space.shape == (38, 9)

    t = solver.getAvgTable()
    tput = np.asarray(t['Tput'], dtype=float)
    assert tput[0] == pytest.approx(2.2988, abs=5e-5)
    assert tput[1] == pytest.approx(0.42355, abs=5e-6)

    # The unweighted phase sums the old defect returned were trace(D1k).
    assert tput[0] != pytest.approx(float(np.trace(D11)), abs=1e-2)
    assert tput[1] != pytest.approx(float(np.trace(D12)), abs=1e-2)

    # ... and the truncation approaches the analytic per-mark rates.
    pi = np.linalg.solve(np.vstack([(D0 + D11 + D12).T[:-1], np.ones(2)]), [0, 1])
    lam = [float(pi @ D11 @ np.ones(2)), float(pi @ D12 @ np.ones(2))]
    t4 = np.asarray(SolverCTMC(_marked_source_model()[0],
                               SolverCTMCOptions(cutoff=4)).getAvgTable()['Tput'],
                    dtype=float)
    assert abs(t4[0] - lam[0]) < abs(tput[0] - lam[0])


# ---------------------------------------------------------------------------
# 3. CTMC first passage and transient sensitivity
# ---------------------------------------------------------------------------

def _closed_two_station():
    model = Network('fp')
    d = Delay(model, 'Think')
    q = Queue(model, 'Q1', SchedStrategy.FCFS)
    cl = ClosedClass(model, 'C1', 3, d, 0)
    d.setService(cl, Exp(1.0))
    q.setService(cl, Exp(2.0))
    model.link(Network.serialRouting(d, q))
    return model


def test_ctmc_first_passage_moments_match_matlab():
    solver = SolverCTMC(_closed_two_station())
    space, _ = solver.getStateSpace()
    n = np.asarray(space).shape[0]

    m, mall = solver.getFirstPassTMoments(1, n, 3)
    np.testing.assert_allclose(np.atleast_1d(m), [2.5, 10.25, 60.75], rtol=1e-12)
    np.testing.assert_allclose(np.asarray(mall)[0, :], [2.5, 10.25, 60.75], rtol=1e-12)


def test_ctmc_first_passage_cdf_matches_matlab():
    solver = SolverCTMC(_closed_two_station())
    space, _ = solver.getStateSpace()
    n = np.asarray(space).shape[0]

    RD, out = solver.getCdfFirstPassT(1, n)
    assert RD.shape == (1000, 2)
    assert RD[99, 0] == pytest.approx(0.991432904835315, rel=1e-12)
    assert RD[99, 1] == pytest.approx(9.90990990990991, rel=1e-12)
    assert RD[-1, 0] == pytest.approx(1.0, abs=1e-12)
    assert np.all(np.diff(RD[:, 0]) >= -1e-12)   # a CDF cannot decrease


def test_ctmc_first_passage_refuses_an_empty_target():
    solver = SolverCTMC(_closed_two_station())
    with pytest.raises(ValueError):
        solver.getCdfFirstPassT(1, [])
    with pytest.raises(ValueError):
        solver.getFirstPassTMoments(1, [])


def test_ctmc_transient_sens_tracks_the_augmented_exponential():
    """The exact solution is [pi, dpi](t) = [pi0, 0] expm(A t) with
    A = [[Q, dQ], [0, Q]], which is what the augmented ODE integrates."""
    from scipy.linalg import expm
    from line_solver.api.mc import ctmc_transient_sens

    Q = np.array([[-2.0, 2, 0], [1, -3, 2], [0, 1, -1]])
    dQ = np.array([[-1.0, 1, 0], [0, 0, 0], [0, 0, 0]])
    dpi, pi, t = ctmc_transient_sens(Q, dQ, [1, 0, 0], 0, 2.0)

    n = 3
    A = np.zeros((2 * n, 2 * n))
    A[:n, :n] = Q
    A[:n, n:] = dQ
    A[n:, n:] = Q
    v0 = np.zeros(2 * n)
    v0[0] = 1.0
    v = v0 @ expm(A * float(t[-1]))

    # ode23/RK23 default tolerances are rtol 1e-3
    np.testing.assert_allclose(pi[-1], v[:n], atol=1e-3)
    np.testing.assert_allclose(dpi[-1], v[n:], atol=1e-3)
    np.testing.assert_allclose(pi.sum(axis=1), 1.0, atol=1e-6)
    # probability mass is conserved, so its sensitivity sums to zero
    np.testing.assert_allclose(dpi.sum(axis=1), 0.0, atol=1e-6)


# ---------------------------------------------------------------------------
# 4. retrieval_mva
# ---------------------------------------------------------------------------

def test_retrieval_mva_matches_the_normalizing_constant_route():
    from line_solver.api.retrieval import retrieval_mva, retrieval_metrics

    m = [2, 1]
    lam = [0.7, 0.4, 1.1, 0.25]
    eta = [[0.30, 0.20, 0.05], [0.10, 0.35, 0.15],
           [0.22, 0.08, 0.30], [0.40, 0.12, 0.02]]
    gam = [[1.0, 0.6], [0.5, 1.3], [0.9, 0.7], [1.4, 0.3]]

    pmiss, phit, pdh = retrieval_mva(m, lam, eta, gam)

    want_miss = [0.194737944942428, 0.17262232326757,
                 0.205666676275919, 0.154034874700374]
    np.testing.assert_allclose(pmiss, want_miss, rtol=1e-12)
    np.testing.assert_allclose(phit[0], [0.563251323387614, 0.216468479430536,
                                         0.461437283065759, 0.758842914116091],
                               rtol=1e-12)
    assert pdh.shape == (3, 4)

    ref = retrieval_metrics(m, lam, eta, gam)
    ref_miss = ref[0] if isinstance(ref, tuple) else ref
    np.testing.assert_allclose(pmiss, np.asarray(ref_miss).ravel(), rtol=1e-9)


def test_retrieval_mva_refuses_a_system_it_cannot_enumerate():
    """It is an exact oracle for SMALL systems; the memo doubles per item."""
    from line_solver.api.retrieval import retrieval_mva

    with pytest.raises(ValueError):
        retrieval_mva([2], np.ones(30), np.ones((30, 2)), np.ones((30, 1)),
                      max_items=20)


# ---------------------------------------------------------------------------
# 5. SolverAUTO method name resolution and the family routes
# ---------------------------------------------------------------------------

def _cqn():
    model = Network('cqn')
    d = Delay(model, 'D')
    q = Queue(model, 'Q', SchedStrategy.PS)
    c = ClosedClass(model, 'C', 4, d, 0)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    model.link(Network.serialRouting(d, q))
    return model


def test_auto_resolves_an_unqualified_algorithm_to_its_owning_family():
    kind, family, submethod = SolverAUTO(_cqn()).resolveMethodToken('comom')
    assert (kind, family, submethod) == ('family', 'nc', 'comom')


def test_auto_keeps_the_submethod_of_a_qualified_token():
    kind, family, submethod = SolverAUTO(_cqn()).resolveMethodToken('nc.comom')
    assert (kind, family, submethod) == ('family', 'nc', 'comom')


def test_auto_refuses_an_unknown_token_instead_of_running_the_heuristic():
    with pytest.raises(ValueError, match='[Uu]nrecognized method'):
        SolverAUTO(_cqn()).resolveMethodToken('nonsense_method')


@pytest.mark.parametrize('method_name,expected', [('comom', 'NC'), ('nc.exact', 'NC'),
                                            ('ctmc', 'CTMC'), ('fluid', 'FLUID')])
def test_auto_runs_the_family_the_token_names(method_name, expected):
    s = SolverAUTO(_cqn(), SolverAUTOOptions(selection_method=method_name, verbose=False))
    qlen = np.asarray(s.getAvgTable()['QLen'], dtype=float)
    assert s.getSelectedSolverName() == expected
    # every family solves the same product-form model
    np.testing.assert_allclose(qlen, [1.80952, 2.19048], atol=5e-4)


def test_auto_routes_a_prior_model_to_uq():
    model = Network('uq')
    d = Delay(model, 'D')
    q = Queue(model, 'Q', SchedStrategy.PS)
    c = ClosedClass(model, 'C', 3, d, 0)
    d.setService(c, Exp(1.0))
    q.setService(c, Prior([Exp(2.0), Exp(3.0)], [0.6, 0.4]))
    model.link(Network.serialRouting(d, q))

    s = SolverAUTO(model)
    assert s._pinned_family == 'uq'
    qlen = np.asarray(s.getAvgTable()['QLen'], dtype=float)
    np.testing.assert_allclose(qlen, [1.73198, 1.26802], atol=5e-5)

    post = s.getPosteriorTable()
    assert len(post) == 4          # two alternatives x two stations


def test_auto_lists_families_and_qualified_methods():
    v = SolverAUTO(_cqn()).listValidMethods()
    for method_name in ('default', 'exact', 'bound', 'mva', 'nc', 'ctmc', 'fluid'):
        assert method_name in v
    assert any(t.startswith('nc.') for t in v)


# ---------------------------------------------------------------------------
# 6. Network initialization wrappers
# ---------------------------------------------------------------------------

def test_init_from_marginal_and_running_sets_every_station():
    model = Network('m')
    d = Delay(model, 'D')
    q = Queue(model, 'Q', SchedStrategy.FCFS)
    q.setNumberOfServers(2)
    c = ClosedClass(model, 'C', 4, d, 0)
    d.setService(c, Exp(1.0))
    q.setService(c, Erlang.fitMeanAndOrder(0.5, 2))
    model.link(Network.serialRouting(d, q))

    model.initFromMarginalAndRunning([[2], [2]], [[2], [2]])

    for station in model.getStations():
        state = station.get_state()
        assert state is not None and np.asarray(state).size > 0

    # a Delay has infinite servers; int(inf) used to abort this path
    t = SolverCTMC(model).getAvgTable()
    np.testing.assert_allclose(np.asarray(t['QLen'], dtype=float),
                               [2.5472, 1.4528], atol=5e-5)


def test_init_from_avg_table_qlen_round_trips_a_solver_result():
    model = _cqn()
    table = SolverMVA(model).getAvgTable()
    model.initFromAvgTableQLen(table)
    assert model._state_marginal is not None
    np.testing.assert_allclose(np.asarray(model._state_marginal, dtype=float).sum(),
                               4.0, atol=1e-9)


def test_init_from_avg_qlen_gives_a_job_back_when_rounding_overshoots():
    """round([1.6, 2.4]) is [2, 2], one job more than the population."""
    model = _cqn()
    model.initFromAvgQLen(np.array([[1.6], [2.4]]))
    np.testing.assert_allclose(np.asarray(model._state_marginal, dtype=float).sum(),
                               4.0, atol=1e-9)


# ---------------------------------------------------------------------------
# 7. Numerical method tails
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('method', [1, 2, 3, 4])
def test_amap2_adjust_gamma_returns_a_feasible_triple(method):
    from line_solver.lib.m3a.amap2 import amap2_adjust_gamma, _nonlcon_theoretical

    M1, M2, M3, GAMMA = 1.0, 3.0, 10.0, -0.9
    M2a, M3a, GAMMAa = amap2_adjust_gamma(M1, M2, M3, GAMMA, method=method)
    assert np.max(_nonlcon_theoretical(M1, M2a, M3a, GAMMAa)) <= 1e-6


def test_amap2_adjust_gamma_method3_is_the_priority_cascade():
    """MATLAB: M2a = 3, M3a = 13.50135, GAMMAa = -0.50045."""
    from line_solver.lib.m3a.amap2 import amap2_adjust_gamma

    M2a, M3a, GAMMAa = amap2_adjust_gamma(1.0, 3.0, 10.0, -0.9, method=3)
    assert M2a == pytest.approx(3.0, rel=1e-12)
    assert M3a == pytest.approx(13.50135, rel=1e-6)
    assert GAMMAa == pytest.approx(-0.50045, rel=1e-6)


def test_amap2_adjust_gamma_rejects_an_unknown_method():
    from line_solver.lib.m3a.amap2 import amap2_adjust_gamma

    with pytest.raises(ValueError):
        amap2_adjust_gamma(1.0, 3.0, 10.0, 0.5, method=9)


@pytest.mark.parametrize('algor', ['FI', 'CR', 'NI', 'IS'])
def test_gim1_r_algorithms_agree_with_matlab(algor):
    from line_solver.lib.thirdparty.smc import gim1_r_dual

    A0 = np.array([[0.10, 0.02, 0.03], [0.04, 0.12, 0.01], [0.02, 0.03, 0.09]])
    A1 = np.array([[0.20, 0.05, 0.04], [0.06, 0.18, 0.05], [0.05, 0.06, 0.22]])
    A2 = np.array([[0.30, 0.10, 0.16], [0.24, 0.14, 0.16], [0.20, 0.13, 0.20]])
    want = np.array([[0.156989337725, 0.0498435210798, 0.0663593763983],
                     [0.0945664287859, 0.169316609178, 0.0483242092906],
                     [0.0601183181013, 0.0611037375552, 0.138373063365]])

    R = gim1_r_dual(np.hstack([A0, A1, A2]), 'A', algor)
    np.testing.assert_allclose(R, want, atol=1e-11)


def test_gim1_r_rejects_an_unknown_algorithm():
    from line_solver.lib.thirdparty.smc import gim1_r_dual

    A = np.hstack([np.full((2, 2), 0.1), np.full((2, 2), 0.2), np.full((2, 2), 0.2)])
    with pytest.raises(ValueError):
        gim1_r_dual(A, 'A', 'ZZ')
