"""Tests for the shortest-job-next (SJN/SJF) station of Kant (1992): pfqn_mvasjn
over the population lattice, pfqn_amvasjn through its Schweitzer fixed point, and
the SolverMVA dispatch.

Mirrors jar/src/test/java/jline/api/PfqnSjnTest.java and
line-test.git/test/testsAPI/test_pfqn_sjn.m one for one.

The mathematics is anchored to results that do not come from that paper: a
network with no SJN station must reproduce exact MVA and Bard-Schweitzer, a
single-job closed network cannot queue at all, Little's law must hold at every
population including where the utilization cap binds, and SJN must beat any
size-blind discipline on mean response time. The literal values are the MATLAB
and JAR outputs, so a divergence in any codebase shows up here.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_mva, pfqn_bs
from line_solver.api.pfqn.sjn import (SjnOptions, SjnStarvationError, pfqn_amvasjn,
                                       pfqn_mvasjn)

L = np.array([[0.125], [0.100], [0.050]])
Z = [1.0]
SJN0 = [0]


def _scv(cv2):
    return np.array([[cv2], [1.0], [1.0]])


def test_no_sjn_matches_exact_mva():
    """With no SJN station the recursion is the Reiser-Lavenberg one."""
    for n in (1, 4, 9):
        X, Q, U, C, _prof, _it = pfqn_mvasjn(L, [n], Z)
        # native pfqn_mva returns (XN, CN, QN, UN, RN, TN, AN), not the MATLAB order
        Xr, _Cr, Qr, _Ur, Rr = pfqn_mva(L, np.array([n]), np.array(Z))[:5]
        assert X[0] == pytest.approx(float(np.ravel(Xr)[0]), rel=1e-11)
        assert np.allclose(Q.ravel(), np.asarray(Qr, dtype=float).ravel(), rtol=1e-11)
        assert np.allclose(C.ravel(), np.asarray(Rr, dtype=float).ravel(), rtol=1e-11)


def test_no_sjn_matches_bard_schweitzer():
    """With no SJN station the fixed point is the Bard-Schweitzer one."""
    X, Q, U, C, _prof, _it = pfqn_amvasjn(L, [8], Z)
    Xr, Qr, Ur, Cr = pfqn_bs(L, np.array([8]), np.array(Z))[:4]
    assert X[0] == pytest.approx(float(np.ravel(Xr)[0]), abs=1e-5)
    assert np.allclose(Q.ravel(), np.asarray(Qr, dtype=float).ravel(), atol=1e-5)


def test_single_job_cannot_queue():
    """One job in a closed network never waits, whatever the size distribution."""
    for cv2 in (0.25, 1.0, 4.0):
        for solver in (pfqn_mvasjn, pfqn_amvasjn):
            X, Q, U, C, _prof, _it = solver(L, [1], [2.0], _scv(cv2), SJN0)
            assert np.allclose(C.ravel(), L.ravel(), rtol=1e-10)
            assert X[0] == pytest.approx(1.0 / (2.0 + float(np.sum(L))), rel=1e-10)


def test_population_conservation():
    """Little's law over the whole network, including where the cap binds.

    The cap acts on the waiting time and the throughput is derived from it, so no
    jobs are lost. Only the fixed point may cap: the lattice refuses the
    starvation regime instead, since the cap would rescale the profile its next
    population step reads back.
    """
    for cv2 in (0.5, 1.0, 2.0):
        for n in (2, 7, 15, 40):
            X, Q, U, C, _prof, _it = pfqn_amvasjn(L, [n], Z, _scv(cv2), SJN0)
            assert float(np.sum(Q)) + X[0] * Z[0] == pytest.approx(n, abs=1e-8)
            if n == 40:
                with pytest.raises(SjnStarvationError):
                    pfqn_mvasjn(L, [n], Z, _scv(cv2), SJN0)
                continue
            X, Q, U, C, _prof, _it = pfqn_mvasjn(L, [n], Z, _scv(cv2), SJN0)
            assert float(np.sum(Q)) + X[0] * Z[0] == pytest.approx(n, abs=1e-8)


def test_multiclass_population_conservation():
    """The same invariant per class in a two-class model, pooled and prioritised."""
    L2 = np.array([[0.125, 0.060], [0.100, 0.080]])
    N2 = [4, 3]
    Z2 = [1.0, 1.0]
    scv2 = np.array([[2.0, 1.0], [1.0, 1.0]])
    for opt in (None, SjnOptions(prio=[1, 2])):
        X, Q, U, C, _prof, _it = pfqn_mvasjn(L2, N2, Z2, scv2, SJN0, None, opt)
        assert np.allclose(np.sum(Q, axis=0) + X * np.asarray(Z2), N2, atol=1e-8)


def test_sjn_beats_size_blind_scheduling():
    """SJN minimises the mean response time among non-preemptive disciplines.

    The comparison is only meaningful at CV^2 = 1, where the product-form solution
    describes the SAME service distribution: product form is insensitive.
    """
    for n in (4, 6, 9):
        X, Q, U, C, _prof, _it = pfqn_mvasjn(L, [n], Z, _scv(1.0), SJN0)
        Xr, _Cr, _Qr, _Ur, Rr = pfqn_mva(L, np.array([n]), np.array(Z))[:5]
        assert X[0] > float(np.ravel(Xr)[0])
        assert C[0, 0] < float(np.asarray(Rr, dtype=float).ravel()[0])


def test_matlab_parity():
    """Literal values of the MATLAB and JAR routines at N = 6, Z = 1."""
    xref = {0.5: 4.298002, 1.0: 4.256655, 2.0: 4.172968, 4.0: 4.021977}
    rref = {0.5: 0.184871, 1.0: 0.199049, 2.0: 0.228550, 4.0: 0.284778}
    for cv2 in xref:
        X, Q, U, C, _prof, _it = pfqn_mvasjn(L, [6], Z, _scv(cv2), SJN0)
        assert X[0] == pytest.approx(xref[cv2], abs=1e-5)
        assert C[0, 0] == pytest.approx(rref[cv2], abs=1e-5)
    X, Q, U, C, _prof, _it = pfqn_amvasjn(L, [6], Z, _scv(2.0), SJN0)
    assert X[0] == pytest.approx(4.162429, abs=1e-5)
    assert C[0, 0] == pytest.approx(0.227866, abs=1e-5)


def test_utilization_cap_is_enforced():
    """A saturated SJN station has no solution to the open-form equation.

    The utilization law caps it in the FIXED POINT, which re-derives the profile
    from the capped state and converges, so it returns instead of failing. The
    lattice cannot: the same factor would rescale the conditional waiting time
    profile that its next population step reads back, the correction compounding
    along the lattice, so it raises 'starvation' and the analyzer re-solves with
    the fixed point.
    """
    for n in (24, 40, 80):
        X, Q, U, C, _prof, _it = pfqn_amvasjn(L, [n], Z, _scv(2.0), SJN0)
        assert U[0, 0] <= 0.999 + 1e-9
        with pytest.raises(SjnStarvationError):
            pfqn_mvasjn(L, [n], Z, _scv(2.0), SJN0)
    X, Q, U, C, _prof, _it = pfqn_amvasjn(L, [24], Z, _scv(2.0), SJN0)
    assert X[0] == pytest.approx(7.992000, abs=1e-5)
    assert U[0, 0] == pytest.approx(0.999000, abs=1e-8)
    X, Q, U, C, _prof, _it = pfqn_amvasjn(L, [40], Z, _scv(2.0), SJN0, None,
                                          SjnOptions(umax=0.9))
    assert U[0, 0] <= 0.9 + 1e-9


def test_fixed_point_tracks_the_lattice():
    """Away from saturation the Schweitzer closure must stay close to the lattice."""
    for cv2 in (0.5, 1.0, 2.0, 4.0):
        for n in (3, 6):
            Xe, Qe, Ue, Ce, _p, _i = pfqn_mvasjn(L, [n], Z, _scv(cv2), SJN0)
            Xf, Qf, Uf, Cf, _p, _i = pfqn_amvasjn(L, [n], Z, _scv(cv2), SJN0)
            assert Xf[0] == pytest.approx(Xe[0], rel=0.01)
            assert Cf[0, 0] == pytest.approx(Ce[0, 0], rel=0.01)


def test_grid_refinement_converges():
    """The quadrature must be converged at the default grid."""
    X32, _q, _u, C32, _p, _i = pfqn_mvasjn(L, [6], Z, _scv(2.0), SJN0, None, SjnOptions(ns=32))
    X64, _q, _u, C64, _p, _i = pfqn_mvasjn(L, [6], Z, _scv(2.0), SJN0, None, SjnOptions(ns=64))
    assert X64[0] == pytest.approx(X32[0], rel=1e-4)
    assert C64[0, 0] == pytest.approx(C32[0, 0], rel=1e-4)


def test_odd_grid_rejected():
    """Composite Simpson integrates over panels of two subdivisions."""
    with pytest.raises(ValueError):
        pfqn_mvasjn(L, [4], Z, _scv(1.0), SJN0, None, SjnOptions(ns=31))


def test_solver_mva_dispatch():
    """SolverMVA must route a closed model with an SJF queue to the SJN analyzer."""
    from line_solver import (ClosedClass, Delay, Exp, HyperExp, Network, Queue,
                             SchedStrategy, SolverMVA)
    model = Network('sjn_dispatch')
    delay = Delay(model, 'Think')
    q1 = Queue(model, 'SJN', SchedStrategy.SJF)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    cl = ClosedClass(model, 'C1', 6, delay)
    delay.setService(cl, Exp(1.0))
    q1.setService(cl, HyperExp.fitMeanAndSCVBalanced(0.125, 2.0))
    q2.setService(cl, Exp(1 / 0.100))
    model.link(Network.serialRouting(delay, q1, q2))
    Qs, Us, Rs, Ts = SolverMVA(model).getAvg()[:4]
    X, Q, U, C, _p, _i = pfqn_mvasjn(np.array([[0.125], [0.100]]), [6], Z,
                                     np.array([[2.0], [1.0]]), SJN0)
    assert float(np.ravel(Ts)[1]) == pytest.approx(X[0], rel=1e-8)
    assert float(np.ravel(Qs)[1]) == pytest.approx(Q[0, 0], rel=1e-8)
    assert float(np.ravel(Us)[1]) == pytest.approx(U[0, 0], rel=1e-8)
    assert float(np.ravel(Rs)[1]) == pytest.approx(C[0, 0], rel=1e-8)
    assert float(np.sum(Qs)) == pytest.approx(6.0, rel=1e-8)


def _dispatch_model():
    """Closed model with one SJF queue, shared by the dispatch tests."""
    from line_solver import (ClosedClass, Delay, Exp, HyperExp, Network, Queue,
                             SchedStrategy)
    model = Network('sjn_dispatch')
    delay = Delay(model, 'Think')
    q1 = Queue(model, 'SJN', SchedStrategy.SJF)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    cl = ClosedClass(model, 'C1', 6, delay)
    delay.setService(cl, Exp(1.0))
    q1.setService(cl, HyperExp.fitMeanAndSCVBalanced(0.125, 2.0))
    q2.setService(cl, Exp(1 / 0.100))
    model.link(Network.serialRouting(delay, q1, q2))
    return model, np.array([[0.125], [0.100]]), np.array([[2.0], [1.0]])


def test_solver_mva_amva_method():
    """Asking for 'amva' must reach the fixed point, not the lattice."""
    from line_solver import SolverMVA
    model, Ld, scvd = _dispatch_model()
    solver = SolverMVA(model, 'amva')
    Qs, Us, Rs, Ts = solver.getAvg()[:4]
    # the analyzer hands options.iter_tol to the fixed point, so the direct call must use it
    opt = SjnOptions(tol=solver.options.iter_tol, iter_max=solver.options.max_iter)
    X, Q, U, C, _p, _i = pfqn_amvasjn(Ld, [6], Z, scvd, SJN0, None, opt)
    assert float(np.ravel(Ts)[1]) == pytest.approx(X[0], rel=1e-8)
    assert float(np.ravel(Qs)[1]) == pytest.approx(Q[0, 0], rel=1e-8)
    assert float(np.ravel(Us)[1]) == pytest.approx(U[0, 0], rel=1e-8)
    assert float(np.ravel(Rs)[1]) == pytest.approx(C[0, 0], rel=1e-8)
    # the two routes must differ, or the test would pass on the lattice
    Xe, _q, _u, Ce, _p, _i = pfqn_mvasjn(Ld, [6], Z, scvd, SJN0)
    assert X[0] != Xe[0]
    assert C[0, 0] != Ce[0, 0]


def test_solver_mva_lattice_threshold():
    """Past config.sjn_lattice_max the default dispatch switches to the fixed point."""
    from line_solver import SolverMVA
    model, Ld, scvd = _dispatch_model()
    solver = SolverMVA(model)
    solver.options.config = {'sjn_lattice_max': 1e5}
    _q, _u, Rbig, Tbig = solver.getAvg()[:4]
    Xe, _q2, _u2, Ce, _p, _i = pfqn_mvasjn(Ld, [6], Z, scvd, SJN0)
    assert float(np.ravel(Tbig)[1]) == pytest.approx(Xe[0], rel=1e-8)
    assert float(np.ravel(Rbig)[1]) == pytest.approx(Ce[0, 0], rel=1e-8)
    solver2 = SolverMVA(model)
    solver2.options.config = {'sjn_lattice_max': 2}   # the 7-state lattice no longer fits
    _q3, _u3, Rsmall, Tsmall = solver2.getAvg()[:4]
    opt2 = SjnOptions(tol=solver2.options.iter_tol, iter_max=solver2.options.max_iter)
    Xf, _q4, _u4, Cf, _p, _i = pfqn_amvasjn(Ld, [6], Z, scvd, SJN0, None, opt2)
    assert float(np.ravel(Tsmall)[1]) == pytest.approx(Xf[0], rel=1e-8)
    assert float(np.ravel(Rsmall)[1]) == pytest.approx(Cf[0, 0], rel=1e-8)
    assert float(np.ravel(Tsmall)[1]) != float(np.ravel(Tbig)[1])


def _multiclass_chain_model():
    """Two classes in one chain, both traversing every station and switching at Q2.

    Disjoint class routing would make ST(i,k) equal STchain(i,c) at every station
    with nonzero alpha, so the two deaggregation routes would coincide and the
    model could not tell them apart.
    """
    from line_solver import (ClosedClass, Delay, Exp, Network, Queue, SchedStrategy)
    model = Network('sjn_multiclass_chain')
    think = Delay(model, 'Think')
    q1 = Queue(model, 'Q1', SchedStrategy.SJF)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    q1.setNumberOfServers(1)
    q2.setNumberOfServers(1)
    ca = ClosedClass(model, 'A', 4, think, 0)
    cb = ClosedClass(model, 'B', 0, think, 0)
    think.setService(ca, Exp.fitMean(1.0))
    think.setService(cb, Exp.fitMean(2.0))
    q1.setService(ca, Exp.fitMean(0.5))
    q1.setService(cb, Exp.fitMean(0.9))
    q2.setService(ca, Exp.fitMean(0.7))
    q2.setService(cb, Exp.fitMean(1.25))
    P = model.init_routing_matrix()
    P.set(ca, ca, think, q1, 1.0)
    P.set(ca, ca, q1, q2, 1.0)
    P.set(ca, cb, q2, think, 1.0)
    P.set(cb, cb, think, q1, 1.0)
    P.set(cb, cb, q1, q2, 1.0)
    P.set(cb, ca, q2, think, 1.0)
    model.link(P)
    return model


def test_multiclass_chain_deaggregation():
    """On a multiclass chain Q and U must be rebuilt from Rchain, not scaled from Qchain.

    solver_mva_sjn_analyzer.m:138 passes [] in the Qchain and Uchain slots, so
    sn_deaggregate_chain_results rebuilds them per class. Passing the chain
    matrices instead scales by alpha alone and drops the ST(i,k)/STchain(i,c)
    weighting, which collapses both classes onto the chain average.

    The literals are the MATLAB ground truth for this model, which has
    STchain = [1.5; 0.7; 0.975] against class demands [1.0 2.0], [0.5 0.9] and
    [0.7 1.25], so no station has ST(i,k) equal to STchain(i,c).
    """
    from line_solver import SolverMVA
    solver = SolverMVA(_multiclass_chain_model(), 'sjn.mva')
    Q, U, R, T = solver.getAvg()[:4]
    Q = np.asarray(Q, dtype=float)
    R = np.asarray(R, dtype=float)
    # a delay station holds each class for its own mean service time, never the
    # chain average, so this alone separates the two deaggregation routes
    assert R[0, 0] == pytest.approx(1.0, rel=1e-9)
    assert R[0, 1] == pytest.approx(2.0, rel=1e-9)
    assert Q[0, 0] == pytest.approx(0.4161596548, rel=1e-8)
    assert Q[0, 1] == pytest.approx(0.8323193095, rel=1e-8)
    assert Q[1, 0] == pytest.approx(0.3480101577, rel=1e-8)
    assert Q[1, 1] == pytest.approx(0.6264182838, rel=1e-8)
    assert Q[2, 0] == pytest.approx(0.6379306748, rel=1e-8)
    assert Q[2, 1] == pytest.approx(1.1391619194, rel=1e-8)
    assert float(np.sum(Q)) == pytest.approx(4.0, rel=1e-8)
    # the classes must not share a queue length, which is what the chain-matrix route gives
    assert Q[1, 0] != pytest.approx(Q[1, 1], rel=1e-6)


def test_infinite_server_queue_is_refused():
    """Only INF makes a delay station; a PS queue with many servers is refused.

    solver_mva_sjn_analyzer.m:32-52 switches on sched(ist) alone, so an
    infinite-server PS station reaches the nservers ~= 1 refusal at :46 instead
    of being folded into the think time.
    """
    from line_solver import (ClosedClass, Delay, Exp, Network, Queue, SchedStrategy,
                             SolverMVA)
    model = Network('sjn_inf_ps')
    think = Delay(model, 'Think')
    q1 = Queue(model, 'Q1', SchedStrategy.SJF)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    q1.setNumberOfServers(1)
    q2.setNumberOfServers(float('inf'))
    cl = ClosedClass(model, 'A', 3, think, 0)
    think.setService(cl, Exp.fitMean(1.0))
    q1.setService(cl, Exp.fitMean(0.5))
    q2.setService(cl, Exp.fitMean(0.7))
    model.link(Network.serialRouting(think, q1, q2))
    with pytest.raises(Exception, match='servers'):
        SolverMVA(model, 'sjn.mva').getAvg()


def test_empty_chain_metrics_are_zero():
    """A chain with no jobs carries no load, per solver_mva_sjn_analyzer.m:131-135.

    This is a regression guard, not evidence for the chain-level zeroing: it
    passes without it too, because an empty chain has Xchain == 0 and the final
    non-finite sweep in sn_deaggregate_chain_results already flattens the NaN
    that Qchain / Tchain leaves behind. The analyzer still zeroes explicitly, to
    match MATLAB and to keep a NaN out of the deaggregation input.
    """
    from line_solver import (ClosedClass, Delay, Exp, Network, Queue, SchedStrategy,
                             SolverMVA)
    model = Network('sjn_empty_chain')
    think = Delay(model, 'Think')
    q1 = Queue(model, 'Q1', SchedStrategy.SJF)
    q1.setNumberOfServers(1)
    ca = ClosedClass(model, 'A', 3, think, 0)
    cb = ClosedClass(model, 'B', 0, think, 0)
    think.setService(ca, Exp.fitMean(1.0))
    think.setService(cb, Exp.fitMean(2.0))
    q1.setService(ca, Exp.fitMean(0.5))
    q1.setService(cb, Exp.fitMean(0.9))
    P = model.init_routing_matrix()
    P.set(ca, ca, think, q1, 1.0)
    P.set(ca, ca, q1, think, 1.0)
    P.set(cb, cb, think, q1, 1.0)
    P.set(cb, cb, q1, think, 1.0)
    model.link(P)
    Q, U, R, T = SolverMVA(model, 'sjn.mva').getAvg()[:4]
    Q = np.asarray(Q, dtype=float)
    U = np.asarray(U, dtype=float)
    T = np.asarray(T, dtype=float)
    assert np.all(Q[:, 1] == 0.0)
    assert np.all(U[:, 1] == 0.0)
    assert np.all(T[:, 1] == 0.0)
    assert float(np.sum(Q[:, 0])) == pytest.approx(3.0, rel=1e-8)


def test_solver_mva_rejects_open_sjn():
    """The population recursion has no open counterpart here."""
    from line_solver import (Exp, Network, OpenClass, Queue, SchedStrategy, Sink,
                             SolverMVA, Source)
    model = Network('sjn_open')
    source = Source(model, 'Source')
    q1 = Queue(model, 'SJN', SchedStrategy.SJF)
    sink = Sink(model, 'Sink')
    cl = OpenClass(model, 'C1')
    source.setArrival(cl, Exp(1.0))
    q1.setService(cl, Exp(4.0))
    model.link(Network.serialRouting(source, q1, sink))
    with pytest.raises(Exception):
        SolverMVA(model).getAvg()


def test_citation_is_registered():
    """A method that reaches a user without a citation entry loses attribution."""
    from line_solver.api.io.citations import citations_for
    e = citations_for(['sjn.mva', 'sjn.amva'])
    assert len(e) == 1
    assert e[0]['key'] == 'Kant92'
    assert len(citations_for(['sjf'])) == 1
