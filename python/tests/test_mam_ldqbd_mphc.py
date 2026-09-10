"""
Native-Python test for the EXACT multiserver phase-type LD-QBD chain.

`solver_mam_ldqbd` used to collapse the c parallel PH servers into one PH
process run at min(n,c) times its speed. That keeps the aggregate service rate
right but forgets which phase each busy server is in, and it cost ~1e-2 relative
against SolverCTMC on Erlang and HyperExp service. Since 2026-08-18 the level
carries the MULTISET of the phases the min(n,c) busy servers sit in
(``api/mam/ldqbd_mphc.py``), which is exact.

The oracle is SolverCTMC, exact for phase-type service on these models, and the
bar is machine precision rather than a percentage: an exact chain against an
exact chain has nothing left to differ by. The c == 1 and exponential cases are
here too, because the multiset construction has to reproduce them -- at c == 1
the multiset IS the phase, and at one phase there is no configuration at all.
"""

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, ClosedClass, Exp, Erlang,
                         HyperExp, Coxian, SchedStrategy, SolverMAM, SolverCTMC)
from line_solver.api.mam import ldqbd_mphc, ph_multisets

TOL = 1e-9          # exact vs exact: only floating-point noise may remain
QI = 1              # station rows: 0 = Delay, 1 = Queue


def _model(N, lam_delay, service, servers, lld=None):
    model = Network('ldqbd_mphc')
    delay = Delay(model, 'Delay')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    queue.setNumberOfServers(servers)
    jobclass = ClosedClass(model, 'Class1', N, delay)
    delay.setService(jobclass, Exp(lam_delay))
    queue.setService(jobclass, service)
    if lld is not None:
        queue.setLoadDependence(np.asarray(lld, dtype=float))
    model.link(model.serialRouting(delay, queue))
    return model


def _mam_vs_ctmc(model):
    mam = SolverMAM(model, method='ldqbd').getAvg()
    ctmc = SolverCTMC(model).getAvg()
    # (QLen, Util, RespT, Tput) at the queue
    return ([np.asarray(m)[QI, 0] for m in mam[:4]],
            [np.asarray(c)[QI, 0] for c in ctmc[:4]])


@pytest.mark.parametrize("servers", [1, 2, 3])
@pytest.mark.parametrize("service", [
    Erlang.fitMeanAndOrder(1.0, 2),
    Erlang.fitMeanAndOrder(1.0, 3),
    HyperExp(0.6, 2.0, 0.5),
    Coxian.fitMeanAndSCV(1.0, 3.0),
])
def test_ph_multiserver_is_exact_against_ctmc(service, servers):
    # the fix: every one of these was ~1e-2 off at c > 1 before the multiset
    got, want = _mam_vs_ctmc(_model(5, 0.8, service, servers))
    for g, w in zip(got, want):
        assert g == pytest.approx(w, abs=TOL, rel=TOL)


@pytest.mark.parametrize("servers", [1, 2, 3, 5])
def test_exponential_is_untouched(servers):
    # the exponential path does not go through the multiset builder at all, and
    # must keep matching the M/M/c boundary exactly
    got, want = _mam_vs_ctmc(_model(6, 0.5, Exp(1.0), servers))
    for g, w in zip(got, want):
        assert g == pytest.approx(w, abs=TOL, rel=TOL)


def test_load_dependence_with_ph_service():
    # load dependence shares the aggregate factor over the busy servers, so at
    # c == 1 it is the single server running sf(n) times faster.
    #
    # UTILIZATION IS COMPARED HERE TOO, and that is the point of the case. MAM
    # used to report P(busy) = 1 - pi(0) under load dependence, which reads a
    # server running alpha(n) times faster as no busier than one at its nominal
    # rate: 0.9587 against CTMC's 0.6612 on exactly this model. It now reports
    # the work-based sum_n pi(n)*sf(n)/max(c, max(alpha)), CTMC's own
    # convention, so all four metrics are exact against the exact chain.
    got, want = _mam_vs_ctmc(_model(4, 1.0, Erlang.fitMeanAndOrder(1.0, 2), 1,
                                    lld=[1.0, 1.5, 2.0, 2.5]))
    for g, w in zip(got, want):
        assert g == pytest.approx(w, abs=TOL, rel=TOL)


def test_multiset_counts_and_order():
    # nchoosek(k+p-1, p-1) configurations, and k == 1 must give the identity
    # rows in phase order -- that ordering is what makes c == 1 coincide with
    # plain phase indexing
    assert ph_multisets(3, 0).tolist() == [[0, 0, 0]]
    assert ph_multisets(3, 1).tolist() == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    assert ph_multisets(2, 2).tolist() == [[2, 0], [1, 1], [0, 2]]
    assert ph_multisets(4, 3).shape[0] == 20          # comb(6,3)
    assert ph_multisets(3, 4).shape[0] == 15          # comb(6,2)


def test_blocks_form_a_generator():
    # every level's three blocks must sum to zero rows: that is what makes the
    # block-tridiagonal matrix a generator, and it catches a miscounted
    # transition anywhere in the construction
    D0 = np.array([[-4.0, 4.0], [0.0, -4.0]])
    D1 = np.array([[0.0, 0.0], [4.0, 0.0]])
    alpha = np.array([1.0, 0.0])
    for c in (1, 2, 3):
        nlev = 6
        arr = np.array([(nlev - n) * 0.7 for n in range(nlev + 1)])
        Q0, Q1, Q2 = ldqbd_mphc(D0, D1, alpha, c, arr)
        assert len(Q0) == nlev and len(Q1) == nlev + 1 and len(Q2) == nlev
        for n in range(nlev + 1):
            rows = Q1[n].sum(axis=1)
            if n < nlev:
                rows = rows + Q0[n].sum(axis=1)
            if n >= 1:
                rows = rows + Q2[n - 1].sum(axis=1)
            assert np.allclose(rows, 0.0, atol=1e-12)
        # level sizes: comb(min(n,c)+p-1, p-1) with p = 2, so min(n,c)+1
        for n in range(nlev + 1):
            assert Q1[n].shape[0] == min(n, c) + 1


def test_load_dependent_speed_reduces_to_the_plain_chain():
    # passing sf(n) == min(n,c) must reproduce the unscaled blocks bit for bit,
    # which is what keeps the no-load-dependence path free of a stray division
    D0 = np.array([[-3.0, 3.0], [0.0, -3.0]])
    D1 = np.array([[0.0, 0.0], [3.0, 0.0]])
    alpha = np.array([1.0, 0.0])
    nlev, c = 5, 2
    arr = np.array([(nlev - n) * 0.9 for n in range(nlev + 1)])
    sf = np.array([min(n, c) for n in range(1, nlev + 1)], dtype=float)
    plain = ldqbd_mphc(D0, D1, alpha, c, arr)
    scaled = ldqbd_mphc(D0, D1, alpha, c, arr, sf)
    for a, b in zip(plain[0] + plain[1] + plain[2], scaled[0] + scaled[1] + scaled[2]):
        assert np.array_equal(a, b)
