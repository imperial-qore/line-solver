"""RCAT (inap/inapplus/inapinf) with phase-type processes.

Each component is a QBD over (queue length, phase) rather than the scalar
birth-death chain the analyzer built before, the phase being the pair (arrival
phase, service phase) in the Kronecker order of qbd_mapmap1.

WHAT MAKES THESE ORACLES. An isolated M/PH/1 has no synchronizing action, so
whatever the reversed-rate iterate does its marginal is the M/G/1 one and its
mean is the Pollaczek-Khinchine value rho + rho^2 (1+scv) / (2 (1-rho)), written
out here from the arrival rate and the SCV alone. The same holds for the first
station of a tandem. So these assert the phase construction itself, not the
fixed point -- and they would ALL have failed before it existed, because the
analyzer then returned the M/M/1 answer for every SCV (measured: 1.000000 at
rho=0.5 and 4.000000 at rho=0.8, whatever the service law).

The PH/M/1 cases are checked against SolverCTMC instead, since a non-Poisson
arrival stream has no P-K formula.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np
import pytest

from line_solver import (Network, Source, Queue, Sink, OpenClass, SchedStrategy,
                         Exp, Erlang, HyperExp, Coxian, Det, SolverAG, SolverCTMC)

METHODS = ('inap', 'inapplus', 'inapinf')


def pk(rho, scv):
    """M/G/1 mean number in system (Pollaczek-Khinchine)."""
    return rho + rho ** 2 * (1.0 + scv) / (2.0 * (1.0 - rho))


def single(arr, svc):
    m = Network('mg1')
    s = Source(m, 'S')
    q = Queue(m, 'Q', SchedStrategy.FCFS)
    k = Sink(m, 'K')
    oc = OpenClass(m, 'C')
    s.setArrival(oc, arr)
    q.setService(oc, svc)
    m.link(Network.serialRouting(s, q, k))
    return m


def tandem(arr, s1, s2):
    m = Network('t')
    s = Source(m, 'S')
    q1 = Queue(m, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Q2', SchedStrategy.FCFS)
    k = Sink(m, 'K')
    oc = OpenClass(m, 'C')
    s.setArrival(oc, arr)
    q1.setService(oc, s1)
    q2.setService(oc, s2)
    m.link(Network.serialRouting(s, q1, q2, k))
    return m


@pytest.mark.parametrize('rho', [0.5, 0.8])
@pytest.mark.parametrize('scv', [0.5, 0.25])
@pytest.mark.parametrize('method', ['inap', 'inapinf'])
def test_mph1_marginal_is_the_exact_pk_mean(rho, scv, method):
    model = single(Exp(rho), Erlang.fitMeanAndSCV(1.0, scv))
    sv = SolverAG(model, method)
    q = np.asarray(sv.getAvgQLen()).ravel()[1]
    u = np.asarray(sv.getAvgUtil()).ravel()[1]
    t = np.asarray(sv.getAvgTput()).ravel()[1]
    assert q == pytest.approx(pk(rho, scv), abs=1e-6)
    assert u == pytest.approx(rho, abs=1e-6)
    assert t == pytest.approx(rho, abs=1e-6)


@pytest.mark.parametrize('method', METHODS)
def test_exponential_is_unchanged_by_the_phase_construction(method):
    """A component with one phase per level IS the old birth-death chain."""
    rho = 0.5
    q_erl1 = np.asarray(SolverAG(single(Exp(rho), Erlang(1.0, 1)), method).getAvgQLen()).ravel()[1]
    q_exp = np.asarray(SolverAG(single(Exp(rho), Exp(1.0)), method).getAvgQLen()).ravel()[1]
    assert q_erl1 == pytest.approx(q_exp, abs=1e-9)


@pytest.mark.parametrize('method', ['inap', 'inapinf'])
@pytest.mark.parametrize('scv', [0.5, 4.0])
def test_ph_arrivals_match_ctmc(method, scv):
    """A non-Poisson Source needs the arrival-phase dimension of the component."""
    arr = Erlang.fitMeanAndSCV(2.0, scv) if scv <= 1 else HyperExp.fitMeanAndSCV(2.0, scv)
    model = single(arr, Exp(1.0))
    ref = np.asarray(SolverCTMC(model, cutoff=60).getAvgQLen()).ravel()[1]
    got = np.asarray(SolverAG(model, method).getAvgQLen()).ravel()[1]
    assert got == pytest.approx(ref, rel=1e-6)


@pytest.mark.parametrize('method', METHODS)
def test_tandem_with_ph_service_conserves_flow(method):
    """The upstream station is an isolated M/PH/1, so the reversed rate must be
    the arrival rate and both stations must carry it."""
    lam = 0.5
    model = tandem(Exp(lam), Erlang.fitMeanAndSCV(1.0, 0.5), Exp(2.0))
    sv = SolverAG(model, method)
    q = np.asarray(sv.getAvgQLen()).ravel()
    t = np.asarray(sv.getAvgTput()).ravel()
    assert q[1] == pytest.approx(pk(lam, 0.5), abs=1e-6)
    assert t[1] == pytest.approx(lam, abs=1e-6)
    assert t[2] == pytest.approx(lam, abs=1e-6)


def test_inapinf_removes_the_truncation_on_a_phase_component():
    """'inapinf' solves the open component on its infinite state space, so on a
    heavy-tailed service law it lands on the exact P-K mean where the 100-level
    truncation of 'inap' cannot."""
    lam = 0.8
    model = tandem(Exp(lam), HyperExp.fitMeanAndSCV(1.0, 4.0), Exp(4.0))
    exact = pk(lam, 4.0)
    q_inap = np.asarray(SolverAG(model, 'inap').getAvgQLen()).ravel()[1]
    q_inf = np.asarray(SolverAG(model, 'inapinf').getAvgQLen()).ravel()[1]
    assert q_inf == pytest.approx(exact, abs=1e-6)
    assert abs(q_inap - exact) > abs(q_inf - exact)


@pytest.mark.parametrize('method', METHODS)
def test_gate_accepts_phase_type_and_refuses_what_has_no_generator(method):
    from line_solver.solvers.solver_ag.algorithms.ag_inap import rcat_supports_processes
    for svc in (Exp(1.0), Erlang(3.0, 3), HyperExp(0.5, 3.0, 10.0),
                Coxian.fitMeanAndSCV(1.0, 2.0), Det(1.0)):
        ok, reason = rcat_supports_processes(single(Exp(0.5), svc).getStruct(), method)
        assert ok, "RCAT %s must accept %s: %s" % (method, type(svc).__name__, reason)

    # build_rcat never reads sn.nservers, so a multiserver station would be
    # driven at rho = lambda/mu instead of lambda/(c*mu).
    m = single(Exp(1.2), Exp(1.0))
    m.getNodeByName('Q').setNumberOfServers(2)
    ok, reason = rcat_supports_processes(m.getStruct(), method)
    assert not ok and 'single-server' in reason
