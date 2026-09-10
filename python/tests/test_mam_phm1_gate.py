"""The PH/M/c fast-path gate of the dec.source decomposition.

`solver_mam_basic` may answer a station with an exact closed form instead of the
decomposition, and the PH/M/c arm (`qsys_phmc`) takes ANOTHER station's process
as the arrival stream. Which station that is decides whether the answer is the
model's: the reference (`solver_mam_basic.m:404`) admits exactly two, a station
whose law is genuinely non-exponential AND renewal, or the Source of a
two-station M/M/c. This port used to scan for neither, keeping whichever station
the loop reached last, so in a tandem it read a DOWNSTREAM SERVICE RATE as the
arrival rate of the queue before it -- Queue2 of the M/M/1 tandem below came back
at 2.0 against the exact 0.5, with `qsys_phmc failed: Load rho=1.5` warned for
Queue1 (rho = mu2/mu1, the ratio of two service rates) and swallowed.

The expected values here are closed forms, not solver output:
  * the tandem is Burke's theorem -- both queues are M/M/1 in isolation;
  * E2/M/1 is the GI/M/1 sigma-root L = rho/(1-sigma) with sigma the root of
    sigma = (2/(2 + mu(1-sigma)))^2;
  * M/M/2 is Erlang-C.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import warnings

import pytest

from line_solver import (Network, Source, Queue, Sink, OpenClass, SchedStrategy,
                         Exp, Erlang, SolverMAM)


def _pick(table, station, column):
    data = getattr(table, 'data', table)
    row = data[(data['Station'] == station) & (data['JobClass'] == 'C')]
    return float(row[column].iloc[0])


def _tandem(lam, mu1, mu2):
    model = Network('Tandem')
    source = Source(model, 'Source')
    queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    cls = OpenClass(model, 'C')
    source.setArrival(cls, Exp(lam))
    queue1.setService(cls, Exp(mu1))
    queue2.setService(cls, Exp(mu2))
    P = model.initRoutingMatrix()
    P.set(cls, cls, source, queue1, 1.0)
    P.set(cls, cls, queue1, queue2, 1.0)
    P.set(cls, cls, queue2, sink, 1.0)
    model.link(P)
    return model


def _single(arrival, mu, servers):
    model = Network('Single')
    source = Source(model, 'Source')
    queue = Queue(model, 'Q', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    cls = OpenClass(model, 'C')
    source.setArrival(cls, arrival)
    queue.setService(cls, Exp(mu))
    queue.setNumberOfServers(servers)
    P = model.initRoutingMatrix()
    P.set(cls, cls, source, queue, 1.0)
    P.set(cls, cls, queue, sink, 1.0)
    model.link(P)
    return model


@pytest.mark.parametrize('lam,mu1,mu2', [(1.0, 2.0, 3.0), (0.5, 1.0, 4.0)])
def test_tandem_is_two_isolated_mm1_queues(lam, mu1, mu2):
    """Burke: the departure stream of an M/M/1 in equilibrium is Poisson at the
    arrival rate, so each queue is M/M/1 with rho = lam/mu and QLen
    rho/(1-rho). The PH/M/c arm must not fire on either station -- neither
    carries a non-exponential law, and the model has three stations."""
    table = SolverMAM(_tandem(lam, mu1, mu2), 'default').getAvgTable()
    for name, mu in (('Queue1', mu1), ('Queue2', mu2)):
        rho = lam / mu
        assert _pick(table, name, 'QLen') == pytest.approx(rho / (1.0 - rho), abs=1e-6)
        assert _pick(table, name, 'Util') == pytest.approx(rho, abs=1e-6)
        assert _pick(table, name, 'Tput') == pytest.approx(lam, abs=1e-6)


def test_tandem_raises_no_spurious_stability_warning():
    """The old scan compared two SERVICE rates, so it reported rho = mu2/mu1 =
    1.5 on a network whose true load is 0.5."""
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter('always')
        SolverMAM(_tandem(1.0, 2.0, 3.0), 'default').getAvgTable()
    assert not [w for w in record if 'qsys_phmc' in str(w.message)]


def test_erlang_arrivals_still_take_the_exact_gim1_root():
    """E2/M/1, the case the PH/M/c arm exists for: a renewal non-exponential
    arrival law at the one other station. Mean interarrival 1, mu = 3."""
    table = SolverMAM(_single(Erlang.fitMeanAndOrder(1.0, 2), 3.0, 1), 'default').getAvgTable()
    # sigma = (2/(2+3(1-sigma)))^2 by fixed point, then L = rho/(1-sigma).
    sigma = 0.2
    for _ in range(200):
        sigma = (2.0 / (2.0 + 3.0 * (1.0 - sigma))) ** 2
    assert _pick(table, 'Q', 'QLen') == pytest.approx((1.0 / 3.0) / (1.0 - sigma), abs=1e-5)


def test_poisson_multiserver_still_takes_erlang_c():
    """M/M/2 with lam = 4, mu = 3: the two-station Source branch of the gate.
    a = 4/3, rho = 2/3, P0 = 1/5, Lq = P0 a^2 rho / (2 (1-rho)^2) = 16/15,
    L = Lq + a = 2.4. The generic single-fast-server surrogate is inexact here,
    which is why the branch exists."""
    table = SolverMAM(_single(Exp(4.0), 3.0, 2), 'default').getAvgTable()
    assert _pick(table, 'Q', 'QLen') == pytest.approx(2.4, abs=1e-6)
