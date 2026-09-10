"""
Regression: MAM must apply the Utilization Law at a multiserver station.

solver_mam_basic computed UN = TN*S without dividing by the number of servers,
even though it scales the service PH by 1/nservers for the queueing analysis.
For an M/M/3 with rho = 0.4 that yields an offered load of 1.2, so MAM declared a
stable station saturated: it reported Util = 1.0 AND clamped Tput to the
single-server capacity. Both throughput and utilization were wrong.

Oracle: the Utilization Law U = X*S/c is exact for an open M/M/c, and MVA, CTMC,
JMT and SSA all agree on it independently.
"""

import numpy as np
import pytest

from line_solver import (Network, Source, Queue, Sink, OpenClass, Exp,
                         SchedStrategy, SolverMAM, SolverMVA)

TOL = 1e-6


def _mmc(lam, mu, c, sched):
    m = Network('mmc')
    s = Source(m, 'S')
    q = Queue(m, 'Q', sched)
    k = Sink(m, 'K')
    cl = OpenClass(m, 'C')
    s.setArrival(cl, Exp(lam))
    q.setService(cl, Exp(mu))
    q.setNumberOfServers(c)
    m.link(Network.serialRouting(s, q, k))
    return m


@pytest.mark.parametrize('lam,mu,c,sched', [
    (1.2, 1.0, 3, SchedStrategy.PS),     # the reported case: rho = 0.4
    (1.2, 1.0, 3, SchedStrategy.FCFS),   # FCFS multiserver was wrong too
    (1.5, 1.0, 2, SchedStrategy.PS),     # rho = 0.75
    (0.6, 1.0, 1, SchedStrategy.PS),     # single server: must be unchanged
    (0.6, 1.0, 1, SchedStrategy.FCFS),
], ids=['ps-c3', 'fcfs-c3', 'ps-c2', 'ps-c1', 'fcfs-c1'])
def test_mam_multiserver_util_follows_utilization_law(lam, mu, c, sched):
    m = _mmc(lam, mu, c, sched)
    t = SolverMAM(m, verbose=False).getAvgTable()
    util = float(np.asarray(t['Util'])[1])
    tput = float(np.asarray(t['Tput'])[1])
    # U = lambda*S/c, exact for an open M/M/c
    assert util == pytest.approx(lam * (1.0 / mu) / c, abs=TOL)
    # a stable station must carry the full offered throughput
    assert tput == pytest.approx(lam, abs=TOL)


@pytest.mark.parametrize('lam,mu,c,sched', [
    (1.2, 1.0, 3, SchedStrategy.PS),
    (1.2, 1.0, 3, SchedStrategy.FCFS),
    (1.5, 1.0, 2, SchedStrategy.PS),
], ids=['ps-c3', 'fcfs-c3', 'ps-c2'])
def test_mam_agrees_with_mva_on_multiserver(lam, mu, c, sched):
    """MVA is exact for the product-form M/M/c and is an independent recursion."""
    a = SolverMAM(_mmc(lam, mu, c, sched), verbose=False).getAvgTable()
    b = SolverMVA(_mmc(lam, mu, c, sched), verbose=False).getAvgTable()
    np.testing.assert_allclose(np.asarray(a['Util'], dtype=float),
                               np.asarray(b['Util'], dtype=float), atol=TOL)
    np.testing.assert_allclose(np.asarray(a['Tput'], dtype=float),
                               np.asarray(b['Tput'], dtype=float), atol=TOL)


def test_mam_still_saturates_a_genuinely_unstable_station():
    """The fix must not defeat the instability cap: rho = 1.2 > 1 on one server."""
    t = SolverMAM(_mmc(1.2, 1.0, 1, SchedStrategy.PS), verbose=False).getAvgTable()
    assert float(np.asarray(t['Util'])[1]) == pytest.approx(1.0, abs=TOL)
