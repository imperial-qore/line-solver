"""SolverNC's mixed load-dependent route must read the whole rate lattice.

`solver_ncld` sends an open or mixed limited-load-dependent model to the
Bruell-Balbo-Afshari effective-capacity MVA (`pfqn_mvaldmx`, method `ncldmx`),
which is exact. Two defects made it answer something else:

- the rate row handed to it was cut at the CLOSED population. `pfqn_ldmx_ec`
  infers the limited-load-dependence level b_i as the first column equal to the
  last one and treats every rate past it as saturated, so a c-server station was
  read as saturated at min(n,c) with n<c whenever c exceeded that population --
  and with no closed class at all the row collapsed to mu(1), i.e. one server.
- `pfqn_mvaldmx` itself answered a purely open model with the single-server law
  lam*D/(1-rho), dropping mu(n) entirely and returning zeros when rho>=1.

Both showed as an open class that a multiserver could not have queued that way.
The oracles are the exact M/M/c law and SolverCTMC.
"""
import math

import numpy as np
import pytest

from line_solver import (Network, Source, Sink, Queue, Delay, OpenClass, ClosedClass,
                         Exp, SchedStrategy, SolverCTMC, SolverNC)

SERVERS = 3
LLD_WIDTH = 40


def _mmc_qlen(lam, mu, c):
    """Exact mean number in an M/M/c queue."""
    a = lam / mu
    rho = a / c
    s = sum(a ** k / math.factorial(k) for k in range(c))
    tail = a ** c / (math.factorial(c) * (1 - rho))
    p0 = 1.0 / (s + tail)
    lq = p0 * a ** c * rho / (math.factorial(c) * (1 - rho) ** 2)
    return lq + a


def _model(lam, closed_jobs=0, servers=SERVERS):
    """Source -> LD queue -> Sink, plus an optional closed chain through a delay.

    The multiserver is expressed as limited load dependence mu(n)=min(n,c),
    which is the only form the load-dependent solver admits.
    """
    net = Network('mixed_ld')
    src = Source(net, 'Src')
    q = Queue(net, 'Q1', SchedStrategy.FCFS)
    snk = Sink(net, 'Snk')
    oc = OpenClass(net, 'oc')
    src.setArrival(oc, Exp(lam))
    q.setService(oc, Exp(1.0))
    if closed_jobs > 0:
        d = Delay(net, 'D')
        cc = ClosedClass(net, 'cc', closed_jobs, d)
        d.setService(cc, Exp(1.0))
        q.setService(cc, Exp(1.0))
    P = net.initRoutingMatrix()
    P.set(oc, oc, src, q, 1.0)
    P.set(oc, oc, q, snk, 1.0)
    if closed_jobs > 0:
        P.set(cc, cc, d, q, 1.0)
        P.set(cc, cc, q, d, 1.0)
    net.link(P)
    q.setLoadDependence(np.minimum(np.arange(1, LLD_WIDTH + 1), servers).astype(float))
    return net


def _qlen(table, station, jobclass):
    for i in range(len(table['Station'])):
        if table['Station'][i] == station and table['JobClass'][i] == jobclass:
            return float(table['QLen'][i])
    raise AssertionError('no %s/%s row' % (station, jobclass))


@pytest.mark.parametrize('method', ['default', 'exact'])
def test_purely_open_multiserver_is_the_mmc_queue(method):
    """No closed class: the rate row must still describe c servers, not one."""
    lam = 1.5
    got = _qlen(SolverNC(_model(lam), method=method).getAvgTable(), 'Q1', 'oc')
    assert got == pytest.approx(_mmc_qlen(lam, 1.0, SERVERS), rel=1e-9)


def test_mixed_model_with_c_above_the_closed_population():
    """c=3 with one closed job: the row used to be cut to mu(1)."""
    lam = 0.4
    nc = SolverNC(_model(lam, closed_jobs=1), method='exact').getAvgTable()
    ctmc = SolverCTMC(_model(lam, closed_jobs=1), cutoff=16).getAvgTable()
    for jobclass in ('oc', 'cc'):
        assert _qlen(nc, 'Q1', jobclass) == pytest.approx(
            _qlen(ctmc, 'Q1', jobclass), rel=1e-3)


def test_pfqn_mvaldmx_reads_mu_on_a_purely_open_network():
    """The API entry point itself, with no closed class to drive the recursion."""
    from line_solver.api.pfqn import pfqn_mvaldmx
    lam = 1.5
    mu = np.minimum(np.arange(1, LLD_WIDTH + 1), SERVERS).astype(float).reshape(1, -1)
    XN, QN = pfqn_mvaldmx(np.array([lam]), np.array([[1.0]]), np.array([np.inf]),
                          np.array([0.0]), mu, np.ones(1))[:2]
    assert XN[0] == pytest.approx(lam)
    assert QN[0, 0] == pytest.approx(_mmc_qlen(lam, 1.0, SERVERS), rel=1e-9)


def test_saturation_level_is_the_start_of_the_trailing_constant_run():
    from line_solver.api.solvers.nc.handler import _lld_saturation_level
    assert _lld_saturation_level(np.array([1.0, 2.0, 3.0, 3.0, 3.0])) == 3
    assert _lld_saturation_level(np.array([1.0, 1.0, 1.0])) == 1
    assert _lld_saturation_level(np.array([1.0, 1.6, 2.0, 2.2, 2.3])) == 5
    assert _lld_saturation_level(np.array([])) == 1


def test_ncldmx_is_classified_as_an_exact_method():
    """The route evaluates the product form itself, so the banner must say so."""
    from line_solver.solvers.base import method_type
    assert method_type('SolverNC', 'ncldmx') == 'exact, deterministic'
    assert method_type('SolverNC', 'default/ncldmx') == 'exact, deterministic'


def test_ncldmx_carries_its_citation():
    from line_solver.api.io.citations import citations_for
    entries = citations_for(['ncldmx'])
    assert len(entries) == 1
    assert entries[0]['key'] == 'AfBaBr84'
    assert 'Bruell' in entries[0]['ref']
