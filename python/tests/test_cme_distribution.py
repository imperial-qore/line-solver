"""Tests for the CME (Concentrated Matrix Exponential) distribution.

Mirrors jar/src/test/java/jline/lang/processes/CMEDistributionTest.java and
line-test.git/test_cme_distribution.m.
"""

import math

import numpy as np
import pytest

from line_solver import CME, ME, Exp, Network, OpenClass, Queue, Sink, Source, SchedStrategy, SolverMAM
from line_solver.distributions.markovian import _cme_representation, _cme_table_entry


def test_is_an_me():
    cme = CME(1.0, 7)
    assert isinstance(cme, ME)
    # The process type stays ME so that every solver gate accepting ME accepts it.
    assert cme.getName() == 'ME'
    assert cme.getNumberOfPhases() == 7
    assert cme.getOrder() == 7


@pytest.mark.parametrize("order", [3, 7, 21, 101])
def test_representation_is_valid(order):
    from line_solver.lib.thirdparty.butools.ph.check import CheckMERepresentation

    alpha, A, _ = _cme_representation(order)
    assert CheckMERepresentation(alpha.reshape(1, -1), A, 1e-14)
    # abs=1e-10: the entries of alpha grow with the order (c alone is 5.2e5 at
    # n = 980) and cancel down to one, so the sum carries a rounding error of
    # about 1e-16 times the largest entry, not 1e-16.
    assert alpha.sum() == pytest.approx(1.0, abs=1e-10)


@pytest.mark.parametrize("order", [3, 7, 21, 101])
def test_mean_and_scv(order):
    cme = CME(2.5, order)
    assert cme.getMean() == pytest.approx(2.5, rel=1e-10)
    # rel=1e-6: the SCV is m2 - m1^2 with m2 = 1 + scv, so at order 101 the
    # subtraction cancels four digits, and MATLAB and the JAR lose about that
    # much through their own linear solves.
    assert cme.getSCV() == pytest.approx(CME.getMinSCV(order), rel=1e-6)


@pytest.mark.parametrize("order", [7, 21, 101])
def test_scv_below_erlang_bound(order):
    # The point of a CME: it beats the phase-type bound 1/order at equal order.
    assert CME.getMinSCV(order) < 1.0 / order


def test_density_matches_tabulated_form():
    # f(x) = mu1*exp(-mu1*x)*(c + sum_k [a_k cos(k w mu1 x) + b_k sin(k w mu1 x)])
    entry = _cme_table_entry(21)
    cme = CME(1.0, 21)
    for x in (0.3, 0.8, 1.0, 1.4):
        closed = sum(entry['a'][k - 1] * math.cos(k * entry['omega'] * entry['mu1'] * x)
                     + entry['b'][k - 1] * math.sin(k * entry['omega'] * entry['mu1'] * x)
                     for k in range(1, entry['n'] + 1))
        closed = entry['mu1'] * math.exp(-entry['mu1'] * x) * (entry['c'] + closed)
        assert cme.evalPDF(x) == pytest.approx(closed, rel=1e-8, abs=1e-12)


def test_density_is_nonnegative():
    cme = CME(1.0, 21)
    for i in range(401):
        assert cme.evalPDF(i * 0.01) > -1e-10


def test_scaling_with_mean():
    # The representation scales as A/mean, so the CDF is a pure time rescaling.
    unit = CME(1.0, 11)
    scaled = CME(3.0, 11)
    for x in (0.2, 0.7, 1.3, 2.0):
        assert unit.evalCDF(x) == pytest.approx(scaled.evalCDF(3.0 * x), rel=1e-9, abs=1e-12)
    assert unit.getSCV() == pytest.approx(scaled.getSCV(), rel=1e-12)


def test_fit_mean_and_scv():
    cme = CME.fitMeanAndSCV(4.0, 1e-3)
    assert cme.getMean() == pytest.approx(4.0, rel=1e-10)
    assert cme.getSCV() <= 1e-3
    # Minimality: no lower tabulated order reaches the target.
    for order in CME.getSupportedOrders():
        if order < cme.getOrder():
            assert CME.getMinSCV(order) > 1e-3


def test_sampling():
    cme = CME(1.0, 11)
    samples = np.asarray(cme.sample(20000))
    assert (samples >= 0).all()
    mean = samples.mean()
    assert mean == pytest.approx(1.0, abs=0.02)
    assert samples.var() / mean ** 2 == pytest.approx(cme.getSCV(), abs=0.01)


def test_invalid_arguments():
    with pytest.raises(ValueError):
        CME(1.0, 4)
    with pytest.raises(ValueError):
        CME(1.0, 1)
    with pytest.raises(ValueError):
        CME(0.0, 7)
    with pytest.raises(ValueError):
        CME(-1.0, 7)
    with pytest.raises(ValueError):
        # 2003 phases is beyond the last tabulated entry (n = 1000).
        CME(1.0, 2003)
    with pytest.raises(ValueError):
        CME.fitMeanAndSCV(1.0, 1e-12)


def test_supported_orders():
    orders = CME.getSupportedOrders()
    assert len(orders) > 100
    assert orders[0] == 3
    assert all(o % 2 == 1 for o in orders)
    assert all(orders[i] > orders[i - 1] for i in range(1, len(orders)))


def test_mg1_queue():
    # M/CME/1 at rho = 0.5 with a near-deterministic service: the Pollaczek
    # -Khinchine mean queue length is rho^2*(1+scv)/(2*(1-rho)) plus rho.
    model = Network('M/CME/1')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1')
    service = CME(1.0, 11)
    source.setArrival(oclass, Exp(0.5))
    queue.setService(oclass, service)
    model.link(Network.serialRouting(source, queue, sink))

    qlen = float(SolverMAM(model).getAvgTable().QLen[1])

    rho = 0.5
    pk = rho + rho ** 2 * (1.0 + service.getSCV()) / (2.0 * (1.0 - rho))
    assert qlen == pytest.approx(pk, rel=1e-3)
