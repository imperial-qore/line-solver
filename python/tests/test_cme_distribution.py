"""Tests for the CME (Concentrated Matrix Exponential) distribution.

Mirrors jar/src/test/java/jline/lang/processes/CMEDistributionTest.java and
line-test.git/test/testsAPI/test_cme_distribution.m.
"""

import math

import numpy as np
import pytest

from line_solver import (CME, ME, ClosedClass, Delay, Erlang, Exp, Network, OpenClass, Queue,
                         Sink, Source, SchedStrategy, SolverCTMC, SolverMAM, fit_me_mean_scv)
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


def test_sn_isph_marks_me_stations():
    # sn.isph marks the CME station as non-Markovian: its (D0,D1) pair has
    # negative off-diagonal entries, so mu/phi/pie carry no probabilistic
    # reading there, while the exponential source stays Markovian.
    model = Network('isph')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1')
    source.setArrival(oclass, Exp(0.5))
    queue.setService(oclass, CME(1.0, 7))
    model.link(Network.serialRouting(source, queue, sink))

    sn = model.getStruct()
    assert bool(sn.isph[0, 0]) is True
    assert bool(sn.isph[1, 0]) is False

    # An Erlang service is a genuine phase-type and stays flagged as such.
    ph_model = Network('isph-ph')
    s2 = Source(ph_model, 'Source')
    q2 = Queue(ph_model, 'Queue', SchedStrategy.FCFS)
    k2 = Sink(ph_model, 'Sink')
    c2 = OpenClass(ph_model, 'Class1')
    s2.setArrival(c2, Exp(0.5))
    q2.setService(c2, Erlang(2, 2))
    ph_model.link(Network.serialRouting(s2, q2, k2))
    assert ph_model.getStruct().isph.all()


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


def _mme1(rho, order, cutoff=30):
    model = Network('M/CME/1')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1')
    source.setArrival(oclass, Exp(rho))
    queue.setService(oclass, CME(1.0, order))
    model.link(Network.serialRouting(source, queue, sink))
    return model


@pytest.mark.parametrize("order,rho", [(3, 0.3), (7, 0.5), (11, 0.5)])
def test_ctmc_mg1_is_exact(order, rho):
    # The ME embeds in the generator exactly as a phase-type does, keeping the
    # negative off-diagonal entries of A. The stationary vector is then a signed
    # measure, but every aggregate over a phase block is exact, so the mean queue
    # length must reproduce Pollaczek-Khinchine to solver precision.
    scv = CME(1.0, order).getSCV()
    avg = SolverCTMC(_mme1(rho, order), cutoff=30, verbose=False).getAvgTable()
    pk = rho + rho ** 2 * (1.0 + scv) / (2.0 * (1.0 - rho))
    assert float(avg.QLen[1]) == pytest.approx(pk, rel=1e-8)
    assert float(avg.Util[1]) == pytest.approx(rho, rel=1e-8)
    assert float(avg.Tput[1]) == pytest.approx(rho, rel=1e-8)


def test_ctmc_refuses_per_state_probabilities():
    # Per-state probabilities and uniformization-based transients do not exist
    # for a signed stationary vector, so they are refused rather than returned.
    solver = SolverCTMC(_mme1(0.5, 3), cutoff=5, verbose=False)
    solver.getAvgTable()
    for query in ('getProbSys', 'getProbSysAggr'):
        with pytest.raises(ValueError):
            getattr(solver, query)()
    with pytest.raises(ValueError):
        solver.getTranProbSysAggr(1.0)


def test_ctmc_closed_model_matches_simulation():
    from line_solver import SolverLDES
    model = Network('closed-cme')
    delay = Delay(model, 'Delay')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    jobclass = ClosedClass(model, 'Class1', 2, delay)
    delay.setService(jobclass, Exp(1.0))
    queue.setService(jobclass, CME(0.5, 7))
    model.link(Network.serialRouting(delay, queue))

    ctmc = SolverCTMC(model, verbose=False).getAvgTable()
    ldes = SolverLDES(model, seed=23, samples=200000, verbose=False).getAvgTable()
    # Two-sided agreement with the sample path, at simulation accuracy.
    assert float(ctmc.QLen[1]) == pytest.approx(float(ldes.QLen[1]), rel=2e-2)
    assert float(ctmc.Tput[1]) == pytest.approx(float(ldes.Tput[1]), rel=2e-2)


@pytest.mark.parametrize("scv", [0.9, 0.5, 0.2, 0.05, 0.01, 1e-3, 1e-4])
def test_fit_me_mean_scv_is_exact(scv):
    # The CME-plus-exponential convolution matches both moments exactly over the
    # whole range (sY/(1+sY), 1), which is where an Erlang needs ceil(1/scv)
    # phases and a two-phase Coxian cannot go at all.
    dist = fit_me_mean_scv(2.0, scv)
    assert dist.getMean() == pytest.approx(2.0, rel=1e-10)
    assert dist.getSCV() == pytest.approx(scv, rel=1e-6)


@pytest.mark.parametrize("scv", [0.05, 0.01, 1e-3, 1e-4])
def test_fit_me_beats_erlang_order(scv):
    # The point of the construction: O(1/n^2) instead of the Erlang O(1/n).
    assert fit_me_mean_scv(1.0, scv).getNumberOfPhases() < math.ceil(1.0 / scv)


def test_fit_me_respects_a_phase_budget():
    # Under a budget the fit returns the closest achievable SCV from below rather
    # than silently truncating an Erlang: 20 phases reach 5.7e-3, where Erlang-20
    # stops at 0.05.
    dist = fit_me_mean_scv(1.0, 1e-4, max_phases=20)
    assert dist.getNumberOfPhases() <= 20
    assert dist.getMean() == pytest.approx(1.0, rel=1e-10)
    assert dist.getSCV() < 1.0 / 20


def test_fit_me_rejects_out_of_range():
    for bad_scv in (0.0, 1.0, 1.5, -0.1):
        with pytest.raises(ValueError):
            fit_me_mean_scv(1.0, bad_scv)
    with pytest.raises(ValueError):
        fit_me_mean_scv(0.0, 0.5)


def test_fit_me_is_solvable_by_ctmc():
    # The fitted ME is a legitimate service process end to end: M/ME/1 at rho=0.5
    # must again reproduce Pollaczek-Khinchine.
    service = fit_me_mean_scv(1.0, 0.05)
    model = Network('M/ME/1 fitted')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1')
    source.setArrival(oclass, Exp(0.5))
    queue.setService(oclass, service)
    model.link(Network.serialRouting(source, queue, sink))
    avg = SolverCTMC(model, cutoff=30, verbose=False).getAvgTable()
    pk = 0.5 + 0.25 * (1.0 + 0.05) / 1.0
    assert float(avg.QLen[1]) == pytest.approx(pk, rel=1e-8)


# An ME with a negative entry in alpha whose density touches zero in the interior
# (f(1.2) = -1.5e-14 against a peak of 2.79), so it admits NO phase-type
# representation of any order. Spectrum {-1, -2 +- 2i}, mean 0.5329, SCV 2.6843.
# Same instance as line-test.git/test/testsAPI/test_mam_me_warning.m.
_NONPH_ALPHA = np.array([0.61058991931158258, -0.15547146730086722, 0.54488154798928464])
_NONPH_A = np.array([[-1.0, 0.0, 0.0], [0.0, -2.0, 2.0], [0.0, -2.0, -2.0]])


def test_nonph_me_has_an_interior_density_zero():
    dist = ME(_NONPH_ALPHA, _NONPH_A)
    # The zero sits at x = 1.2 exactly, so it is resolved on a grid that lands on
    # it; a coarse sweep only gets within 1e-7 of it.
    grid = np.linspace(1.15, 1.25, 2001)
    pdf = np.array([dist.evalPDF(x) for x in grid])
    # The interior zero is what rules out a phase-type representation of ANY
    # order: a PH density is strictly positive throughout the interior of its
    # support. Round-off can make it slightly negative, hence the tolerance.
    assert abs(pdf.min()) < 1e-10
    assert max(dist.evalPDF(x) for x in (0.3, 0.6, 2.0, 3.0)) > 0.1
    assert (_NONPH_ALPHA < 0).any()


@pytest.mark.parametrize("rho,cutoff", [(0.3, 30), (0.5, 60)])
def test_ctmc_solves_a_non_phasetype_me(rho, cutoff):
    # The CTMC branch must not depend on the ME happening to be a phase-type in
    # disguise. Pollaczek-Khinchine still applies, since the service law is a
    # genuine distribution and the arrivals are Poisson.
    service = ME(_NONPH_ALPHA, _NONPH_A)
    mean, scv = service.getMean(), service.getSCV()

    model = Network('M/ME/1 non-PH')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1')
    source.setArrival(oclass, Exp(rho / mean))
    queue.setService(oclass, ME(_NONPH_ALPHA, _NONPH_A))
    model.link(Network.serialRouting(source, queue, sink))

    assert bool(model.getStruct().isph[1, 0]) is False

    avg = SolverCTMC(model, cutoff=cutoff, verbose=False).getAvgTable()
    pk = rho + rho ** 2 * (1.0 + scv) / (2.0 * (1.0 - rho))
    assert float(avg.QLen[1]) == pytest.approx(pk, rel=1e-7)
    assert float(avg.Util[1]) == pytest.approx(rho, rel=1e-6)
