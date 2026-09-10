"""Variational inference for Markovian queueing networks (Perez-Casale, AAP 53(3), 2021).

The estimator carries no random-number stream: the expectations over the other
transitions are taken on a Halton lattice mapped through the inverse marginal
c.d.f. The fixture below must therefore reproduce the MATLAB, Java and C++
implementations digit for digit, and the golden values are the MATLAB
reference's.
"""
import numpy as np

from line_solver import (Network, Delay, Queue, Exp, ClosedClass, SchedStrategy,
                         MetricType)
from line_solver.inference import ParamEstimator, SampledMetric
from line_solver.inference.api import (VariationalSpec, VariationalOptions,
                                       infer_variational)

OBSQ = np.array([1, 2, 3, 2, 4, 3, 5, 4, 3, 4], dtype=float)
N = 10
LAM = 0.5


def _closed_loop_spec():
    """Closed two-station loop of the paper's Section 6.1, on fixed readings."""
    return VariationalSpec(
        arcs=[[1, 2, 1], [2, 1, 1]], x0=[[N], [0]], sched=[0, 1], nservers=[1, 1],
        routeprob=[1, 1], arcparam=[0, 1], arcrate=[LAM, np.nan],
        alpha0=[2.0], beta0=[1.0], obsTimes=np.arange(1, 11, dtype=float),
        obsData=np.column_stack([N - OBSQ, OBSQ]), obsRange=[[N], [N]],
        epsilon=0.1, capacity=[[N], [N]])


def _fixture_options():
    return VariationalOptions(ngrid=51, nsamples=32, ymax=60, iter_max=5, tol=0.0,
                              delta=1e-3)


def test_infer_variational_matches_matlab():
    out = infer_variational(_closed_loop_spec(), _fixture_options())

    assert abs(out.alpha[0] - 20.8738918926) < 1e-8
    assert abs(out.beta[0] - 10.4687500000) < 1e-8
    assert abs(out.rates[0] - 1.9939240017) < 1e-8
    assert out.iter == 5

    expected_bound = [-89.392514468, -164.0463074709, -163.7108480833,
                      -174.6500856726, -163.2336393829]
    assert np.allclose(out.bound, expected_bound, atol=1e-6)

    yv = np.arange(out.Y.shape[2])
    assert abs(out.Y[0][-1].dot(yv) - 22.7916783314) < 1e-8
    assert abs(out.Y[1][-1].dot(yv) - 18.8738918926) < 1e-8
    assert abs(out.qlen[-1, 0] - 6.0822135613) < 1e-8
    assert abs(out.qlen[-1, 1] - 3.9177864387) < 1e-8
    # the two stations hold the whole closed population at every epoch
    assert np.allclose(out.qlen.sum(axis=1), N, atol=1e-9)


def test_infer_variational_single_transition_is_exact():
    """With one transition the mean field is exact, so the marginal must
    reproduce the transient of the underlying birth process."""
    npop, lam, tmax = 50, 0.1, 20.0
    spec = VariationalSpec(
        arcs=[[1, 0, 1]], x0=[[npop]], sched=[0], nservers=[1], routeprob=[1],
        arcparam=[0], arcrate=[lam], alpha0=[], beta0=[], obsTimes=[tmax],
        obsData=[[np.nan]], obsRange=[[npop]], epsilon=0.2, capacity=[[npop]])
    out = infer_variational(spec, VariationalOptions(ngrid=201, nsamples=50,
                                                     iter_max=4, tmax=tmax))
    yv = np.arange(out.Y.shape[2])
    assert abs(out.Y[0][-1].dot(yv) - npop * (1 - np.exp(-lam * tmax))) < 2e-2
    assert out.tailmass < 1e-6


def test_estimator_reaches_the_api_fixture():
    """The same fixture driven through ParamEstimator, so that the translation
    from the network to the transition set is pinned to the api golden."""
    model = Network('vi')
    think = Delay(model, 'Think')
    queue = Queue(model, 'Q', SchedStrategy.FCFS)
    cl = ClosedClass(model, 'C', N, think, 0)
    think.setService(cl, Exp(LAM))
    queue.setService(cl, Exp(2.0))
    model.link(Network.serialRouting(think, queue))

    t = np.arange(1, 11, dtype=float)
    pe = ParamEstimator(model)
    pe.addSamples(SampledMetric(MetricType.QLen, t, N - OBSQ, think, cl))
    pe.addSamples(SampledMetric(MetricType.QLen, t, OBSQ, queue, cl))
    pe.options['method'] = 'vi'
    pe.options['epsilon'] = 0.1
    pe.options['prior_shape'] = 2.0  # with the model rate 2.0 this gives Gamma(2,1)
    pe.options['ngrid'] = 51
    pe.options['nsamples'] = 32
    pe.options['ymax'] = 60
    pe.options['iter_max'] = 5
    pe.options['tol'] = 0.0
    pe.options['delta'] = 1e-3
    pe.options['verbose'] = 0

    est = np.ravel(pe.estimate_at([queue]))
    post = np.atleast_2d(pe.options['posterior'])
    assert abs(post[0, 0] - 20.8738918926) < 1e-8
    assert abs(post[0, 1] - 10.4687500000) < 1e-8
    assert abs(est[0] - 10.4687500000 / 20.8738918926) < 1e-9
