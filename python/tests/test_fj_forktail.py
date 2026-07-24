"""Regression tests for the native-python ForkTail port
(line_solver/api/fjnative/fj_tail_forktail.py) and for the
getPerctRespT(..., method='forktail') entry point.

Golden values come from the MATLAB implementation
(matlab/src/api/fj/fj_tail_forktail.m); the JAR reproduces them too.

Reference: Nguyen, Alesawi, Li, Che, Jiang, "ForkTail", ACM HPDC 2018."""

import math

import numpy as np
import pytest

from line_solver import (
    Network, Source, Sink, Fork, Join, Queue, OpenClass,
    Exp, Erlang, HyperExp, SchedStrategy, SolverMVA,
)
from line_solver.api.fjnative import fj_tail_forktail, fj_mg1_respt_moments, ge_fit


def _three_branch_fj():
    model = Network('fj3')
    source = Source(model, 'Source')
    fork = Fork(model, 'Fork')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
    q3 = Queue(model, 'Q3', SchedStrategy.FCFS)
    join = Join(model, 'Join', fork)
    sink = Sink(model, 'Sink')
    jobclass = OpenClass(model, 'C1')
    source.setArrival(jobclass, Exp(0.8))
    q1.setService(jobclass, Exp(1.0))
    q2.setService(jobclass, Erlang.fitMeanAndOrder(1.0, 2))
    q3.setService(jobclass, HyperExp.fitMeanAndSCV(1.0, 4))
    P = model.initRoutingMatrix()
    for a, b in ((source, fork), (fork, q1), (fork, q2), (fork, q3),
                 (q1, join), (q2, join), (q3, join), (join, sink)):
        P.set(jobclass, jobclass, a, b, 1.0)
    model.link(P)
    return model


def test_exponential_branch_is_exact_fit():
    ET = 2.5
    xp, alpha, beta = fj_tail_forktail(ET, ET**2, 1, 99)
    assert alpha == pytest.approx(1.0, abs=1e-10)
    assert beta == pytest.approx(ET, rel=1e-10)
    assert xp == pytest.approx(-ET*math.log(1-0.99), rel=1e-10)


def test_homogeneous_closed_form():
    ET, K, p = 1.3, 5, 0.95
    xp, _, _ = fj_tail_forktail(ET, ET**2, K, p)
    assert xp == pytest.approx(-ET*math.log(1-p**(1/K)), rel=1e-8)


def test_vector_path_matches_scalar_path():
    ET, VT, K = 0.8, 1.7, 4
    xs, _, _ = fj_tail_forktail(ET, VT, K, 99)
    xv, _, _ = fj_tail_forktail([ET]*K, [VT]*K, None, 99)
    assert xv == pytest.approx(xs, rel=1e-6)


def test_shape_tracks_variability():
    a_low, _ = ge_fit(1.0, 0.25)
    a_exp, _ = ge_fit(1.0, 1.00)
    a_high, _ = ge_fit(1.0, 4.00)
    assert a_low > a_exp > a_high


def test_random_fanout_mixture():
    ET, VT, p = 1.4, 2.6, 99
    x_degenerate, _, _ = fj_tail_forktail(ET, VT, [3, 5], p, [1.0, 0.0])
    x3, _, _ = fj_tail_forktail(ET, VT, 3, p)
    x9, _, _ = fj_tail_forktail(ET, VT, 9, p)
    x_mix, _, _ = fj_tail_forktail(ET, VT, [3, 9], p, [0.5, 0.5])
    assert x_degenerate == pytest.approx(x3, rel=1e-6)
    assert x3 < x_mix < x9
    with pytest.raises(ValueError):
        fj_tail_forktail(ET, VT, [3, 9], p, [0.5, 0.6])


def test_mg1_moments_match_mm1():
    mu, lam = 1.0, 0.7
    ET, VT = fj_mg1_respt_moments(lam, 1/mu, 2/mu**2, 6/mu**3)
    assert ET == pytest.approx(1/(mu-lam), rel=1e-10)
    assert VT == pytest.approx(ET**2, rel=1e-10)
    with pytest.raises(ValueError):
        fj_mg1_respt_moments(1.2, 1.0, 2.0, 6.0)


def test_hyperexp_third_moment_is_not_symmetric():
    """HyperExp.getSkew used to return the base-class 0, which silently
    understated E[S^3] for every consumer that reconstructs it."""
    h = HyperExp.fitMeanAndSCV(1.0, 4)
    ES, VS = h.getMean(), h.getVar()
    ES2 = VS + ES**2
    ES3 = h.getSkew()*VS**1.5 + 3*ES*ES2 - 2*ES**3
    exact = 6*float(np.sum(h._probs*(1.0/h._rates)**3))
    assert ES3 == pytest.approx(exact, rel=1e-9)
    assert ES3 > 100


def test_getperctrespt_forktail_matches_matlab():
    """MATLAB and the JAR both give 45.838545 and 80.521692 on this model."""
    model = _three_branch_fj()
    PercRT, PercTable = SolverMVA(model).getPerctRespT([95, 99], None, 'forktail')
    assert PercRT[0]['class'] == 'C1'
    assert PercRT[0]['method'] == 'forktail'
    assert PercRT[0]['values'][0] == pytest.approx(45.838545, abs=1e-4)
    assert PercRT[0]['values'][1] == pytest.approx(80.521692, abs=1e-4)
    assert len(PercTable) == 2


def test_getperctrespt_forktail_rejects_model_without_fork():
    model = Network('noFork')
    source = Source(model, 'Source')
    queue = Queue(model, 'Q1', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    jobclass = OpenClass(model, 'C1')
    source.setArrival(jobclass, Exp(0.5))
    queue.setService(jobclass, Exp(1.0))
    P = model.initRoutingMatrix()
    P.set(jobclass, Network.serialRouting(source, queue, sink))
    model.link(P)
    with pytest.raises(ValueError):
        SolverMVA(model).getPerctRespT([99], None, 'forktail')
