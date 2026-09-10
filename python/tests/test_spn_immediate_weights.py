"""
Regression tests for GSPN immediate-transition semantics and for the phase-type
expansion of non-Markovian firing times.

Both defects these cover were invisible to the existing SPN tests because those
tests only ever exercise a single immediate mode (no conflict) and only
exponential firing times.

1. FREE-CHOICE CONFLICT. When two immediate modes are enabled by the same
   marking, GSPN semantics require the branching to follow the ratio of their
   firing weights. LINE latches a transition's servers with an ENABLE event
   emitted at GlobalConstants.Immediate, which is the same scale as an immediate
   firing, so a mode that happened to be latched first could fire before its
   competitor was latched at all. Gating the firing on the latched server count
   then made the branching follow the latching order: with weights 1 and 3 the
   split came out at 0.34375 instead of 0.25. The fix gates an immediate firing
   on the marking (afterGlobalEvent.m:315-330 and its Java/Python ports).

2. NON-MARKOVIAN FIRING TIMES. A Transition mode carrying Det/Gamma/Weibull/
   Pareto/Uniform/Lognormal is expanded to an acyclic PH by sn_nonmarkov_toph.
   Skipping that expansion leaves the single nominal phase that the node layer
   installs, i.e. the model is silently solved as if the firing time were
   exponential with the same mean.
"""

import os

import numpy as np
import pytest

import line_solver
# The worktree copy of line_solver must be the one under test; a global install
# would silently validate the wrong code.
assert os.path.realpath(__file__).rsplit('/python/', 1)[0] in os.path.realpath(
    line_solver.__file__), (
    "test must import the worktree line_solver, got %s" % line_solver.__file__)

from line_solver import (Network, Place, Transition, Source, Sink, ClosedClass,
                         OpenClass, Exp, Det, Immediate, SolverCTMC,
                         TimingStrategy)


def _free_choice(weight1, weight2):
    """One token cycling through a free-choice conflict between two immediate
    modes. P1 -> (T1 | T2) -> P2 | P3 -> (T3 | T4) -> P1.

    In steady state every cycle passes through exactly one of T1, T2, so the
    throughputs of P2 and P3 are the branch probabilities scaled by the common
    cycle rate. Their ratio is therefore the firing-weight ratio, independently
    of which immediate server happens to be latched first.
    """
    m = Network('spn_free_choice')
    P1 = Place(m, 'P1'); P2 = Place(m, 'P2'); P3 = Place(m, 'P3')
    T1 = Transition(m, 'T1'); T2 = Transition(m, 'T2')
    T3 = Transition(m, 'T3'); T4 = Transition(m, 'T4')
    jc = ClosedClass(m, 'C', 1, P1, 0)

    a = T1.addMode('m1'); T1.setDistribution(a, Immediate())
    T1.setTimingStrategy(a, TimingStrategy.IMMEDIATE)
    T1.setFiringWeights(a, weight1)
    T1.setEnablingConditions(a, jc, P1, 1); T1.setFiringOutcome(a, jc, P2, 1)

    b = T2.addMode('m2'); T2.setDistribution(b, Immediate())
    T2.setTimingStrategy(b, TimingStrategy.IMMEDIATE)
    T2.setFiringWeights(b, weight2)
    T2.setEnablingConditions(b, jc, P1, 1); T2.setFiringOutcome(b, jc, P3, 1)

    c = T3.addMode('m3'); T3.setDistribution(c, Exp(1.0))
    T3.setEnablingConditions(c, jc, P2, 1); T3.setFiringOutcome(c, jc, P1, 1)

    d = T4.addMode('m4'); T4.setDistribution(d, Exp(1.0))
    T4.setEnablingConditions(d, jc, P3, 1); T4.setFiringOutcome(d, jc, P1, 1)

    rm = m.initRoutingMatrix()
    rm.set(jc, jc, P1, T1, 1.0); rm.set(jc, jc, P1, T2, 1.0)
    rm.set(jc, jc, T1, P2, 1.0); rm.set(jc, jc, T2, P3, 1.0)
    rm.set(jc, jc, P2, T3, 1.0); rm.set(jc, jc, T3, P1, 1.0)
    rm.set(jc, jc, P3, T4, 1.0); rm.set(jc, jc, T4, P1, 1.0)
    m.link(rm)
    P1.setState(1); P2.setState(0); P3.setState(0)
    return m


@pytest.mark.parametrize('w1,w2', [(1.0, 3.0), (1.0, 1.0), (2.0, 1.0), (1.0, 9.0)])
def test_immediate_branching_follows_firing_weights(w1, w2):
    """The branch split must equal w1/(w1+w2), not the latching order."""
    tput = np.asarray(SolverCTMC(_free_choice(w1, w2)).getAvgTput()).ravel()
    # Rows follow node creation order: P1, P2, P3.
    branch1, branch2 = tput[1], tput[2]
    assert branch1 + branch2 > 0, 'net is deadlocked'
    observed = branch1 / (branch1 + branch2)
    expected = w1 / (w1 + w2)
    assert abs(observed - expected) < 1e-6, (
        'branch probability %g, expected %g (weights %g:%g)'
        % (observed, expected, w1, w2))


def test_immediate_branching_is_not_the_latching_artifact():
    """Guard the specific regression value.

    Gating on the latched server count instead of on the marking produced
    0.34375 for weights 1:3. Assert we are not back at that value, so a
    reintroduction cannot hide behind a loose tolerance.
    """
    tput = np.asarray(SolverCTMC(_free_choice(1.0, 3.0)).getAvgTput()).ravel()
    observed = tput[1] / (tput[1] + tput[2])
    assert abs(observed - 0.34375) > 1e-3, (
        'branching reverted to the latching-order artifact (%g)' % observed)


_LAMBDA = 0.5
_MU = 1.0
_RHO = _LAMBDA / _MU
_CAP = 12


def _mg1k(firing_dist):
    """Open SPN forming an M/G/1/K queue: Source -> P1 -> T1 -> Sink.

    The firing distribution of T1 is the service time, so the mean marking of P1
    is the M/G/1 queue length and depends on the service SCV, not only its mean.
    """
    m = Network('spn_mg1k')
    source = Source(m, 'source'); sink = Sink(m, 'sink')
    p1 = Place(m, 'P1'); t1 = Transition(m, 'T1')
    jc = OpenClass(m, 'jobs')
    source.setArrival(jc, Exp.fitMean(1.0 / _LAMBDA))
    p1.setClassCapacity(jc, _CAP)
    mode = t1.addMode('fire')
    t1.setDistribution(mode, firing_dist)
    t1.setEnablingConditions(mode, jc, p1, 1)
    rm = m.initRoutingMatrix()
    rm.set(jc, jc, source, p1, 1.0)
    rm.set(jc, jc, p1, t1, 1.0)
    rm.set(jc, jc, t1, sink, 1.0)
    m.link(rm)
    return m


def _solve_p1(model):
    solver = SolverCTMC(model)
    solver.options.cutoff = _CAP + 4
    qn = solver.getAvg()[0]
    pidx = model.getStationNames().index('P1')
    return float(np.sum(qn[pidx, :]))


def _pollaczek_khinchine(scv):
    """Mean number in an M/G/1 system: rho + rho^2 (1 + SCV) / (2 (1 - rho))."""
    return _RHO + _RHO ** 2 * (1.0 + scv) / (2.0 * (1.0 - _RHO))


def test_deterministic_firing_matches_md1_not_mm1():
    """A Det firing time must expand to phase-type, not collapse to Exp.

    Det(1) and Exp(1) share a mean but differ in SCV (0 vs 1), so the M/G/1
    queue length differs: 0.75 for M/D/1 against 1.0 for M/M/1. The conversion
    approximates Det by an Erlang of nonmkvorder phases (20 by default), whose
    SCV is 1/20, so the reachable target is the Pollaczek-Khinchine value at
    SCV = 0.05, i.e. 0.7625. Equality with the M/M/1 value would mean the firing
    distribution was ignored, which is what happened while the SPN branch of the
    non-Markovian-to-PH conversion was missing from the JAR.
    """
    # The PH expansion this test is about announces itself; assert it happened.
    # Under lang='java' the solve is delegated to jline.jar over JSON, so the
    # notice is cast by snNonmarkovToPh INSIDE the JAR (on its own stream) and no
    # python warning is raised; the numeric assertions below still cover the
    # expansion there, since collapsing Det to Exp would move q_det onto q_exp.
    if os.environ.get('LINE_SOLVER_LANG') == 'java':
        q_det = _solve_p1(_mg1k(Det(1.0 / _MU)))
    else:
        with pytest.warns(UserWarning, match='non-Markovian and will be converted to PH'):
            q_det = _solve_p1(_mg1k(Det(1.0 / _MU)))
    q_exp = _solve_p1(_mg1k(Exp.fitMean(1.0 / _MU)))

    assert np.isfinite(q_det), 'Det firing produced a non-finite marking'
    # Exponential firing reproduces M/M/1 (truncation at K=12, rho=0.5, is ~1e-4).
    assert q_exp == pytest.approx(_pollaczek_khinchine(1.0), rel=2e-3)
    # Det firing must land on the Erlang-20 approximation of M/D/1, well away
    # from the M/M/1 value.
    assert q_det == pytest.approx(_pollaczek_khinchine(1.0 / 20.0), rel=5e-3), (
        'Det firing gave %g, expected ~%g (M/D/1 via Erlang-20); M/M/1 is %g'
        % (q_det, _pollaczek_khinchine(1.0 / 20.0), _pollaczek_khinchine(1.0)))
    assert q_det < q_exp, 'lower service variability must not raise the queue'
