"""Solver-level validation of Rational Arrival Process (RAP) components.

Covers RAP arrival with exponential service, exponential arrival with RAP
service, and RAP arrival with RAP service, solved through SolverMAM. Mirrors
jar/src/test/java/jline/solvers/mam/RapQueueingModelTest.java and
line-test.git/test/testsAdvFeatures/mam/test_rap_queueing_models.m; the
reference values are the MATLAB ones, MATLAB being the ground truth of this
codebase.

Two oracles are used, and neither is a tolerance chosen to fit.

1. The MAP control is an exact identity. A RAP whose (H0, H1) happen to be
   nonnegative IS a Markovian arrival process, so SolverMAM on the RAP object
   must return what SolverMAM returns on the equivalent MAP object. Both reach
   the exact MAP/MAP/1 queue, so what this pins down is the RAP side of the sn
   marshalling - that a RAP object delivers its (H0, H1) to the solver
   unaltered and is recognised as Markovian when it is - rather than an
   agreement between two algorithms. Being an identity, the tolerance is
   numerical only.

2. For a genuinely non-Markovian RAP (H0 carries a negative off-diagonal, so
   the pair is not a MAP) there is no closed form, and SolverLDES is the
   stochastic reference: its RAP sampler is an independent conditional-vector
   recursion, so a comparison against it is a real cross-check and not a
   restatement of the matrix-analytic result. Where a value below is pinned, it
   is one all three codebases agree on AND that sits close to that reference;
   where SolverLDES contradicts the analytic answer the model is covered by
   inequalities only and the discrepancy is recorded in the comment, because a
   golden regenerated from a wrong answer hides the defect instead of holding
   the line on it.

The non-Markovian RAP used here is
    H0 = [[-9.9, -0.2], [0.1, -1.0]] * 30/119,  H1 = [[9.7, 0.4], [0, 0.9]] * 30/119
with mean 1, SCV 4.5448028674 and lag-1 autocorrelation 0.3432018450. Its
dominant left eigenvector is strictly positive ([1, 88.98]), which keeps the
conditional density of one sign over the reachable set and so makes the
representation a genuine point process, not merely an algebraically admissible
(H0, H1) pair. That distinction matters: the nearby family obtained by driving
the same off-diagonal negative from an MMPP2 has a sign-mixed dominant
eigenvector, is not a point process at all, and no sampler can reproduce its
nominal moments.
"""
import numpy as np
import pytest

from line_solver import (MAP, RAP, Exp, Network, OpenClass, Queue,
                         SchedStrategy, Sink, SolverMAM, Source)

# Identity between a RAP object and the equivalent MAP object: numerical only.
IDENTITY_TOL = 1e-9

# Agreement with the MATLAB reference. Both codebases run the same algorithm on
# the same matrices; the margin absorbs the QBD level-truncation cutoff, the
# only place the ports may legitimately differ.
MATLAB_TOL = 1e-6

NM_SCALE = 30.0 / 119.0

# MMPP2 scaled to unit rate. Nonnegative, so this pair is a MAP.
MAP_H0 = np.array([[-10.1, 0.1], [0.1, -1.1]]) / 5.5
MAP_H1 = np.array([[10.0, 0.0], [0.0, 1.0]]) / 5.5

# Genuinely non-Markovian: H0[0, 1] < 0.
NM_H0 = np.array([[-9.9, -0.2], [0.1, -1.0]]) * NM_SCALE
NM_H1 = np.array([[9.7, 0.4], [0.0, 0.9]]) * NM_SCALE

# Partner of (NM_H0, NM_H1) with the SAME H0 and H1 = (-H0 e) pie, where pie is
# the arrival-embedded equilibrium vector of the correlated process. The
# embedded matrix (-H0)^-1 H1 is then rank one, so the stream is a renewal
# process with zero autocorrelation at every lag while the marginal - hence the
# mean and the SCV - is unchanged.
NM_H1_UNCORR = np.array([[8.232773109243697, 1.867226890756303],
                         [0.7336134453781513, 0.16638655462184874]]) * NM_SCALE


def _build(arrival, service):
    m = Network('rapqn')
    src = Source(m, 'Src')
    q = Queue(m, 'Q1', SchedStrategy.FCFS)
    sink = Sink(m, 'Snk')
    cls = OpenClass(m, 'C1')
    src.setArrival(cls, arrival)
    q.setService(cls, service)
    m.link(Network.serialRouting(src, q, sink))
    return m


def _qlen(model):
    return float(np.asarray(SolverMAM(model).getAvgTable().QLen)[1])


class TestMapControlIdentity:
    """A RAP with nonnegative matrices is a MAP; the solver must not care which
    object carried it."""

    def test_rap_arrival_matches_map_arrival(self):
        q_map = _qlen(_build(MAP(MAP_H0, MAP_H1), Exp(2.0)))
        q_rap = _qlen(_build(RAP(MAP_H0, MAP_H1), Exp(2.0)))
        assert q_rap == pytest.approx(q_map, abs=IDENTITY_TOL)
        # Guard against the identity holding on a degenerate value: a correlated
        # arrival must push the queue well past the M/M/1 value of 1.0 at the
        # same utilization.
        assert q_map > 3.0

    def test_rap_service_matches_map_service(self):
        q_map = _qlen(_build(Exp(1.0), MAP(MAP_H0 * 2, MAP_H1 * 2)))
        q_rap = _qlen(_build(Exp(1.0), RAP(MAP_H0 * 2, MAP_H1 * 2)))
        assert q_rap == pytest.approx(q_map, abs=IDENTITY_TOL)
        assert q_map > 8.0

    def test_rap_rap_matches_map_map(self):
        q_map = _qlen(_build(MAP(MAP_H0, MAP_H1), MAP(MAP_H0 * 2, MAP_H1 * 2)))
        q_rap = _qlen(_build(RAP(MAP_H0, MAP_H1), RAP(MAP_H0 * 2, MAP_H1 * 2)))
        assert q_rap == pytest.approx(q_map, abs=IDENTITY_TOL)
        assert q_map > 16.0


class TestNonMarkovianRap:
    """Genuinely non-Markovian RAPs against the MATLAB reference."""

    def test_rap_service_matches_matlab(self):
        # MATLAB SolverMAM = 6.36312174319464, equal to a direct
        # qbd_raprap1({-1,1},{H0*2,H1*2}) to 12 digits, so the model really does
        # reach the RAP/RAP/1 QBD rather than a phase-type surrogate.
        q = _qlen(_build(Exp(1.0), RAP(NM_H0 * 2, NM_H1 * 2)))
        assert q == pytest.approx(6.36312174319464, abs=MATLAB_TOL)

    # The RAP-arrival-with-RAP-service case is deliberately NOT pinned to a
    # value. Python returns 11.58346574053015 for it, MATLAB and the JAR return
    # 10.80992052415996 (they clip the negative entries of the arrival
    # representation on the way in, which replaces the process), and SolverLDES
    # puts the truth at 16.13581322 +/- 0.10309226 over 20 replications of 2e6
    # samples. All three analytic answers are far outside that interval, so a
    # golden here would enshrine a defect. The autocorrelation test below is
    # what covers this model until the dispatch and the accuracy of
    # qbd_raprap1 for genuinely non-Markovian blocks are resolved.

    def test_rap_service_beats_the_exponential_baseline(self):
        # Sanity floor: a correlated non-Markovian service at rho = 0.5 must be
        # far above the M/M/1 value of 1.0. Without this an implementation that
        # collapsed the RAP to its rate would still satisfy the equality above
        # if the reference were ever regenerated from it.
        assert _qlen(_build(Exp(1.0), RAP(NM_H0 * 2, NM_H1 * 2))) > 5.0


class TestAutocorrelationPropagates:
    """The lag-1 autocorrelation must survive the sn marshalling."""

    def test_matched_marginals_differ_only_in_autocorrelation(self):
        corr = RAP(NM_H0, NM_H1)
        uncorr = RAP(NM_H0, NM_H1_UNCORR)
        assert uncorr.getMean() == pytest.approx(corr.getMean(), abs=1e-12)
        assert uncorr.getSCV() == pytest.approx(corr.getSCV(), abs=1e-12)
        assert corr.getACF(1)[0] == pytest.approx(0.34320184504427, abs=1e-12)
        assert uncorr.getACF(1)[0] == pytest.approx(0.0, abs=1e-12)

    def test_autocorrelation_changes_the_queue_length(self):
        corr = RAP(NM_H0, NM_H1)
        uncorr = RAP(NM_H0, NM_H1_UNCORR)
        service = RAP(NM_H0 * 2, NM_H1 * 2)
        q_corr = _qlen(_build(corr, service))
        q_uncorr = _qlen(_build(uncorr, service))
        # SolverMAM gives a ratio of 1.4873 here (11.5835 against 7.7885) and
        # SolverLDES puts it at 2.0013 (16.13581322 +/- 0.10309226 against
        # 8.06276774 +/- 0.05791355), so the solver understates the effect but
        # does carry it. The threshold is 1.25: far enough above 1 that it
        # cannot be met by numerical drift, which is at the 1e-6 level here, and
        # far enough below the observed 1.49 to be stable, while still failing
        # outright if the correlation were dropped in the marshalling, in which
        # case the two models would return equal values. The same threshold is
        # used by the JAR and MATLAB tests, where the ratio is 1.3905 because
        # those two clip the arrival representation.
        assert q_corr / q_uncorr > 1.25
