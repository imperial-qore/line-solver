"""
SolverMAM must not answer a correlated arrival with a closed form that reads
only the arrival's renewal marginal.

qsys_phmc solves PH/M/c from the pair (pie, D0), i.e. from the interarrival
distribution alone. It used to be reached whenever the source process was merely
non-exponential, which let a correlated MAP, RAP or ME arrival onto it and
discarded the autocorrelation: for the RAP below the answer came out at
1.98523360058856, close to what Kingman gives from the arrival SCV alone, against
a true value near 6.844. The gate is now RENEWAL, so a correlated arrival falls
through to MMAPPH1FCFS, which is handed the arrival (D0, D1) and therefore
carries its correlation. A renewal arrival, phase-type or matrix-exponential,
still takes the fast path.

Mirrors jar/src/test/java/jline/solvers/mam/MamRapArrivalTest.java and
line-test.git/test_mam_rap_arrival.m.
"""

import numpy as np

from line_solver import (Network, Source, Queue, Sink, OpenClass, Exp, Erlang,
                         HyperExp, ME, RAP, SchedStrategy, SolverMAM)

TOL = 1e-8

# Non-phase-type ME: alpha carries a negative entry and the density has an
# interior zero, so no phase-type representation of any order exists. It is
# nonetheless a RENEWAL process.
_W = 2 * np.pi
NON_PH_ALPHA = np.array([0.984694494294579, -0.040430911430916,
                         0.0557364171363366])
NON_PH_A = np.array([[-0.5, 0, 0],
                     [0, -1, _W],
                     [0, -_W, -1]])


def _correlated_rap():
    """
    Rate-1 RAP with a negative off-diagonal entry in H0, so it is a genuine
    rational arrival process and not a MAP. Mean 1, SCV 4.5448028674, lag-1
    autocorrelation 0.3432018450.
    """
    s = 30 / 119
    H0 = np.array([[-9.9, -0.2], [0.1, -1.0]]) * s
    H1 = np.array([[9.7, 0.4], [0.0, 0.9]]) * s
    return RAP(H0, H1)


def _queue_length(arrival, service, nservers=1):
    model = Network('mam_rap_arrival')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1')
    source.setArrival(oclass, arrival)
    queue.setService(oclass, service)
    queue.setNumberOfServers(nservers)
    model.link(Network.serialRouting(source, queue, sink))
    return SolverMAM(model).getAvgTable().QLen[1]


def test_correlated_rap_arrival_carries_autocorrelation():
    # 6.84449125861694 is MMAPPH1FCFS applied to the true (H0, H1); an
    # independent SolverLDES run of 60 replications of 4e6 samples brackets it at
    # [6.83540045, 6.86599427], so the analytic value is pinned against
    # simulation and not merely against itself.
    qlen = _queue_length(_correlated_rap(), Exp(2))
    assert abs(qlen - 6.84449125861694) < TOL
    assert 6.83540045 < qlen < 6.86599427


def test_moment_identical_arrivals_with_different_acf_differ():
    # The discrimination that makes this class of bug visible: two arrival
    # processes with the SAME mean and SAME SCV but different lag-1
    # autocorrelation must give materially different queue lengths. A solver
    # that reads only the marginal cannot tell them apart.
    correlated = _correlated_rap()
    mean = correlated.getMean()
    scv = correlated.getSCV()
    renewal = HyperExp.fitMeanAndSCV(mean, scv)
    assert abs(renewal.getMean() - mean) < 1e-9
    assert abs(renewal.getSCV() - scv) < 1e-9

    q_correlated = _queue_length(correlated, Exp(2))
    q_renewal = _queue_length(renewal, Exp(2))
    assert q_correlated > 3.0 * q_renewal, (
        f"moment-identical arrivals differing only in autocorrelation gave "
        f"{q_correlated} and {q_renewal}; the correlation is not reaching the "
        f"queue-length computation")


def test_renewal_me_arrival_unchanged_on_fast_path():
    # A matrix-exponential arrival is RENEWAL and must keep the fast path.
    assert abs(_queue_length(ME(NON_PH_ALPHA, NON_PH_A), Exp(1 / 0.977420))
               - 1.0370455525) < 1e-9


def test_renewal_ph_arrivals_unchanged_on_fast_path():
    # PH/M/2 at 1.1867726850 is inside the SolverLDES interval
    # [1.18614125, 1.18750740], so the fast path retained here is the accurate
    # one.
    assert abs(_queue_length(Erlang.fitMeanAndOrder(1.0, 2), Exp(2)) - 0.8090169944) < 1e-9
    assert abs(_queue_length(Erlang.fitMeanAndOrder(0.5, 2), Exp(2), 2) - 1.1867726850) < 1e-9
    assert abs(_queue_length(Exp(1), Exp(2)) - 1.0) < 1e-9
    assert abs(_queue_length(Exp(3), Exp(2), 2) - 3.4285714286) < 1e-9


def test_me_service_anchors_unchanged():
    # The gate added here is on the ARRIVAL process only.
    assert abs(_queue_length(Exp(0.5), ME.fromErlang(2, 2.0)) - 0.8749999987) < 1e-9
    assert abs(_queue_length(Exp(0.5), ME.fromHyperExp([0.6, 0.4], [2.0, 0.5]))
               - 1.5222222174) < 1e-9
    assert abs(_queue_length(Exp(0.255775446238906), ME(NON_PH_ALPHA, NON_PH_A))
               - 1.0152146751) < 1e-9


def test_ph_service_anchors_unchanged():
    assert abs(_queue_length(Exp(0.5), Erlang.fitMeanAndOrder(0.5, 2)) - 0.3125000000) < 1e-9
    assert abs(_queue_length(Exp(1.2), Erlang.fitMeanAndOrder(0.5, 2), 2) - 0.6964285714) < 1e-9
