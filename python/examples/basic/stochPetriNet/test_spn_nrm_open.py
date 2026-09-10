"""
test_spn_nrm_open

Validates the SSA Next-Reaction-Method (NRM) stochastic-Petri-net path on OPEN
nets: a Source feeds a Place whose tokens drain through a Transition to a Sink.
Before the Source-arrival reaction was added the fed Place stayed empty and the
run threw "Deadlock: no transition is enabled".

Each net asserts (a) the solver actually ran method 'nrm' (never a silent serial
fallback) and (b) the simulated marking mean and throughput match the analytic
M/M/1 result. A Source Exp(lambda) feeding a single-server Transition Exp(mu) is
an M/M/1 queue at the Place: mean tokens = rho/(1-rho), throughput = lambda
(rho = lambda/mu < 1). The canonical net is cross-checked against JMT.
"""

import numpy as np

from line_solver import (Exp, GlobalConstants, JMT, Network, OpenClass, Place,
                         SSA, Sink, Source, Transition, VerboseLevel)

RTOL = 0.04
SAMPLES = int(3e5)
SEED = 23000


def mm1spn(lam, mu):
    model = Network('mm1spn')
    source = Source(model, 'Source')
    sink = Sink(model, 'Sink')
    P1 = Place(model, 'P1')
    T1 = Transition(model, 'T1')
    jobclass = OpenClass(model, 'Class1', 0)
    source.setArrival(jobclass, Exp(lam))
    mode = T1.addMode('Mode1')
    T1.setDistribution(mode, Exp(mu))
    T1.setEnablingConditions(mode, jobclass, P1, 1)
    T1.setFiringOutcome(mode, jobclass, sink, 1)
    R = model.initRoutingMatrix()
    R.set(jobclass, jobclass, source, P1, 1.0)
    R.set(jobclass, jobclass, P1, T1, 1.0)
    R.set(jobclass, jobclass, T1, sink, 1.0)
    model.link(R)
    return model


def tandemspn(lam, mu1, mu2):
    model = Network('tandemspn')
    source = Source(model, 'Source')
    sink = Sink(model, 'Sink')
    P1 = Place(model, 'P1')
    P2 = Place(model, 'P2')
    T1 = Transition(model, 'T1')
    T2 = Transition(model, 'T2')
    jobclass = OpenClass(model, 'Class1', 0)
    source.setArrival(jobclass, Exp(lam))
    m1 = T1.addMode('Mode1')
    T1.setDistribution(m1, Exp(mu1))
    T1.setEnablingConditions(m1, jobclass, P1, 1)
    T1.setFiringOutcome(m1, jobclass, P2, 1)
    m2 = T2.addMode('Mode1')
    T2.setDistribution(m2, Exp(mu2))
    T2.setEnablingConditions(m2, jobclass, P2, 1)
    T2.setFiringOutcome(m2, jobclass, sink, 1)
    R = model.initRoutingMatrix()
    R.set(jobclass, jobclass, source, P1, 1.0)
    R.set(jobclass, jobclass, P1, T1, 1.0)
    R.set(jobclass, jobclass, T1, P2, 1.0)
    R.set(jobclass, jobclass, P2, T2, 1.0)
    R.set(jobclass, jobclass, T2, sink, 1.0)
    model.link(R)
    return model


def _station_totals(mat):
    """Row sums of a (station x class) matrix, so a one-class net reads as a vector."""
    return np.asarray(mat, dtype=float).sum(axis=1)


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    # --- Net 1: M/M/1 SPN, Source Exp(0.5) -> P1 -> T1 Exp(1.0) -> Sink ----
    lam, mu = 0.5, 1.0
    rho = lam / mu
    q_exact = rho / (1 - rho)          # = 1.0
    solver = SSA(mm1spn(lam, mu), method='nrm', samples=SAMPLES, seed=SEED)
    avg = solver.getAvg()
    Qn, Tn = _station_totals(avg[0]), _station_totals(avg[3])
    assert 'nrm' in str(solver.result.method), 'net1 did not run NRM'
    # stations: Source(0), P1(1)
    assert abs(Qn[1] - q_exact) / q_exact < RTOL, \
        'net1 P1 tokens %g vs exact %g' % (Qn[1], q_exact)
    assert abs(Tn[0] - lam) / lam < RTOL, 'net1 Source tput %g vs %g' % (Tn[0], lam)
    assert abs(Tn[1] - lam) / lam < RTOL, 'net1 P1 tput %g vs %g' % (Tn[1], lam)

    # Cross-check mean tokens against JMT's simulation of the same net.
    Qj = _station_totals(JMT(mm1spn(lam, mu), samples=SAMPLES, seed=SEED).getAvg()[0])
    assert abs(Qn[1] - Qj[1]) / max(Qj[1], 1e-9) < RTOL, \
        'net1 P1 tokens NRM %g vs JMT %g' % (Qn[1], Qj[1])

    # --- Net 2: open tandem, two places in series -------------------------
    lam, mu1, mu2 = 0.5, 1.0, 2.0
    q1 = (lam / mu1) / (1 - lam / mu1)   # = 1.0
    q2 = (lam / mu2) / (1 - lam / mu2)   # = 1/3
    solver = SSA(tandemspn(lam, mu1, mu2), method='nrm', samples=SAMPLES, seed=SEED)
    avg = solver.getAvg()
    Qn, Tn = _station_totals(avg[0]), _station_totals(avg[3])
    assert 'nrm' in str(solver.result.method), 'net2 did not run NRM'
    # stations: Source(0), P1(1), P2(2)
    assert abs(Qn[1] - q1) / q1 < RTOL, 'net2 P1 tokens %g vs exact %g' % (Qn[1], q1)
    assert abs(Qn[2] - q2) / q2 < RTOL, 'net2 P2 tokens %g vs exact %g' % (Qn[2], q2)
    assert abs(Tn[1] - lam) / lam < RTOL, 'net2 P1 tput %g vs %g' % (Tn[1], lam)
    assert abs(Tn[2] - lam) / lam < RTOL, 'net2 P2 tput %g vs %g' % (Tn[2], lam)

    print('test_spn_nrm_open passed')
