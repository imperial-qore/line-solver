"""
Native-Python tests for OPEN stochastic Petri nets (Source -> Place ->
Transition -> Sink) on the SERIAL (Gillespie) SSA engine
(api/solvers/ssa/serial.py). Mirror of the JAR SolverSSASpnOpenSerialTest.

Before the fix the serial engine never fired an open-SPN transition: the gsync
mapping was iterated by key (skipping every firing), the Transition state was
not initialised to [idle, phase, fired], and the Place enabling/consume logic
read the (empty) buffer slot instead of the buffer+server total. Tokens piled
up to the cutoff with zero throughput.

Ground truth is the exact closed form, which SolverJMT reproduces (the JMT
cross-check is in the MATLAB test): an M/M/1-as-SPN has mean tokens rho/(1-rho);
an M/M/inf-as-SPN has mean tokens lambda/mu; a tandem of two single-server
transitions is a Jackson network whose places are each M/M/1 at the shared
arrival rate.

Each test asserts (a) the run used the serial engine (result.method contains
'serial', never a silent NRM fallback) and (b) the marking mean and transition
throughput match the closed form.
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

from line_solver import (Network, Source, Sink, Place, Transition, OpenClass,
                         Exp, SolverSSA)

SAMPLES = 80000
SEED = 23000
RTOL = 0.06


def _single_queue(lam, mu, nserv):
    m = Network('spn_open')
    src = Source(m, 'Source'); snk = Sink(m, 'Sink')
    P1 = Place(m, 'P1'); T1 = Transition(m, 'T1')
    jc = OpenClass(m, 'Class1', 0)
    src.setArrival(jc, Exp(lam))
    a = T1.addMode('Mode1')
    T1.setNumberOfServers(a, nserv)
    T1.setDistribution(a, Exp(mu))
    T1.setEnablingConditions(a, jc, P1, 1)
    T1.setFiringOutcome(a, jc, snk, 1)
    m.link(Network.serialRouting(src, P1, T1, snk))
    return m


def _tandem(lam, mu1, mu2):
    m = Network('spn_tandem')
    src = Source(m, 'Source'); snk = Sink(m, 'Sink')
    P1 = Place(m, 'P1'); P2 = Place(m, 'P2')
    T1 = Transition(m, 'T1'); T2 = Transition(m, 'T2')
    jc = OpenClass(m, 'Class1', 0)
    src.setArrival(jc, Exp(lam))
    a = T1.addMode('Mode1'); T1.setNumberOfServers(a, 1)
    T1.setDistribution(a, Exp(mu1)); T1.setEnablingConditions(a, jc, P1, 1)
    T1.setFiringOutcome(a, jc, P2, 1)
    b = T2.addMode('Mode2'); T2.setNumberOfServers(b, 1)
    T2.setDistribution(b, Exp(mu2)); T2.setEnablingConditions(b, jc, P2, 1)
    T2.setFiringOutcome(b, jc, snk, 1)
    m.link(Network.serialRouting(src, P1, T1, P2, T2, snk))
    return m


def _serial(model):
    return SolverSSA(model, 'serial', seed=SEED, samples=SAMPLES, cutoff=60,
                     verbose=False)


def _assert_ran_serial(solver):
    method = getattr(solver._result, 'method', None)
    assert method is not None and 'serial' in str(method), (
        "open SPN must run the serial method, got: %s" % method)


def _close(expected, got, what):
    denom = max(abs(expected), 1e-9)
    err = abs(got - expected) / denom
    assert err < RTOL, (
        "%s: serial SSA %g deviates from exact %g by %.2f%%"
        % (what, got, expected, 100.0 * err))


def test_mm1_single_server():
    lam, mu = 0.5, 1.0
    rho = lam / mu
    qexact = rho / (1.0 - rho)   # 1.0
    s = _serial(_single_queue(lam, mu, 1))
    at = s.getAvgTable()
    _assert_ran_serial(s)
    _close(qexact, at.QLen.iloc[1], "M/M/1 P1 mean tokens")
    _close(lam, at.Tput.iloc[1], "M/M/1 P1 throughput")


def test_mm1_heavier_load():
    lam, mu = 0.8, 1.0
    rho = lam / mu
    qexact = rho / (1.0 - rho)   # 4.0
    s = _serial(_single_queue(lam, mu, 1))
    at = s.getAvgTable()
    _assert_ran_serial(s)
    _close(qexact, at.QLen.iloc[1], "M/M/1 rho=0.8 P1 mean tokens")
    _close(lam, at.Tput.iloc[1], "M/M/1 rho=0.8 P1 throughput")


def test_mm_infinite_server():
    lam, mu = 2.0, 1.0
    qexact = lam / mu            # 2.0
    s = _serial(_single_queue(lam, mu, float('inf')))
    at = s.getAvgTable()
    _assert_ran_serial(s)
    _close(qexact, at.QLen.iloc[1], "M/M/inf P1 mean tokens")
    _close(lam, at.Tput.iloc[1], "M/M/inf P1 throughput")


def test_tandem_jackson():
    lam, mu1, mu2 = 0.5, 1.0, 0.8
    q1 = (lam / mu1) / (1.0 - lam / mu1)   # 1.0
    q2 = (lam / mu2) / (1.0 - lam / mu2)   # 1.5
    s = _serial(_tandem(lam, mu1, mu2))
    at = s.getAvgTable()
    _assert_ran_serial(s)
    _close(q1, at.QLen.iloc[1], "tandem P1 mean tokens")
    _close(q2, at.QLen.iloc[2], "tandem P2 mean tokens")
    _close(lam, at.Tput.iloc[1], "tandem P1 throughput")
    _close(lam, at.Tput.iloc[2], "tandem P2 throughput")
