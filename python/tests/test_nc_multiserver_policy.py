"""Regression tests for SolverNC's ``config.multiserver`` handling.

SolverNC represents a finite multiserver station in one of two ways: Seidmann's
approximation (demand ``L/c`` plus a delay ``L(c-1)/c``, applied in
``solver_nc``) or the exact load-dependent lattice ``mu(n)=min(n,c)``, which
routes the model to ``solver_ncld``. Which one it uses used to be decided by
the method name alone -- Seidmann on ``default``, the lattice on ``exact`` --
with no way to ask for either, even though ``config.multiserver`` exists on the
shared solver options and SolverMVA honours it.

The load-bearing property here is that the SHIPPED DEFAULT MOVES NO RESULT.
``config.multiserver='default'``, and an options object that never mentions the
field, must both reproduce the historical dispatch exactly; only an explicit
``'seidmann'`` or ``'lld'`` changes anything. The tests below pin that, using
Bolch, Greiner, de Meer & Trivedi Example 8.1 (three FCFS stations with 2, 3
and 1 servers, one closed class of 3 jobs), whose exact load-dependent answer
is known in closed form and is what SolverMVA('exact') returns.
"""
import numpy as np
import pytest

from line_solver import (ClosedClass, Exp, Network, Queue, SchedStrategy,
                         SolverMVA, SolverNC)

# Seidmann's approximation on this model, i.e. what the default path returns.
SEIDMANN_QLEN = [1.35037, 1.05217, 0.59745]
# The exact load-dependent product form, mu(n)=min(n,c). Reproduced independently
# by enumerating the population lattice and by SolverMVA('exact').
EXACT_QLEN = [1.28981, 1.04988, 0.66031]


def _multiserver_model():
    """Bolch Example 8.1: three FCFS stations, 2/3/1 servers, K=3."""
    model = Network('Bolch_08_01')
    node1 = Queue(model, 'Node1', SchedStrategy.FCFS)
    node1.setNumberOfServers(2)
    node2 = Queue(model, 'Node2', SchedStrategy.FCFS)
    node2.setNumberOfServers(3)
    node3 = Queue(model, 'Node3', SchedStrategy.FCFS)
    jobclass = ClosedClass(model, 'Class1', 3, node1, 0)
    node1.setService(jobclass, Exp(0.8))
    node2.setService(jobclass, Exp(0.6))
    node3.setService(jobclass, Exp(0.4))
    P = model.initRoutingMatrix()
    P[jobclass, jobclass] = np.array([[0.133, 0.667, 0.2], [1.0, 0, 0], [1.0, 0, 0]])
    model.link(P)
    return model


def _single_server_model():
    """Bolch Example 7.5: the same shape with one server everywhere."""
    model = Network('Bolch_07_05')
    nodes = [Queue(model, 'Node%d' % i, SchedStrategy.FCFS) for i in (1, 2, 3)]
    jobclass = ClosedClass(model, 'Class1', 3, nodes[0], 0)
    for node, rate in zip(nodes, (0.8, 0.6, 0.4)):
        node.setService(jobclass, Exp(rate))
    P = model.initRoutingMatrix()
    P[jobclass, jobclass] = np.array([[0.6, 0.3, 0.1], [0.2, 0.3, 0.5], [0.4, 0.0, 0.6]])
    model.link(P)
    return model


def _qlen(model, method=None, multiserver=None):
    options = SolverNC.defaultOptions()
    if method is not None:
        options['method'] = method
    if multiserver is not None:
        options['config'] = {'multiserver': multiserver}
    return np.asarray(SolverNC(model, options).getAvgTable().QLen, dtype=float)


@pytest.mark.parametrize('multiserver', [None, 'default'])
def test_default_policy_is_unchanged(multiserver):
    """The shipped default must not move: still Seidmann on 'default'."""
    got = _qlen(_multiserver_model(), multiserver=multiserver)
    np.testing.assert_allclose(got, SEIDMANN_QLEN, rtol=1e-4)


@pytest.mark.parametrize('multiserver', [None, 'default'])
def test_exact_method_is_unchanged(multiserver):
    """'exact' keeps converting to the load-dependent lattice."""
    got = _qlen(_multiserver_model(), method='exact', multiserver=multiserver)
    np.testing.assert_allclose(got, EXACT_QLEN, rtol=1e-4)


def test_lld_makes_the_default_exact():
    """config.multiserver='lld' reaches the exact answer from the default method."""
    got = _qlen(_multiserver_model(), multiserver='lld')
    np.testing.assert_allclose(got, EXACT_QLEN, rtol=1e-4)


def test_seidmann_opts_the_exact_method_out():
    """config.multiserver='seidmann' keeps Seidmann even when 'exact' is asked for."""
    got = _qlen(_multiserver_model(), method='exact', multiserver='seidmann')
    np.testing.assert_allclose(got, SEIDMANN_QLEN, rtol=1e-4)


def test_lld_matches_mva_exact():
    """The opt-in answer is SolverMVA's exact multiserver answer, not a near miss."""
    got = _qlen(_multiserver_model(), multiserver='lld')
    mva = np.asarray(SolverMVA(_multiserver_model(), 'exact').getAvgTable().QLen, dtype=float)
    np.testing.assert_allclose(got, mva, rtol=1e-6)


def test_unimplemented_mva_value_falls_back_to_default():
    """A SolverMVA-only rule warns and is not silently honoured as something else."""
    got = _qlen(_multiserver_model(), multiserver='softmin')
    np.testing.assert_allclose(got, SEIDMANN_QLEN, rtol=1e-4)


@pytest.mark.parametrize('method', [None, 'exact'])
@pytest.mark.parametrize('multiserver', [None, 'default', 'lld', 'seidmann'])
def test_single_server_model_is_policy_independent(method, multiserver):
    """With one server everywhere the policy has nothing to select, under any method."""
    baseline = _qlen(_single_server_model())
    got = _qlen(_single_server_model(), method=method, multiserver=multiserver)
    np.testing.assert_allclose(got, baseline, rtol=1e-9)
