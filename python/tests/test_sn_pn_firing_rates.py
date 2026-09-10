"""Decline paths of the Place firing-rate recovery.

The departure equation may only reference TIMED modes, because an immediate
firing takes zero time and is never counted in TN. A Place drained only by
immediate modes therefore contributes no departure row. When NO Place
contributes one the system is homogeneous, pinv returns the zero vector, and
reporting it would mark every Place idle; the recovery must decline instead so
the caller keeps the throughputs it already had.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver import (Network, Place, Transition, ClosedClass, Exp, Immediate,
                         TimingStrategy)
from line_solver.api.sn import sn_pn_firing_rates


def _cycle(t1_immediate):
    """Two Places in a cycle. T2 is always immediate; T1 follows the flag."""
    model = Network('spn_all_immediate' if t1_immediate else 'spn_one_timed')
    P1 = Place(model, 'P1')
    P2 = Place(model, 'P2')
    T1 = Transition(model, 'T1')
    T2 = Transition(model, 'T2')
    jobclass = ClosedClass(model, 'Class1', 1, P1, 0)

    m1 = T1.addMode('M1')
    if t1_immediate:
        T1.setDistribution(m1, Immediate())
        T1.setTimingStrategy(m1, TimingStrategy.IMMEDIATE)
        T1.setFiringPriorities(m1, 1)
        T1.setFiringWeights(m1, 1.0)
    else:
        T1.setDistribution(m1, Exp(3))
    T1.setEnablingConditions(m1, jobclass, P1, 1)
    T1.setFiringOutcome(m1, jobclass, P2, 1)

    m2 = T2.addMode('M2')
    T2.setDistribution(m2, Immediate())
    T2.setTimingStrategy(m2, TimingStrategy.IMMEDIATE)
    T2.setFiringPriorities(m2, 1)
    T2.setFiringWeights(m2, 1.0)
    T2.setEnablingConditions(m2, jobclass, P2, 1)
    T2.setFiringOutcome(m2, jobclass, P1, 1)

    rm = model.initRoutingMatrix()
    rm.set(jobclass, jobclass, P1, T1, 1.0)
    rm.set(jobclass, jobclass, P2, T2, 1.0)
    rm.set(jobclass, jobclass, T1, P2, 1.0)
    rm.set(jobclass, jobclass, T2, P1, 1.0)
    model.link(rm)

    P1.setState(1)
    P2.setState(0)
    return model


def test_declines_when_every_place_is_drained_only_by_immediate_modes():
    sn = _cycle(t1_immediate=True).getStruct()
    # Whatever the analyzer measured is uninformative here: an immediate firing
    # is not a timed event, so these entries carry no firing rate.
    TN = np.zeros((sn.nstations, sn.nclasses))
    x, _consumed, _produced, _place_nodes = sn_pn_firing_rates(sn, TN, False)
    assert x is None or len(x) == 0, (
        'recovery must decline when no timed mode drains any Place, otherwise '
        'the homogeneous system returns all-zero rates, got %r' % (x,))


def test_recovers_when_at_least_one_timed_mode_drains_a_place():
    sn = _cycle(t1_immediate=False).getStruct()
    TN = np.full((sn.nstations, sn.nclasses), 1.5)
    x, _consumed, _produced, _place_nodes = sn_pn_firing_rates(sn, TN, False)
    assert x is not None and len(x) == 2, (
        'one timed departure row is enough to determine the rates, got %r' % (x,))
    # T1 is the only measured mode and drains P1 at 1.5; the balance equations
    # carry that through to the immediate T2.
    assert np.allclose(x, [1.5, 1.5], atol=1e-9), (
        'timed rate must equal the measured throughput and the immediate rate '
        'is pinned by token balance, got %r' % (x,))
