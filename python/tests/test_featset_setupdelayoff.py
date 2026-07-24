"""
Regression test: the SetupDelayOff feature gate.

Setup/delay-off is honoured only by MAM, JMT and LDES. No featset entry
declared it, so CTMC/SSA/MVA/NC/FLD accepted a model built with setDelayOff and
solved it setup-free, silently returning the plain M/M/1 answer. The feature is
now registered and detected in getUsedLangFeatures, so those solvers reject
instead.

SolverSSA.supports was additionally a stub that checked only the station and
class counts, so it accepted every model regardless of the features used and
could not be gated on anything.
"""
import warnings

import pytest

from line_solver import (Exp, Network, OpenClass, Queue, SchedStrategy, Sink,
                         SolverCTMC, SolverMVA, SolverNC, SolverSSA, Source)


def _build(with_setup):
    model = Network('setup_featset')
    source = Source(model, 'Source')
    queue = Queue(model, 'Q', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'C', 0)
    source.setArrival(oclass, Exp(0.5))
    queue.setService(oclass, Exp(1.0))
    if with_setup:
        queue.setDelayOff(oclass, Exp(2.0), Exp(4.0))
    model.link(Network.serialRouting(source, queue, sink))
    return model


def _supports(solver, model):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return solver.supports(model)


def test_feature_is_detected():
    used = _build(True).getUsedLangFeatures()
    assert used.list['SetupDelayOff'] is True
    assert _build(False).getUsedLangFeatures().list['SetupDelayOff'] is False


@pytest.mark.parametrize('solver', [SolverCTMC, SolverSSA, SolverMVA, SolverNC])
def test_solvers_without_setup_support_reject(solver):
    assert _supports(solver, _build(True)) is False


@pytest.mark.parametrize('solver', [SolverCTMC, SolverSSA, SolverMVA, SolverNC])
def test_same_solvers_still_accept_a_plain_model(solver):
    # The gate must reject setup, not the whole model class.
    assert _supports(solver, _build(False)) is True
