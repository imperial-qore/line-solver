"""
Regression tests: every native solver's supports() must consult its feature set.

MATLAB gates each solver's supports() against getFeatureSet(). Several Python
solvers did not: JMT returned a literal True ("for now return True"), QNS and
SSA checked only the station and class counts, ENV guarded on
isinstance(used, set) which is never true (getUsedLangFeatures returns a
SolverFeatureSet) and so always fell through to True, and LDES/LQNS defined no
supports() anywhere in the MRO, so calling it raised AttributeError. AUTO
returned True whenever any candidate merely existed, overriding the verdict of
its own loop.

A gate must reject the FEATURE, not the model class, so each solver is also
checked against a model it does support.
"""
import warnings

import pytest

from line_solver import (Exp, Network, OpenClass, Queue, SchedStrategy, Sink,
                         SolverAuto, SolverJMT, SolverLDES, SolverLQNS,
                         SolverQNS, SolverSSA, Source)


def _plain():
    model = Network('plain')
    source = Source(model, 'S')
    queue = Queue(model, 'Q', SchedStrategy.FCFS)
    sink = Sink(model, 'K')
    oclass = OpenClass(model, 'C', 0)
    source.setArrival(oclass, Exp(0.5))
    queue.setService(oclass, Exp(1.0))
    model.link(Network.serialRouting(source, queue, sink))
    return model


def _fcfspr():
    model = Network('fcfspr')
    source = Source(model, 'S')
    queue = Queue(model, 'Q', SchedStrategy.FCFSPR)
    sink = Sink(model, 'K')
    oclass = OpenClass(model, 'C', 0)
    source.setArrival(oclass, Exp(0.5))
    queue.setService(oclass, Exp(1.0))
    model.link(Network.serialRouting(source, queue, sink))
    return model


def _supports(solver, model):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return solver.supports(model)


@pytest.mark.parametrize('solver', [SolverJMT, SolverQNS, SolverSSA, SolverLQNS])
def test_supports_rejects_an_unsupported_feature(solver):
    # None of these declare SchedStrategy_FCFSPR.
    assert 'SchedStrategy_FCFSPR' not in solver.getFeatureSet()
    assert _supports(solver, _fcfspr()) is False


@pytest.mark.parametrize('solver', [SolverJMT, SolverQNS, SolverSSA])
def test_supports_accepts_a_plain_model(solver):
    # The gate must reject the feature, not the model class.
    assert _supports(solver, _plain()) is True


def test_ldes_supports_exists_and_honours_its_featureset():
    # LDES had no supports() in the MRO at all -> AttributeError.
    assert hasattr(SolverLDES, 'supports')
    assert 'SchedStrategy_FCFSPR' in SolverLDES.getFeatureSet()
    assert _supports(SolverLDES, _fcfspr()) is True


def test_auto_follows_its_candidate_loop():
    # AUTO returned True whenever a candidate existed, regardless of support.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert SolverAuto(_plain()).supports(_plain()) is True
        # LDES is a candidate and supports FCFSPR, so AUTO must accept it. This
        # also pins the candidate construction: LDES/JMT moved under wrappers/
        # and AUTO kept the old import path, so _create_solver raised
        # ModuleNotFoundError and every LDES/JMT candidate was silently skipped.
        assert SolverAuto(_fcfspr()).supports(_fcfspr()) is True


def test_auto_can_construct_its_candidates():
    auto = SolverAuto(_plain())
    for candidate in auto._candidate_solvers:
        auto._create_solver(candidate)  # must not raise ModuleNotFoundError
