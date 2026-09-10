"""Regression tests for the language-feature registry entries no model could set.

``supports`` iterates the registry, so a name ``get_used_lang_features`` never
emits gates nothing: a solver declaring it says nothing about what it accepts,
and a solver omitting it refuses nothing. ``Cox2`` and ``Trace`` were two such
names in MATLAB, the JAR and C++ (see _kb/06-solver-catalog.md), shadowed by
the parent distribution's name.

Native Python was the codebase where the OPPOSITE half of the defect had already
bitten: it marked ``Cox2``, but no generalization fallback existed, so
``SolverMVA`` REFUSED a Cox2 model it solved happily when the same distribution
was spelled ``Coxian``. Both halves are pinned here.
"""

import numpy as np
import pytest

from line_solver import (Coxian, Cox2, Exp, Network, OpenClass, Queue, Replayer,
                         SchedStrategy, Sink, Source, SolverCTMC, SolverFluid,
                         SolverMVA, SolverNC, Trace)
from line_solver.solvers.base import SolverFeatureSet

TRACE_SAMPLES = np.array([0.5, 1.0, 1.5, 2.0])


def _model_with(dist):
    """Source -> FCFS queue -> Sink, one open class, service set to dist."""
    model = Network('featreach')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    jobclass = OpenClass(model, 'Class1')
    source.setArrival(jobclass, Exp(0.2))
    queue.setService(jobclass, dist)
    model.link(Network.serialRouting(source, queue, sink))
    return model


def _marks(dist):
    return {name for name, on in _model_with(dist).getUsedLangFeatures().list.items() if on}


# --------------------------------------------------------------------------
# emission
# --------------------------------------------------------------------------

def test_two_phase_coxian_marks_cox2():
    # MATLAB has no Cox2 OBJECT -- Cox2.fitMeanAndSCV returns a Coxian -- so the
    # phase count is the only reading of the entry the four codebases share.
    marks = _marks(Cox2.fitMeanAndSCV(1.0, 2.0))
    assert 'Cox2' in marks
    assert 'Coxian' not in marks

    # The class the user happened to call is not the test: a generic Coxian with
    # two phases is the same distribution and marks the same entry.
    marks = _marks(Coxian([0.5, 0.5], [0.4, 1.0]))
    assert 'Cox2' in marks
    assert 'Coxian' not in marks


def test_coxian_beyond_two_phases_marks_coxian():
    marks = _marks(Coxian([0.4, 0.3, 0.3], [0.3, 0.5, 1.0]))
    assert 'Coxian' in marks
    assert 'Cox2' not in marks


def test_trace_marks_trace_and_replayer_marks_replayer():
    assert 'Trace' in _marks(Trace(TRACE_SAMPLES))
    assert 'Replayer' not in _marks(Trace(TRACE_SAMPLES))
    assert 'Replayer' in _marks(Replayer(TRACE_SAMPLES))
    assert 'Trace' not in _marks(Replayer(TRACE_SAMPLES))


def test_name_is_unchanged_so_the_wire_type_is_unchanged():
    # get_feature_name exists precisely so that name can stay put: it selects the
    # ProcessType and the JSON wire type, and moving it would change what a saved
    # model reloads as.
    assert Cox2.fitMeanAndSCV(1.0, 2.0).name == 'Cox2'
    assert Coxian([0.5, 0.5], [0.4, 1.0]).name == 'Coxian'
    assert Coxian([0.5, 0.5], [0.4, 1.0]).get_feature_name() == 'Cox2'
    assert Trace(TRACE_SAMPLES).name == 'Trace'


# --------------------------------------------------------------------------
# the gate
# --------------------------------------------------------------------------

def test_specialization_falls_back_to_its_generalization():
    assert SolverFeatureSet.GENERALIZATION_OF == {'Cox2': 'Coxian', 'Trace': 'Replayer'}, \
        'a fallback WIDENS what a declared set accepts; none may be added silently'

    used = SolverFeatureSet()
    used.set_true('Cox2')
    general = SolverFeatureSet()
    general.set_true('Coxian')
    assert SolverFeatureSet.unsupported_features(general, used) == []
    assert SolverFeatureSet.supports_with_reason(general, used)[0]

    # One way only: declaring the SPECIAL case does not buy the general one.
    specific = SolverFeatureSet()
    specific.set_true('Cox2')
    used_general = SolverFeatureSet()
    used_general.set_true('Coxian')
    assert SolverFeatureSet.unsupported_features(specific, used_general) == ['Coxian']

    # Neither declared is still a refusal, naming the specific feature.
    assert SolverFeatureSet.unsupported_features(SolverFeatureSet(), used) == ['Cox2']


@pytest.mark.parametrize('solver', [SolverMVA, SolverNC, SolverCTMC, SolverFluid])
def test_solvers_still_solve_a_two_phase_coxian(solver):
    # The defect this pins: SolverMVA declares Coxian and not Cox2, so before the
    # fallback existed it refused this exact model with
    # "features not supported by the chosen solver (feature: Cox2)".
    assert solver(_model_with(Cox2.fitMeanAndSCV(1.0, 2.0))).getAvgTable() is not None


def test_every_coxian_declaring_featset_accepts_cox2():
    cox2 = SolverFeatureSet()
    cox2.set_true('Cox2')
    trace = SolverFeatureSet()
    trace.set_true('Trace')
    for declared in (SolverMVA.getFeatureSet(), SolverNC.getFeatureSet(),
                     SolverCTMC.getFeatureSet(), SolverFluid.getFeatureSet()):
        # getFeatureSet returns a name set for some solvers here and a
        # SolverFeatureSet for others; both are accepted.
        featset = declared
        if not isinstance(declared, SolverFeatureSet):
            featset = SolverFeatureSet()
            featset.set_true(sorted(declared))
        if featset.list.get('Coxian', False):
            assert SolverFeatureSet.unsupported_features(featset, cox2) == []
        if featset.list.get('Replayer', False):
            assert SolverFeatureSet.unsupported_features(featset, trace) == []


def test_the_lqn_cache_names_are_registered():
    # CacheTask, ItemEntry and ActivityPrecedence_POST_CACHE are declared by the
    # JAR SolverLDES.getLNFeatureSet, where an unregistered name is a hard error:
    # that whole feature set threw before it could compare anything. The four
    # registries must agree, so they are registered here too.
    registry = SolverFeatureSet()
    for name in ('CacheTask', 'ItemEntry', 'ActivityPrecedence_POST_CACHE'):
        assert name in registry.list
        registry.set_true(name)  # set_true raises on an unregistered name


# --------------------------------------------------------------------------
# the last two names the C++ port carried alone (aligned 2026-08-22)
# --------------------------------------------------------------------------

def test_heterogeneous_server_pools_are_marked_and_gated():
    # A pooled station is NOT a station of the same total size: the pools carry
    # their own class compatibilities and their own rates. Every solver that
    # reads only sn.nservers therefore has to refuse rather than flatten them,
    # which it can only do once the pools are marked.
    from line_solver import ServerType, SolverJMT

    model = Network('hetfeat')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    jobclass = OpenClass(model, 'Class1')
    source.setArrival(jobclass, Exp(0.2))
    queue.setService(jobclass, Exp(1.0))
    model.link(Network.serialRouting(source, queue, sink))
    assert not model.getUsedLangFeatures().list['HeteroServers']

    queue.addServerType(ServerType('fast', 1, [jobclass]))
    queue.addServerType(ServerType('slow', 2, [jobclass]))
    used = model.getUsedLangFeatures()
    assert used.list['HeteroServers']

    # JMT serialises the pools and the LDES engine simulates them; the
    # analytical solvers read nservers and must refuse.
    assert 'HeteroServers' in set(SolverJMT.getFeatureSet())
    assert not _declared(SolverMVA).list['HeteroServers']
    assert not _declared(SolverNC).list['HeteroServers']
    assert not _declared(SolverCTMC).list['HeteroServers']

    assert SolverFeatureSet.unsupported_features(_declared(SolverMVA), used) == \
        ['HeteroServers']


def test_non_normal_departure_discipline_is_refused_by_every_solver():
    # A FIFO depository releases a served token only after the tokens that
    # entered service before it, so which output transitions are enabled depends
    # on the arrival order and not only on the marking. NO engine in ANY
    # codebase implements it, so no feature set declares it and the model is
    # refused rather than served as if it were NORMAL.
    from line_solver import ClosedClass, DepartureDiscipline, Place, Transition
    from line_solver import SolverJMT, SolverSSA

    model = Network('depfeat')
    place = Place(model, 'P1')
    Transition(model, 'T1')
    jobclass = ClosedClass(model, 'Class1', 1, place)
    place.setService(jobclass, Exp(1.0))
    assert not model.getUsedLangFeatures().list['DepartureDiscipline']

    place.setDepartureDiscipline(jobclass, DepartureDiscipline.FIFO)
    used = model.getUsedLangFeatures()
    assert used.list['DepartureDiscipline']

    for solver in (SolverJMT, SolverSSA, SolverCTMC, SolverMVA):
        featset = _declared(solver)
        assert not featset.list['DepartureDiscipline'], \
            'no solver implements a FIFO depository, so none may declare it'
        assert 'DepartureDiscipline' in \
            SolverFeatureSet.unsupported_features(featset, used)


def _declared(solver):
    """A solver's declaration as a SolverFeatureSet, whichever shape it returns."""
    declared = solver.getFeatureSet()
    if isinstance(declared, SolverFeatureSet):
        return declared
    featset = SolverFeatureSet()
    featset.set_true(sorted(declared))
    return featset
