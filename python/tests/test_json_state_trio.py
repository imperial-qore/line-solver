"""A declared state is a TRIO on the wire, and a marginal may be FRACTIONAL.

Two contracts that the JSON round-trip of a random environment depends on, and
that nothing else exercised:

1. `state`, `space` and `statePrior` are set together by `init_default`,
   `init_from_marginal` and `init_from_marginal_and_started`. The writer omits a
   trivial `[1]` prior over one row because `initialState` already carries that
   row, so the READER owes the pair back. A node holding a state over an EMPTY
   space is not one any solver can start from -- MATLAB's fluid solver indexes
   `sn.space` directly, and the missing space returned an all-zero table for
   every stage of a reloaded `Environment`.

2. `init_from_marginal` is NOT `init_from_marginal_and_started(n, 0)`. A
   fractional marginal is a fluid initial condition and is kept verbatim;
   routing it through the discrete state-space generator truncates it per
   station and silently drops jobs, so a closed model of 5 was written out
   holding 4.

3. `stages[].type` is written as well as read, so an `Environment` does not come
   back with its stage types blanked.
"""
import json
import os
import tempfile

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Environment, Exp, Network, Queue,
                         SchedStrategy)
from line_solver.io.linemodel_io import load_model, save_model


def _closed_model(name='trio'):
    model = Network(name)
    delay = Delay(model, 'ThinkTime')
    queue = Queue(model, 'Server', SchedStrategy.FCFS)
    jobclass = ClosedClass(model, 'Jobs', 5, delay)
    delay.setService(jobclass, Exp(1.0))
    queue.setService(jobclass, Exp(2.0))
    model.link(Network.serialRouting(delay, queue))
    return model


def _node(model, name):
    for nd in model.get_nodes():
        if nd.get_name() == name:
            return nd
    raise AssertionError('node %s not found' % name)


def _roundtrip(model, tag):
    path = os.path.join(tempfile.mkdtemp(), tag + '.json')
    save_model(model, path)
    return load_model(path), path


def test_a_declared_state_comes_back_with_its_space_and_prior():
    model = _closed_model()
    model.init_from_marginal(np.array([[2.0], [3.0]]))
    back, _ = _roundtrip(model, 'state_trio')
    for name in ('ThinkTime', 'Server'):
        nd = _node(back, name)
        state = np.atleast_1d(np.asarray(nd.state, dtype=float)).ravel()
        space = np.atleast_2d(np.asarray(nd.get_state_space(), dtype=float))
        prior = np.atleast_1d(np.asarray(nd.get_state_prior(), dtype=float)).ravel()
        assert space.shape[0] == 1, '%s: a declared state is a ONE-ROW space' % name
        assert np.allclose(space[0], state), '%s: space row is not the state' % name
        assert np.allclose(prior, [1.0]), '%s: prior is not [1]' % name


def test_an_explicit_prior_is_not_overwritten_by_the_trivial_one():
    model = _closed_model('explicit')
    model.init_from_marginal(np.array([[3.0], [2.0]]))
    queue = _node(model, 'Server')
    space = np.array([[1.0, 1.0], [0.0, 1.0]])
    queue.set_state_space(space)
    queue.setStatePrior(np.array([0.25, 0.75]))
    queue.set_state(space[0])
    back, _ = _roundtrip(model, 'explicit_prior')
    nd = _node(back, 'Server')
    assert np.allclose(np.asarray(nd.get_state_space(), dtype=float), space)
    assert np.allclose(np.asarray(nd.get_state_prior(), dtype=float).ravel(),
                       [0.25, 0.75])


def test_a_fractional_marginal_is_kept_and_conserves_the_population():
    model = _closed_model('fluid')
    model.init_from_marginal(np.array([[2.5], [2.5]]))
    total = 0.0
    for name in ('ThinkTime', 'Server'):
        state = np.atleast_1d(np.asarray(_node(model, name).state, dtype=float)).ravel()
        assert state.size == 1 and not float(state[0]).is_integer(), \
            '%s: a purposely fractional marginal was rounded away' % name
        total += float(state[0])
    assert total == pytest.approx(5.0), 'the closed population was not conserved'


def test_a_marginal_that_is_no_state_is_refused():
    model = _closed_model('bad')
    with pytest.raises(ValueError):
        model.init_from_marginal(np.array([[2.0], [2.0]]))


def _environment():
    env = Environment('ServerModes', 2)
    stages = (('Fast', 'operational', 4.0), ('Slow', 'degraded', 1.0))
    for idx, (stage, stype, rate) in enumerate(stages):
        model = _closed_model(stage)
        _node(model, 'Server').setService(model.get_classes()[0], Exp(rate))
        env.add_stage(idx, stage, stype, model)
    env.add_transition(0, 1, Exp(0.5))
    env.add_transition(1, 0, Exp(1.0))
    return env


def test_environment_stage_type_survives_the_round_trip():
    env = _environment()
    path = os.path.join(tempfile.mkdtemp(), 'env_stage_type.json')
    save_model(env, path)
    with open(path) as fh:
        doc = json.load(fh)
    stages = doc['model']['stages']
    assert [s.get('type') for s in stages] == ['operational', 'degraded'], \
        'the writer dropped stages[].type'
    back = load_model(path)
    assert back._stage_types == ['operational', 'degraded'], \
        'the reader dropped stages[].type'
