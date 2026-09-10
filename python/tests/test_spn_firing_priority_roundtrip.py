"""
Regression test for the JSON round trip of a Transition's firing priority.

Zero is a legal JMT firing priority, and the builder default is 1. The three
model.json writers used to guard the key with ``> 0``, so an explicitly set 0 was
omitted from the file and every reader restored the default: the guard dropped a
meaningful value and preserved the redundant one. The correct guard is ``!= 1``,
matching the neighbouring firingWeight. See BUGS.md BUG-90.
"""

import json
import os

import pytest

import line_solver
# The worktree copy of line_solver must be the one under test; a global install
# would silently validate the wrong code.
assert os.path.realpath(__file__).rsplit('/python/', 1)[0] in os.path.realpath(
    line_solver.__file__), (
    'line_solver resolves outside this worktree: %s' % line_solver.__file__)

from line_solver import (Exp, GlobalConstants, Network, OpenClass, Place, Sink,
                         Source, Transition)
from line_solver.io.linemodel_io import load_model, save_model


def _build(priority):
    model = Network('spnprio')
    source = Source(model, 'Source')
    sink = Sink(model, 'Sink')
    place = Place(model, 'P1')
    trans = Transition(model, 'T1')
    jobclass = OpenClass(model, 'Class1', 0)
    source.set_arrival(jobclass, Exp(1.0))
    mode = trans.add_mode('Mode1')
    trans.set_number_of_servers(mode, GlobalConstants.MaxInt)
    trans.set_distribution(mode, Exp(4.0))
    trans.set_enabling_conditions(mode, jobclass, place, 1)
    trans.set_firing_outcome(mode, jobclass, sink, 1)
    trans.set_firing_priorities(mode, priority)
    routing = model.init_routing_matrix()
    routing.set(jobclass, jobclass, source, place, 1.0)
    routing.set(jobclass, jobclass, place, trans, 1.0)
    routing.set(jobclass, jobclass, trans, sink, 1.0)
    model.link(routing)
    return model


@pytest.mark.parametrize('priority', [0, 1, 3])
def test_firing_priority_survives_the_round_trip(tmp_path, priority):
    path = str(tmp_path / ('spnprio_%s.json' % priority))
    save_model(_build(priority), path)

    doc = json.load(open(path))
    mode = [n for n in doc['model']['nodes'] if n['name'] == 'T1'][0]['modes'][0]
    # The key travels iff it differs from the builder default of 1.
    assert ('firingPriority' in mode) == (priority != 1)

    back = load_model(path)
    trans = [n for n in back.get_nodes() if n.get_name() == 'T1'][0]
    assert float(trans._firing_priorities[0]) == float(priority)
