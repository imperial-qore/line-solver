"""
Wrapper plumbing for the discrete-time LDES features.

MATLAB and Python reach the LDES engine only through model.json plus CLI flags,
so anything that does not serialize or is not emitted is a feature those
codebases cannot use. These tests cover both channels:

  - the arrivalBatch key survives save_model / load_model;
  - the --slotted / --slotlength flags are emitted from LDESOptions.

The numerical validation of the slotted engine itself lives in the JAR suites
(SolverLDESGeoGeo1Test, SolverLDESGeoXGeo1Test); here the concern is only that
the request reaches it intact.
"""

import json
import os
import tempfile

import pytest

from line_solver import Network, Source, Queue, Sink, OpenClass, SchedStrategy
from line_solver.distributions.discrete import Geometric
from line_solver.io.linemodel_io import save_model, load_model
from line_solver.solvers.wrappers.solver_ldes.ldes_options import LDESOptions


def _geo_x_model(a=0.1, beta=0.5, s=0.9):
    m = Network('geox')
    src = Source(m, 'S')
    q = Queue(m, 'Q', SchedStrategy.FCFS)
    sk = Sink(m, 'K')
    c = OpenClass(m, 'C')
    src.setArrival(c, Geometric(a))
    src.setArrivalBatch(c, Geometric(beta))
    q.setService(c, Geometric(s))
    P = m.initRoutingMatrix()
    P.set(c, c, src, q, 1.0)
    P.set(c, c, q, sk, 1.0)
    m.link(P)
    return m, src, c


def _find_source_json(obj):
    if isinstance(obj, dict):
        if obj.get('type') == 'Source':
            return obj
        for v in obj.values():
            found = _find_source_json(v)
            if found is not None:
                return found
    if isinstance(obj, list):
        for v in obj:
            found = _find_source_json(v)
            if found is not None:
                return found
    return None


# ----------------------------------------------------------------------
# arrivalBatch over the model JSON
# ----------------------------------------------------------------------

def test_arrival_batch_is_serialized():
    m, _, _ = _geo_x_model()
    path = tempfile.mktemp(suffix='.json')
    try:
        save_model(m, path)
        with open(path) as fh:
            doc = json.load(fh)
        src_json = _find_source_json(doc)
        assert src_json is not None
        assert src_json.get('arrivalBatch') == {
            'C': {'type': 'Geometric', 'params': {'p': 0.5}}
        }
    finally:
        if os.path.exists(path):
            os.remove(path)


def test_arrival_batch_survives_round_trip():
    m, _, _ = _geo_x_model()
    path = tempfile.mktemp(suffix='.json')
    try:
        save_model(m, path)
        reloaded = load_model(path)
    finally:
        if os.path.exists(path):
            os.remove(path)

    src = [n for n in reloaded.getNodes() if isinstance(n, Source)][0]
    cls = [c for c in reloaded.getClasses() if c.name == 'C'][0]
    batch = src.getArrivalBatch(cls)
    assert batch is not None, "batch law must survive the round trip"
    assert batch.getMean() == pytest.approx(2.0), "E[X] = 1/beta"
    # The interarrival law must be untouched by the batch round trip.
    assert src.getArrival(cls).getMean() == pytest.approx(10.0)


def test_model_without_batch_omits_the_key():
    # A key that is always emitted would change every existing model.json.
    m = Network('plain')
    src = Source(m, 'S')
    q = Queue(m, 'Q', SchedStrategy.FCFS)
    sk = Sink(m, 'K')
    c = OpenClass(m, 'C')
    src.setArrival(c, Geometric(0.1))
    q.setService(c, Geometric(0.9))
    P = m.initRoutingMatrix()
    P.set(c, c, src, q, 1.0)
    P.set(c, c, q, sk, 1.0)
    m.link(P)

    path = tempfile.mktemp(suffix='.json')
    try:
        save_model(m, path)
        with open(path) as fh:
            doc = json.load(fh)
    finally:
        if os.path.exists(path):
            os.remove(path)
    assert 'arrivalBatch' not in _find_source_json(doc)


# ----------------------------------------------------------------------
# Batch-law validation
# ----------------------------------------------------------------------

def test_batch_law_supported_below_one_is_rejected():
    from line_solver.distributions.discrete import Bernoulli
    m, src, c = _geo_x_model()
    # Bernoulli(0.5) has mean 0.5, so some "batches" would carry no job.
    with pytest.raises(ValueError):
        src.setArrivalBatch(c, Bernoulli(0.5))


def test_batch_law_must_be_discrete():
    from line_solver.distributions.continuous import Exp
    m, src, c = _geo_x_model()
    with pytest.raises(ValueError):
        src.setArrivalBatch(c, Exp(0.5))


def test_none_restores_single_arrivals():
    m, src, c = _geo_x_model()
    assert src.getArrivalBatch(c) is not None
    src.setArrivalBatch(c, None)
    assert src.getArrivalBatch(c) is None


# ----------------------------------------------------------------------
# Slotted CLI flags
# ----------------------------------------------------------------------

def _cli_args(**kwargs):
    from line_solver.solvers.wrappers.solver_ldes.solver_ldes import SolverLDES
    m, _, _ = _geo_x_model()
    opts = LDESOptions(**kwargs)
    solver = SolverLDES(m, opts)
    return solver._build_cli_args('model.json', 'result.json')


def test_slotted_flag_absent_by_default():
    args = _cli_args()
    assert '--slotted' not in args
    assert '--slotlength' not in args


def test_slotted_flag_emitted():
    args = _cli_args(slotted=True)
    assert '--slotted' in args
    # Unit slot length is the CLI default, so it is not restated.
    assert '--slotlength' not in args


def test_slot_length_emitted_when_non_default():
    args = _cli_args(slotted=True, slot_length=0.5)
    assert '--slotted' in args
    assert '--slotlength' in args
    assert args[args.index('--slotlength') + 1] == repr(0.5)


def test_slotted_survives_options_copy():
    # LDESOptions.copy() enumerates every field by hand, so an omitted field is
    # silently dropped rather than failing loudly.
    opts = LDESOptions(slotted=True, slot_length=0.25)
    clone = opts.copy()
    assert clone.slotted is True
    assert clone.slot_length == 0.25
