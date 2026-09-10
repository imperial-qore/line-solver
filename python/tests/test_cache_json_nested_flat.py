"""A cache model.json may split its keys between `cache` and the node itself.

MATLAB's `linemodel_save` writes hitClass / missClass / popularity NESTED under
a `cache` object and retrievalSystem FLAT on the node. The reader used to pick
one source per node (`cache_src = cc if cc else nd`), so the flat keys were
dropped whenever the nested object existed: a MATLAB-exported delayed-hit model
loaded here as a PLAIN cache, hit + miss summed to 1, and both `lang='python'`
and `lang='cpp'` reported it without complaint.

The read is per KEY now. These tests pin that, in both shapes.
"""
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from line_solver import (Cache, DiscreteSampler, Exp, Network, OpenClass,
                         Queue, ReplacementStrategy, SchedStrategy, Sink,
                         Source, MVA, load_model, save_model)

NESTED_KEYS = ('hitClass', 'missClass', 'popularity')


def _retrieval_model():
    access_prob = [0.6, 0.3, 0.1]
    model = Network('DelayedHits')
    source = Source(model, 'Source')
    cache_node = Cache(model, 'Cache', len(access_prob), [1],
                       ReplacementStrategy.FIFO)
    queue = Queue(model, 'Queue', SchedStrategy.INF)
    sink = Sink(model, 'Sink')
    job_class = OpenClass(model, 'InitClass', 0)
    hit_class = OpenClass(model, 'HitClass', 0)
    miss_class = OpenClass(model, 'MissClass', 0)
    source.set_arrival(job_class, Exp(1))
    queue.set_service(job_class, Exp(2.0))
    cache_node.set_read(job_class, DiscreteSampler(access_prob))
    cache_node.set_hit_class(job_class, hit_class)
    cache_node.set_miss_class(job_class, miss_class)
    cache_node.set_retrieval_system(job_class, miss_class, queue)
    P = model.init_routing_matrix()
    P.set(job_class, job_class, source, cache_node, 1.0)
    P.set(job_class, job_class, cache_node, queue, 1.0)
    P.set(job_class, job_class, queue, cache_node, 1.0)
    P.set(hit_class, hit_class, cache_node, sink, 1.0)
    P.set(miss_class, miss_class, cache_node, sink, 1.0)
    model.link(P)
    return model


def _to_matlab_shape(doc):
    """Move the nested-in-MATLAB keys under `cache`, leave retrievalSystem flat."""
    for nd in doc['model']['nodes']:
        if nd.get('type') != 'Cache' and 'retrievalSystem' not in nd:
            continue
        cache = nd.setdefault('cache', {})
        for key in NESTED_KEYS:
            if key in nd:
                cache[key] = nd[key]
    return doc


def _cache_row(model):
    return MVA(model).get_avg_cache_table().iloc[0]


@pytest.fixture(scope='module')
def exported(tmp_path_factory):
    path = str(tmp_path_factory.mktemp('cachejson') / 'model.json')
    save_model(_retrieval_model(), path)
    with open(path) as f:
        return json.load(f)


def test_flat_shape_keeps_the_retrieval_system(exported, tmp_path):
    """The shape python writes itself: every cache key flat on the node."""
    path = str(tmp_path / 'flat.json')
    with open(path, 'w') as f:
        json.dump(exported, f)
    node = [n for n in load_model(path).get_nodes() if isinstance(n, Cache)][0]
    assert node._retrieval_system_capacity > 0


def test_matlab_shape_keeps_the_retrieval_system(exported, tmp_path):
    """MATLAB's shape: SOME keys nested under `cache`, retrievalSystem flat.

    This is the regression. With a per-node source choice the nested object
    wins and the flat retrievalSystem is silently lost.
    """
    path = str(tmp_path / 'matlab.json')
    with open(path, 'w') as f:
        json.dump(_to_matlab_shape(json.loads(json.dumps(exported))), f)
    node = [n for n in load_model(path).get_nodes() if isinstance(n, Cache)][0]
    assert node._retrieval_system_capacity > 0, \
        'retrievalSystem dropped when a nested `cache` object is present'


def test_both_shapes_give_the_same_solved_cache_table(exported, tmp_path):
    """The shape of the file must not change the answer.

    A dropped retrieval system does not raise; it reports a plain cache whose
    hit and miss sum to 1 and whose delayed-hit column is zero. Asserting the
    3-way split is what makes the loss visible.
    """
    flat_path = str(tmp_path / 'flat2.json')
    matlab_path = str(tmp_path / 'matlab2.json')
    with open(flat_path, 'w') as f:
        json.dump(exported, f)
    with open(matlab_path, 'w') as f:
        json.dump(_to_matlab_shape(json.loads(json.dumps(exported))), f)

    flat = _cache_row(load_model(flat_path))
    matlab = _cache_row(load_model(matlab_path))

    for col in ('HitProb', 'DelayedHitProb', 'MissProb', 'ArvR'):
        assert flat[col] == pytest.approx(matlab[col], rel=1e-12), col
    assert matlab['DelayedHitProb'] > 0, 'retrieval system produced no delayed hits'
    assert (matlab['HitProb'] + matlab['DelayedHitProb']
            + matlab['MissProb']) == pytest.approx(1.0, abs=1e-9)
