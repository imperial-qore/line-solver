"""
Native-Python tests for the exact delayed-hit queue length in SolverCTMC.

Mirror of the JAR ``SolverCTMCDelayedHitTest`` and of the MATLAB
``test_ctmc_delayed_hit``. Block A of the cache local state marks the items being
fetched and block B counts the secondary requests merged onto each fetch, so per
item i, phi_i = P(block A bit i set), d1_i = E[block B slots of i] and
dfull_i = d1_i + phi_i are state rewards of the stationary distribution. The
delayed-hit RATE is a transition reward over the generator (the transitions that
clear block A bit i); it splits the hit-class rate into true hits and delayed hits.

Fixture: open, 3 items, cache capacity 1, FIFO, Exp(1) arrivals, Exp(2) fetch.
The cutoff truncates block B at maxPending = cutoff - 1, so every finite cutoff
gives a LOWER bound on the delayed-hit mass and approaches SolverNC (exact for
this product-form retrieval cache) from below.

The BINDING control is ``test_delayed_hit_mass_increases_with_cutoff``: an
implementation that drops merged requests instead of parking them in block B
reports zero delayed hits at every cutoff and fails it. The golden values are
shared verbatim with the JAR and MATLAB twins, so they also pin cross-codebase
parity.
"""

import numpy as np
import pytest

from line_solver import (Network, Source, Cache, Queue, Sink, OpenClass, Exp,
                         DiscreteSampler, ReplacementStrategy, SchedStrategy,
                         CTMC, NC)

TOL = 1e-4


def _model():
    model = Network('DelayedHits')
    n, capacity = 3, [1]
    source = Source(model, 'Source')
    cache_node = Cache(model, 'Cache', n, capacity, ReplacementStrategy.FIFO)
    queue = Queue(model, 'Queue', SchedStrategy.INF)
    sink = Sink(model, 'Sink')
    job_class = OpenClass(model, 'InitClass', 0)
    hit_class = OpenClass(model, 'HitClass', 0)
    miss_class = OpenClass(model, 'MissClass', 0)
    source.set_arrival(job_class, Exp(1))
    queue.set_service(job_class, Exp(2.0))
    cache_node.set_read(job_class, DiscreteSampler([0.6, 0.3, 0.1]))
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
    return model, cache_node


def _solve(cutoff):
    model, cache = _model()
    solver = CTMC(model, cutoff=cutoff)
    solver._table_silent = True
    table = solver.get_avg_cache_table()
    d1, dfull = cache.get_delayed_hit_qlen()
    return table, np.ravel(d1), np.ravel(dfull)


@pytest.fixture(scope='module')
def cutoff2():
    return _solve(2)


@pytest.fixture(scope='module')
def cutoff3():
    return _solve(3)


def test_delayed_hit_qlen_cutoff2(cutoff2):
    _, d1, dfull = cutoff2
    assert d1.size == 3
    np.testing.assert_allclose(d1, [0.037531, 0.020324, 0.003728], atol=TOL)
    np.testing.assert_allclose(dfull, [0.138703, 0.110729, 0.047859], atol=TOL)


def test_dfull_minus_d1_is_fetch_probability(cutoff2):
    # dfull_i - d1_i = phi_i = P(a fetch of item i is in flight), by construction:
    # dfull counts the triggering request, d1 does not.
    _, d1, dfull = cutoff2
    phi = dfull - d1
    assert np.all(phi > 0) and np.all(phi < 1)


def test_hit_delayed_miss_split_sums_to_one(cutoff2):
    # Delayed hits depart in the hit class; the exact transition reward splits the
    # hit-class rate so the three reported fractions partition the arrivals.
    table, _, _ = cutoff2
    h = float(table.HitProb[0])
    d = float(table.DelayedHitProb[0])
    m = float(table.MissProb[0])
    assert h == pytest.approx(0.467000, abs=TOL)
    assert d == pytest.approx(0.061583, abs=TOL)
    assert m == pytest.approx(0.471417, abs=TOL)
    assert h + d + m == pytest.approx(1.0, abs=1e-9)


def test_delayed_hit_mass_increases_with_cutoff(cutoff2, cutoff3):
    # maxPending = cutoff - 1 truncates block B, so a coarser cutoff can only
    # UNDERSTATE the delayed-hit mass and the delayed-hit queue length. An
    # implementation that absorbs merged requests reports 0 at every cutoff.
    t2, d1_2, _ = cutoff2
    t3, d1_3, _ = cutoff3
    assert float(t3.DelayedHitProb[0]) > float(t2.DelayedHitProb[0])
    assert d1_3[0] > d1_2[0]


def test_converges_toward_exact_normalizing_constant_solver(cutoff2, cutoff3):
    # SolverNC solves this retrieval cache exactly (product form), so the truncated
    # CTMC must approach it from below as the cutoff grows.
    model, _ = _model()
    nc = NC(model)
    nc._table_silent = True
    nct = nc.get_avg_cache_table()
    exact = float(nct.HitProb[0]) + float(nct.DelayedHitProb[0])

    t2, _, _ = cutoff2
    t3, _, _ = cutoff3
    s2 = float(t2.HitProb[0]) + float(t2.DelayedHitProb[0])
    s3 = float(t3.HitProb[0]) + float(t3.DelayedHitProb[0])
    assert s2 < s3
    assert s3 < exact + TOL
    assert (exact - s3) < (exact - s2)
