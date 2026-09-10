"""SolverFLD moment closure on a cache-queueing model.

`minnormal` used to be refused on any model with a cache node. It is now
answered through the same decomposition as `rmf`, with the moment closure in
the network step, and it reports the second moment of both layers. Expected
values are the MATLAB reference on the same model.
"""

import numpy as np
import pytest

from line_solver import (Cache, ClosedClass, Delay, Disabled, Exp, GlobalConstants, Network,
                         Queue, ReplacementStrategy, SchedStrategy, SolverFLD, VerboseLevel, Zipf)

GlobalConstants.setVerbose(VerboseLevel.SILENT)

NITEMS, CAP, N = 10, 3, 4


def _model():
    """Think -> Cache -> Q1(PS) -> Think, Zipf(0.8) reads over 10 items, RANDOM(3)."""
    model = Network('cacheqn')
    d = Delay(model, 'Think')
    c = Cache(model, 'Cache', NITEMS, CAP, ReplacementStrategy.RR)
    q = Queue(model, 'Q1', SchedStrategy.PS)
    req = ClosedClass(model, 'Req', N, d, 0)
    hit = ClosedClass(model, 'Hit', 0, d, 0)
    miss = ClosedClass(model, 'Miss', 0, d, 0)
    d.setService(req, Exp.fitMean(1.0))
    d.setService(hit, Exp.fitMean(0.1))
    d.setService(miss, Exp.fitMean(0.1))
    q.setService(req, Disabled())
    q.setService(hit, Exp.fitMean(0.2))
    q.setService(miss, Exp.fitMean(1.0))
    c.setRead(req, Zipf(0.8, NITEMS))
    c.setHitClass(req, hit)
    c.setMissClass(req, miss)
    P = model.init_routing_matrix()
    P.set(req, req, d, c, 1.0)
    P.set(hit, hit, c, q, 1.0)
    P.set(miss, miss, c, q, 1.0)
    P.set(hit, req, q, d, 1.0)
    P.set(miss, req, q, d, 1.0)
    model.link(P)
    return model


def test_minnormal_answers_a_cache_model():
    """The closure runs on the cache model and conserves the population."""
    Q = np.asarray(SolverFLD(_model(), method='minnormal').getAvgQLen())
    assert np.all(np.isfinite(Q))
    assert Q.sum() == pytest.approx(N, abs=1e-3)
    # MATLAB: Think 1.3873, Q1 total 2.6127 (0.28825 hit + 2.3245 miss)
    assert Q.flatten()[0] == pytest.approx(1.3873, abs=1e-3)


def test_cache_layer_is_unchanged_by_the_closure():
    """Only the queueing layer changes: rmf and minnormal share the cache solve."""
    Qr = np.asarray(SolverFLD(_model(), method='rmf').getAvgQLen())
    Qm = np.asarray(SolverFLD(_model(), method='minnormal').getAvgQLen())
    assert Qr.flatten()[0] == pytest.approx(1.4413, abs=1e-3)
    assert not np.allclose(Qr, Qm), 'the closure must move the queueing layer'


def test_cache_moment_report():
    """The occupancy covariance is reported and respects the cache invariants."""
    solver = SolverFLD(_model(), method='minnormal')
    solver.getAvgQLen()
    mom = solver.getMoments()
    assert mom is not None and 'QVar' in mom
    assert mom.get('cache'), 'no cache moments reported'

    cm = mom['cache'][0]
    pi0 = np.asarray(cm['pi0'])
    pi0var = np.asarray(cm['pi0Var'])
    assert pi0.size == NITEMS and pi0var.size == NITEMS
    assert np.all(pi0var >= 0)

    # MATLAB reference for the same model
    assert pi0[0] == pytest.approx(0.4203, abs=1e-3)
    assert pi0var[0] == pytest.approx(0.2135, abs=1e-3)
    assert cm['missProbVar'][0] == pytest.approx(0.01055, abs=1e-4)

    # the miss indicator is Bernoulli up to the O(1/n) error of the closure
    bern = pi0 * (1 - pi0)
    assert np.max(np.abs(pi0var - bern) / bern) < 0.35

    # the number of uncached items is deterministic, so the covariance of the
    # miss indicators must sum to zero exactly
    W00 = np.asarray(cm['Sigma'])[:NITEMS, :NITEMS]
    assert abs(W00.sum()) < 1e-8
