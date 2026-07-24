"""Regression: a Cache node's restored initial state must not suppress the
default placement of the closed-class population.

The Cache state vector is [class counts | contents | retrieval bitmap], not a
per-class marginal. _refresh_state used to read its first nclasses entries as an
explicit marginal and mark every class as already placed, so a closed network
whose model.json carried a cache initialState (which linemodel_save writes once
getStruct has run) was solved with zero jobs and returned an empty average
table. MATLAB and the JAR never had this: MNetwork.getState / Network.getState
call initDefault whenever hasInitState is false, and a state on a strict subset
of the stateful nodes leaves hasInitState false.
"""
import numpy as np

from line_solver import (Network, Delay, Cache, ClosedClass, ReplacementStrategy,
                         Zipf, Exp, Immediate, SolverSSA)


def _build():
    nitems, cachesize = 1000, 50
    model = Network('Model')
    client = Delay(model, 'Client')
    cache = Cache(model, 'Cache', nitems, cachesize, ReplacementStrategy.LRU)
    cachedelay = Delay(model, 'CacheDelay')

    cclass = ClosedClass(model, 'ClientClass', 1, client, 0)
    hclass = ClosedClass(model, 'HitClass', 0, client, 0)
    mclass = ClosedClass(model, 'MissClass', 0, client, 0)

    client.setService(cclass, Immediate())
    cachedelay.setService(hclass, Exp.fitMean(0.2))
    cachedelay.setService(mclass, Exp.fitMean(1.0))

    cache.setRead(cclass, Zipf(1.4, nitems))
    cache.setHitClass(cclass, hclass)
    cache.setMissClass(cclass, mclass)

    P = model.initRoutingMatrix()
    P.set(cclass, cclass, client, cache, 1.0)
    P.set(hclass, hclass, cache, cachedelay, 1.0)
    P.set(mclass, mclass, cache, cachedelay, 1.0)
    P.set(hclass, cclass, cachedelay, client, 1.0)
    P.set(mclass, cclass, cachedelay, client, 1.0)
    model.link(P)
    return model, cache


def test_cache_initstate_preserves_closed_population():
    model, cache = _build()
    ref = SolverSSA(model, samples=20000, seed=23000).getAvg()[0]

    # Reproduce what load_model does for a model.json exported after getStruct:
    # only the Cache node carries a state, every other stateful node does not.
    model2, cache2 = _build()
    nclasses = 3
    cachesize = 50
    state = np.concatenate([np.zeros(nclasses), np.arange(1, cachesize + 1)])
    cache2.setState(state)
    got = SolverSSA(model2, samples=20000, seed=23000).getAvg()[0]

    assert np.asarray(ref).sum() > 0
    np.testing.assert_allclose(np.asarray(got), np.asarray(ref), rtol=1e-12)


def test_partial_state_does_not_zero_the_population():
    model, _ = _build()
    sn = model.getStruct()
    assert float(np.sum([np.sum(np.asarray(s)) for s in sn.state])) > 0

    model2, cache2 = _build()
    cache2.setState(np.concatenate([np.zeros(3), np.arange(1, 51)]))
    sn2 = model2.getStruct()
    assert float(np.sum([np.sum(np.asarray(s)) for s in sn2.state])) > 0
