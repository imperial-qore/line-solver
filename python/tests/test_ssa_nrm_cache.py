"""
Native-Python tests for NRM cache support (core, no retrieval system).

Mirror of the JAR ``SolverSSANrmCacheTest``. A Cache node is simulated by the
NRM as an immediate state-dependent class switch: a read-class job draws an item
from ``pread``, the cache contents decide a hit (``hitClass``) or a miss
(``missClass``), and the replacement policy updates the contents -- a faithful
port of ``State.afterEventCache`` (READ, isSimulation). Retrieval/delayed-hit is
NOT covered here (it routes to the serial engine).

Oracle: the same-codebase CTMC. Fixture: closed 1-class
Delay(Exp 1.0) -> Cache -> Delay with read/hit/miss classes. The cache hit ratio
is the Cache HitClass throughput / (Hit + Miss) from the node table.

The BINDING control is the Zipf policy divergence: under a skewed popularity the
recency policies (LRU/HLRU) keep the hot items and beat FIFO/SFIFO/RR, so a
no-op or serial-fallback implementation (which returns a policy-independent
ratio) fails ``test_zipf_policy_divergence``. Uniform access is the degenerate
control where every policy holds m of n items and the ratio is m/n for all.
"""

import numpy as np
import pandas as pd
import pytest

from line_solver import (
    Network, Delay, Cache, ClosedClass, Exp, DiscreteSampler,
    ReplacementStrategy, SolverSSA, SolverCTMC,
)

SAMPLES = 400000
SEED = 24000
REL_TOL = 0.03


def _build(n, cap, strat, pvec):
    model = Network('cache')
    d = Delay(model, 'Delay')
    cache = Cache(model, 'Cache', n, cap, strat)
    job = ClosedClass(model, 'JobClass', 1, d, 0)
    hit = ClosedClass(model, 'HitClass', 0, d, 0)
    miss = ClosedClass(model, 'MissClass', 0, d, 0)
    d.set_service(job, Exp(1))
    cache.set_read(job, DiscreteSampler(np.array(pvec)))
    cache.set_hit_class(job, hit)
    cache.set_miss_class(job, miss)
    P = model.init_routing_matrix()
    P.set(job, job, d, cache, 1.0)
    P.set(hit, job, cache, d, 1.0)
    P.set(miss, job, cache, d, 1.0)
    model.link(P)
    return model


def _hit_ratio(solver):
    t = solver.get_avg_node_table()
    df = t if isinstance(t, pd.DataFrame) else t.tabulate()
    tcol = [c for c in df.columns if 'Tput' in c][0]
    ncol = [c for c in df.columns if c.lower() in ('node', 'station')][0]
    ccol = [c for c in df.columns if 'class' in c.lower() or c in ('JobClass', 'Class')][0]
    h = m = 0.0
    for _, row in df.iterrows():
        if str(row[ncol]) == 'Cache' and str(row[ccol]) == 'HitClass':
            h = float(row[tcol])
        if str(row[ncol]) == 'Cache' and str(row[ccol]) == 'MissClass':
            m = float(row[tcol])
    assert (h + m) > 0, 'no cache hit/miss throughput measured'
    return h / (h + m)


def _nrm_ratio(model):
    s = SolverSSA(model, 'nrm', samples=SAMPLES, seed=SEED, verbose=False)
    r = _hit_ratio(s)
    # the dispatch must actually use NRM, not silently fall back to serial
    assert 'nrm' in str(getattr(s, 'method', 'nrm')).lower()
    return r


def _zipf(n, a=1.2):
    w = (np.arange(1, n + 1)) ** (-a)
    return list(w / w.sum())


POLICIES = [
    ('LRU', ReplacementStrategy.LRU), ('FIFO', ReplacementStrategy.FIFO),
    ('RR', ReplacementStrategy.RR), ('SFIFO', ReplacementStrategy.SFIFO),
    ('HLRU', ReplacementStrategy.HLRU), ('QLRU', ReplacementStrategy.QLRU),
]


@pytest.mark.parametrize('name,strat', POLICIES, ids=[p[0] for p in POLICIES])
def test_uniform_matches_ctmc(name, strat):
    # uniform access: every policy holds m of n items, ratio = m/n = 0.4
    pv = [0.2] * 5
    r_nrm = _nrm_ratio(_build(5, 2, strat, pv))
    r_ctmc = _hit_ratio(SolverCTMC(_build(5, 2, strat, pv), cutoff=20))
    assert r_nrm == pytest.approx(r_ctmc, rel=REL_TOL), f'{name} uniform: {r_nrm} vs {r_ctmc}'


@pytest.mark.parametrize('name,strat', POLICIES, ids=[p[0] for p in POLICIES])
def test_zipf_matches_ctmc(name, strat):
    pv = _zipf(6)
    r_nrm = _nrm_ratio(_build(6, 2, strat, pv))
    r_ctmc = _hit_ratio(SolverCTMC(_build(6, 2, strat, pv), cutoff=20))
    assert r_nrm == pytest.approx(r_ctmc, rel=REL_TOL), f'{name} zipf: {r_nrm} vs {r_ctmc}'


def test_zipf_policy_divergence():
    # binding control: under Zipf, LRU (recency) must beat FIFO; an exp-collapsed
    # or serial-fallback impl would return a policy-independent ratio.
    pv = _zipf(6)
    r_lru = _nrm_ratio(_build(6, 2, ReplacementStrategy.LRU, pv))
    r_fifo = _nrm_ratio(_build(6, 2, ReplacementStrategy.FIFO, pv))
    assert r_lru > r_fifo + 0.01, f'LRU {r_lru} should exceed FIFO {r_fifo} under Zipf'


def test_multilist_matches_ctmc():
    pv = _zipf(6, a=1.0)
    for name, strat in [('HLRU', ReplacementStrategy.HLRU), ('FIFO', ReplacementStrategy.FIFO)]:
        r_nrm = _nrm_ratio(_build(6, [2, 2], strat, pv))
        r_ctmc = _hit_ratio(SolverCTMC(_build(6, [2, 2], strat, pv), cutoff=20))
        assert r_nrm == pytest.approx(r_ctmc, rel=REL_TOL), f'{name} multilist: {r_nrm} vs {r_ctmc}'
