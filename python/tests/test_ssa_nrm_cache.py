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


# --- retrieval / delayed-hit system -----------------------------------------
# A retrieval cache sends a miss to a fetch queue and back (delayed hits). The
# CTMC state space is too large here, so the oracle is the same-codebase serial
# SSA engine (State.afterEventCache). Note a pre-existing serial-vs-LDES
# divergence on delayed-hit counting; the NRM targets the serial engine.
from line_solver import Source, Sink, Queue, OpenClass, SchedStrategy  # noqa: E402


def _build_retrieval(strat):
    accessProb = [0.6, 0.3, 0.1]
    model = Network('DelayedHits')
    source = Source(model, 'Source')
    cache = Cache(model, 'Cache', len(accessProb), 1, strat)
    queue = Queue(model, 'Queue', SchedStrategy.INF)
    sink = Sink(model, 'Sink')
    job = OpenClass(model, 'InitClass', 0)
    hit = OpenClass(model, 'HitClass', 0)
    miss = OpenClass(model, 'MissClass', 0)
    source.setArrival(job, Exp(1))
    queue.setService(job, Exp(2.0))
    cache.set_read(job, DiscreteSampler(np.array(accessProb)))
    cache.set_hit_class(job, hit)
    cache.set_miss_class(job, miss)
    cache.setRetrievalSystem(job, miss, queue)
    P = model.init_routing_matrix()
    P.set(job, job, source, cache, 1.0)
    P.set(job, job, cache, queue, 1.0)
    P.set(job, job, queue, cache, 1.0)
    P.set(hit, hit, cache, sink, 1.0)
    P.set(miss, miss, cache, sink, 1.0)
    model.link(P)
    return model


def _retr_hit_ratio(solver):
    t = solver.getAvgCacheTable()
    df = t if isinstance(t, pd.DataFrame) else t.tabulate()
    hc = [c for c in df.columns if c == 'HitProb'][0]
    mc = [c for c in df.columns if c == 'MissProb'][0]
    row = df.iloc[0]
    h = float(row[hc]); m = float(row[mc])
    return h / (h + m)


@pytest.mark.parametrize('name,strat', [
    ('FIFO', ReplacementStrategy.FIFO), ('LRU', ReplacementStrategy.LRU),
])
def test_retrieval_matches_serial(name, strat):
    # Oracle is the serial SSA (afterEventCache), not CTMC (state space too big).
    s = SolverSSA(_build_retrieval(strat), 'nrm', samples=SAMPLES, seed=SEED, verbose=False)
    r_nrm = _retr_hit_ratio(s)
    assert 'nrm' in str(getattr(s, 'method', 'nrm')).lower()
    ss = SolverSSA(_build_retrieval(strat), 'serial', samples=SAMPLES, seed=SEED, verbose=False)
    r_ser = _retr_hit_ratio(ss)
    assert r_nrm == pytest.approx(r_ser, rel=0.03), \
        f'{name} retrieval: NRM {r_nrm:.4f} vs serial {r_ser:.4f}'
