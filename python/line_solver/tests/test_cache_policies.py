"""
Regression tests for list-based cache replacement policies.

Validates the CTMC (exact) hit ratios of HLRU (h-LRU / LRU(m)), CLIMB
(transposition) and QLRU (q-LRU probabilistic admission) against MATLAB
SolverCTMC reference values (2026-07-12), and the MVA characteristic-time
(TTL) approximation for HLRU (Gast and Van Houdt, SIGMETRICS 2015).
"""

import unittest
import os
import sys

import numpy as np

# Add the python root for direct execution
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from line_solver import (
    Network, Source, Sink, Cache, OpenClass, Exp, Zipf,
    ReplacementStrategy, GlobalConstants, VerboseLevel, CTMC, MVA,
    MarkedMAP, DiscreteSampler,
)


def _build_open_cache(n, caps, strategy, alpha, q=None):
    model = Network('cache_policy_test')
    source = Source(model, 'Source')
    cache_node = Cache(model, 'Cache', n, caps, strategy)
    if q is not None:
        cache_node.set_admission_prob(q)
    sink = Sink(model, 'Sink')
    job_class = OpenClass(model, 'InitClass', 0)
    hit_class = OpenClass(model, 'HitClass', 0)
    miss_class = OpenClass(model, 'MissClass', 0)
    source.set_arrival(job_class, Exp(1))
    cache_node.set_read(job_class, Zipf(alpha, n))
    cache_node.set_hit_class(job_class, hit_class)
    cache_node.set_miss_class(job_class, miss_class)
    P = model.init_routing_matrix()
    P.set(job_class, job_class, source, cache_node, 1.0)
    P.set(hit_class, hit_class, cache_node, sink, 1.0)
    P.set(miss_class, miss_class, cache_node, sink, 1.0)
    model.link(P)
    return model


def _hit_ratio(table):
    cols = list(table.columns)
    h = float(table[(table[cols[0]] == 'Cache')
                    & (table[cols[1]] == 'HitClass')]['Tput'].iloc[0])
    m = float(table[(table[cols[0]] == 'Cache')
                    & (table[cols[1]] == 'MissClass')]['Tput'].iloc[0])
    return h / (h + m)


class TestCachePolicies(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        GlobalConstants.set_verbose(VerboseLevel.SILENT)

    def test_hlru_ctmc_exact(self):
        # MATLAB SolverCTMC reference (2026-07-12)
        model = _build_open_cache(6, [2, 1], ReplacementStrategy.HLRU, 1.2)
        hr = _hit_ratio(CTMC(model, cutoff=1).avg_node_table())
        self.assertAlmostEqual(hr, 0.700849020321130, places=6)

    def test_hlru_mva_ttl(self):
        # MATLAB SolverMVA TTL (LRU(m) characteristic-time) reference
        model = _build_open_cache(6, [2, 1], ReplacementStrategy.HLRU, 1.2)
        hr = _hit_ratio(MVA(model).avg_node_table())
        self.assertAlmostEqual(hr, 0.701566, places=4)

    def test_hlru_h1_reduces_to_lru(self):
        # single-list h-LRU must coincide with the LRU (Che) approximation
        m_hlru = _build_open_cache(6, 2, ReplacementStrategy.HLRU, 1.2)
        m_lru = _build_open_cache(6, 2, ReplacementStrategy.LRU, 1.2)
        hr_hlru = _hit_ratio(MVA(m_hlru).avg_node_table())
        hr_lru = _hit_ratio(MVA(m_lru).avg_node_table())
        self.assertAlmostEqual(hr_hlru, hr_lru, places=6)

    def test_climb_ctmc_exact(self):
        # MATLAB SolverCTMC reference (2026-07-12)
        model = _build_open_cache(5, 2, ReplacementStrategy.CLIMB, 1.2)
        hr = _hit_ratio(CTMC(model, cutoff=1).avg_node_table())
        self.assertAlmostEqual(hr, 0.597633165363664, places=6)

    def test_qlru_ctmc_exact(self):
        # MATLAB SolverCTMC reference (2026-07-12), admission prob q=0.5
        model = _build_open_cache(5, 2, ReplacementStrategy.QLRU, 1.2, q=0.5)
        hr = _hit_ratio(CTMC(model, cutoff=1).avg_node_table())
        self.assertAlmostEqual(hr, 0.575143053074336, places=6)

    def test_lru_marked_map_mva(self):
        # Marked MMAP source driving per-item classes of an LRU [2,1] cache
        # (phase 1 prefers items 1-2, phase 2 items 3-4, non-IRM). MVA
        # dispatches to cache_ttl_lrum_map; MATLAB/JAR reference 0.8400817539
        # (exact CTMC 0.8108225223), 2026-07-12.
        n, m = 4, [2, 1]
        D0 = np.array([[-2.1, 0.1], [0.2, -0.7]])
        D1 = np.array([[2.0, 0.0], [0.0, 0.5]])
        pfast = [0.4, 0.3, 0.2, 0.1]
        pslow = [0.1, 0.2, 0.3, 0.4]
        Dk = [D1 @ np.diag([pfast[k], pslow[k]]) for k in range(n)]
        model = Network('m')
        source = Source(model, 'Source')
        cache_node = Cache(model, 'Cache', n, m, ReplacementStrategy.LRU)
        sink = Sink(model, 'Sink')
        rc = [OpenClass(model, 'Item%d' % (k + 1), 0) for k in range(n)]
        hit_class = OpenClass(model, 'HitClass', 0)
        miss_class = OpenClass(model, 'MissClass', 0)
        source.set_marked_arrival(MarkedMAP([D0] + Dk), rc)
        P = model.init_routing_matrix()
        for k in range(n):
            onehot = np.zeros(n)
            onehot[k] = 1.0
            cache_node.set_read(rc[k], DiscreteSampler(onehot, np.arange(1, n + 1)))
            cache_node.set_hit_class(rc[k], hit_class)
            cache_node.set_miss_class(rc[k], miss_class)
            P.set(rc[k], rc[k], source, cache_node, 1.0)
        P.set(hit_class, hit_class, cache_node, sink, 1.0)
        P.set(miss_class, miss_class, cache_node, sink, 1.0)
        model.link(P)
        hr = _hit_ratio(MVA(model).avg_node_table())
        self.assertAlmostEqual(hr, 0.8400817539, places=6)


if __name__ == '__main__':
    unittest.main()
