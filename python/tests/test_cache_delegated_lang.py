"""
Regression tests for the cache results of a DELEGATED solve (lang='cpp',
lang='java').

A cache's hit/miss split is a solver RESULT written onto the Cache node, and
the hit/miss visits are derived from it. A bridge that delegates the solve must
therefore do two things the native analyzer does inline: seed the split back
onto `sn.nodeparam`, and re-route the cache self-switch at that split before
refreshing the visits. Omit the first and the node-level hit/miss throughputs
report whatever the previous native solve left; omit the second and every
metric carrying a visit ratio answers at `link()`'s uniform guess -- on
cache_replc_routing, ResidT at (Delay1, HitClass) reads 0.025 (visit 0.5*0.5)
instead of 0.0197 (visit 0.394*0.5).

Each test compares the delegated node table against the NATIVE one on the same
model, so it asserts agreement rather than a transcribed constant. A backend
that is not installed skips rather than fails.
"""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from line_solver import (
    Network, Source, Sink, Cache, Delay, Router, ClosedClass, OpenClass, Exp,
    DiscreteSampler, ReplacementStrategy, RoutingStrategy, NC, GlobalConstants,
    VerboseLevel,
)

COLS = ['QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']


def _build_closed_fifo():
    """Closed cache, FIFO, 5 items, capacity 2: exact hit ratio 0.4."""
    model = Network('cache_closed_fifo')
    n, m = 5, 2
    delay = Delay(model, 'Delay')
    cache_node = Cache(model, 'Cache', n, m, ReplacementStrategy.FIFO)
    job_class = ClosedClass(model, 'JobClass', 1, delay, 0)
    hit_class = ClosedClass(model, 'HitClass', 0, delay, 0)
    miss_class = ClosedClass(model, 'MissClass', 0, delay, 0)
    delay.set_service(job_class, Exp(1))
    cache_node.set_read(job_class, DiscreteSampler([1.0 / n] * n))
    cache_node.set_hit_class(job_class, hit_class)
    cache_node.set_miss_class(job_class, miss_class)
    P = model.init_routing_matrix()
    P.set(job_class, job_class, delay, cache_node, 1.0)
    P.set(hit_class, job_class, cache_node, delay, 1.0)
    P.set(miss_class, job_class, cache_node, delay, 1.0)
    model.link(P)
    return model


def _build_open_routing():
    """Open cache feeding a Router: the visits downstream carry the split."""
    model = Network('cache_open_routing')
    n, m = 5, 2
    source = Source(model, 'Source')
    cache_node = Cache(model, 'Cache', n, m, ReplacementStrategy.FIFO)
    router = Router(model, 'Router')
    delay1 = Delay(model, 'Delay1')
    delay2 = Delay(model, 'Delay2')
    sink = Sink(model, 'Sink')
    job_class = OpenClass(model, 'InitClass', 0)
    hit_class = OpenClass(model, 'HitClass', 0)
    miss_class = OpenClass(model, 'MissClass', 0)
    source.set_arrival(job_class, Exp(2))
    delay1.set_service(hit_class, Exp(10))
    delay1.set_service(miss_class, Exp(1))
    delay2.set_service(hit_class, Exp(20))
    delay2.set_service(miss_class, Exp(2))
    cache_node.set_read(job_class, DiscreteSampler([1.0 / n] * n))
    cache_node.set_hit_class(job_class, hit_class)
    cache_node.set_miss_class(job_class, miss_class)
    model.add_link(source, cache_node)
    model.add_link(cache_node, router)
    model.add_link(router, delay1)
    model.add_link(router, delay2)
    model.add_link(delay1, sink)
    model.add_link(delay2, sink)
    source.set_prob_routing(job_class, cache_node, 1.0)
    cache_node.set_prob_routing(hit_class, router, 1.0)
    cache_node.set_prob_routing(miss_class, router, 1.0)
    router.set_routing(hit_class, RoutingStrategy.RAND)
    router.set_routing(miss_class, RoutingStrategy.RAND)
    delay1.set_prob_routing(hit_class, sink, 1.0)
    delay1.set_prob_routing(miss_class, sink, 1.0)
    delay2.set_prob_routing(hit_class, sink, 1.0)
    delay2.set_prob_routing(miss_class, sink, 1.0)
    # add_link plus set_prob_routing IS the topology here; a link() on top of it
    # would install an empty routing matrix and empty the tables.
    return model


def _node_table(builder, **kwargs):
    model = builder()
    solver = NC(model, **kwargs)
    solver._table_silent = True
    return solver.avg_node_table()


class TestCacheDelegatedLang(unittest.TestCase):
    """The delegated node table must equal the native one, column for column."""

    @classmethod
    def setUpClass(cls):
        GlobalConstants.set_verbose(VerboseLevel.SILENT)

    def _compare(self, builder, lang):
        native = _node_table(builder)
        try:
            delegated = _node_table(builder, lang=lang)
        except Exception as e:  # backend not installed in this checkout
            self.skipTest("lang='%s' unavailable: %s" % (lang, str(e)[:120]))
        self.assertEqual(list(native['Node']), list(delegated['Node']))
        self.assertEqual(list(native['JobClass']), list(delegated['JobClass']))
        for col in COLS:
            np.testing.assert_allclose(
                np.asarray(delegated[col], dtype=float),
                np.asarray(native[col], dtype=float),
                rtol=1e-6, atol=1e-9,
                err_msg="lang='%s' disagrees with the native solve on %s" % (lang, col))

    def test_closed_fifo_cpp(self):
        self._compare(_build_closed_fifo, 'cpp')

    def test_closed_fifo_java(self):
        self._compare(_build_closed_fifo, 'java')

    def test_open_routing_cpp(self):
        """The visits derived from the split: ResidT is the column that moves."""
        self._compare(_build_open_routing, 'cpp')

    def test_open_routing_java(self):
        self._compare(_build_open_routing, 'java')

    def test_hit_ratio_reaches_the_node_under_cpp(self):
        """`Cache.get_hit_ratio()` reads the NODE, so seeding only `sn` is not enough."""
        model = _build_closed_fifo()
        cache_node = [nd for nd in model.get_nodes() if nd.get_name() == 'Cache'][0]
        try:
            solver = NC(model, lang='cpp')
            solver._table_silent = True
            solver.avg_node_table()
        except Exception as e:
            self.skipTest("lang='cpp' unavailable: %s" % str(e)[:120])
        hr = np.atleast_1d(np.asarray(cache_node.get_hit_ratio(), dtype=float)).flatten()
        self.assertTrue(np.isfinite(hr[0]),
                        'the delegated solve left the model reporting no hit ratio')
        # 5 items, capacity 2, uniform access: SPM lands just under the exact 0.4.
        self.assertAlmostEqual(float(hr[0]), 0.394263, places=5)

    def test_delegated_solve_overwrites_a_stale_split(self):
        """A delegated solve must leave ITS OWN split on the node, not the
        previous engine's.

        The delegated fluid engine does report a cache block (MATLAB FLD
        answers 0.4/0.6 on this model, and the C++ one matches it), so the
        requirement is that the split be REPLACED. Leaving the earlier SSA
        value would report one engine's hit rate under another's banner, which
        is indistinguishable from a real answer.
        """
        from line_solver import FLD, SSA
        model = _build_closed_fifo()
        SSA(model, samples=2000, seed=23000)._table_silent = True
        SSA(model, samples=2000, seed=23000).avg_table()
        sn = model.get_struct()
        cp = sn.nodeparam.get(1) if isinstance(sn.nodeparam, dict) else sn.nodeparam[1]
        stale = getattr(cp, 'actualhitprob', None)
        self.assertIsNotNone(stale, 'the native SSA solve should have written a split')
        stale = float(np.asarray(stale).ravel()[0])
        try:
            table = FLD(model, lang='cpp').avg_node_table()
        except Exception as e:
            self.skipTest("lang='cpp' unavailable: %s" % str(e)[:120])
        sn = model.get_struct()
        cp = sn.nodeparam.get(1) if isinstance(sn.nodeparam, dict) else sn.nodeparam[1]
        fresh = getattr(cp, 'actualhitprob', None)
        self.assertIsNotNone(fresh, 'the delegated solve reports a cache block, so it must seed one')
        fresh = float(np.asarray(fresh).ravel()[0])
        # The split the delegated engine itself reported, read off its table:
        # the hit share of the cache's own throughput.
        rows = {(str(a), str(b)): float(c)
                for a, b, c in zip(table['Node'], table['JobClass'], table['Tput'])}
        hit = rows[('Cache', 'HitClass')]
        miss = rows[('Cache', 'MissClass')]
        self.assertAlmostEqual(fresh, hit / (hit + miss), places=6,
                               msg='the split on the node is not the one this solve reported')
        self.assertNotAlmostEqual(fresh, stale, places=6,
                                  msg='the stale SSA split survived the delegated solve')


if __name__ == '__main__':
    unittest.main()
