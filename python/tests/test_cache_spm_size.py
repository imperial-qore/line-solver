"""Ray (WKB) expansion of the cost-capped cache normalizing constant.

The profile is the one the derivation note validates against its own dynamic
program: gamma_i = 0.3 + 2.7 i/n, sizes {1,2,3} in equal thirds of the item
index (so gcd = 1 and the sizes are genuinely diverse), m = n/4, k = 1.8 m.
"""
import math
import warnings

import numpy as np
import pytest
from scipy.special import gammaln

from line_solver.api.cache import cache_erec, cache_spm_size


def profile(n):
    return (0.3 + 2.7 * (np.arange(1, n + 1) / n))[:, None]


def sizes(n):
    return 1 + np.arange(n) * 3 // n


@pytest.mark.parametrize("n,m,want", [(50, 12, 27.234323), (100, 25, 57.620291),
                                      (200, 50, 118.438233), (400, 100, 240.843731)])
def test_exact_cost_mode_reproduces_the_notes_table(n, m, want):
    k = math.ceil(1.8 * m)
    _, logE, _ = cache_spm_size(profile(n), [m], sizes(n), [k], 'exact')
    assert logE - gammaln(m + 1) == pytest.approx(want, abs=1e-5)


def test_cumulative_caps_track_cache_erec_and_degenerate_to_size_free():
    n, m = 120, 30
    g, sg = profile(n), sizes(n)

    # binding cap: the expansion sits close to the exact recursion
    _, logE, out = cache_spm_size(g, [m], sg, [54])
    exact = math.log(cache_erec(g, np.array([float(m)]), sg, np.array([54.0])))
    assert out.binding[0]
    assert out.zeta[0] < 1.0
    assert abs(logE - exact) < 0.15
    assert out.method == 'spm-size'
    assert out.iter < 30      # the saddle converges in a handful of Newton steps

    # slack cap: zeta returns to 1 and the cost coordinate leaves the saddle
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _, logE, out = cache_spm_size(g, [m], sg, [200])
    exact = math.log(cache_erec(g, np.array([float(m)]), sg, np.array([200.0])))
    assert not out.binding[0]
    assert out.zeta[0] == pytest.approx(1.0, abs=1e-12)
    assert abs(logE - exact) < 0.05
    assert out.iter < 30

    # a cap below the cheapest m items admits no state at all
    E, logE, out = cache_spm_size(g, [m], sg, [10])
    assert E == 0.0
    assert out.method == 'boundary'


def test_size_lattice_is_divided_out_and_uniform_sizes_are_exact():
    n, m = 120, 30
    g, sg = profile(n), sizes(n)
    _, base, _ = cache_spm_size(g, [m], sg, [54])
    # sizes 2/4/6 at cap 108 is the same problem; 109 snaps down to the same lattice point
    _, twice, out2 = cache_spm_size(g, [m], 2 * sg, [108])
    _, odd, _ = cache_spm_size(g, [m], 2 * sg, [109])
    assert out2.span == 2
    assert twice == pytest.approx(base, abs=1e-12)
    assert odd == pytest.approx(base, abs=1e-12)

    # a single item size: the cost of the list is sigma*m identically, so the cap
    # carries no information and the 2h saddle would be singular (note Sec. 5)
    flat = np.full(n, 2)
    _, logE, out = cache_spm_size(g, [m], flat, [60])
    assert out.method == 'uniform-size'
    assert out.zeta[0] == pytest.approx(1.0, abs=1e-12)
    exact = math.log(cache_erec(g, np.array([float(m)]), flat, np.array([60.0])))
    assert abs(logE - exact) < 1e-2
    E, _, _ = cache_spm_size(g, [m], flat, [59])
    assert E == 0.0           # 2*30 = 60 > 59, no state fits


def test_occupancy_and_mean_cost_satisfy_the_saddle_conditions():
    n = 80
    sg = sizes(n)
    _, _, out = cache_spm_size(profile(n), [20], sg, [36])
    assert np.all(out.pij >= 0.0)
    assert out.pij.sum(axis=1) == pytest.approx(np.ones(n), abs=1e-12)
    # the saddle conditions are exactly sum_i pi_i = m and sum_i sigma_i pi_i = k
    assert out.pij[:, 1].sum() == pytest.approx(20.0, abs=1e-9)
    assert (sg * out.pij[:, 1]).sum() == pytest.approx(36.0, abs=1e-9)
    assert out.K[0] == pytest.approx(36.0, abs=1e-9)


def test_sizes_are_required():
    n = 40
    with pytest.raises(ValueError, match="both required"):
        cache_spm_size(profile(n), [10], None, None)
    with pytest.raises(ValueError, match="positive integers"):
        cache_spm_size(profile(n), [10], np.full(n, 1.5), [20])
    with pytest.raises(ValueError, match="'atmost' or 'exact'"):
        cache_spm_size(profile(n), [10], sizes(n), [20], 'bogus')


def _cache_model(n, m, sizes, caps):
    """Source -> Cache -> Sink, Zipf-ish popularity, optional sizes and caps."""
    from line_solver import (Network, Source, Cache, Sink, OpenClass, Exp,
                             DiscreteSampler, ReplacementStrategy)
    model = Network('model')
    source = Source(model, 'Source')
    cache_node = Cache(model, 'Cache', n, m, ReplacementStrategy.RR)
    sink = Sink(model, 'Sink')
    job_class = OpenClass(model, 'InitClass', 0)
    hit_class = OpenClass(model, 'HitClass', 0)
    miss_class = OpenClass(model, 'MissClass', 0)
    source.set_arrival(job_class, Exp(2))
    pop = np.arange(1, n + 1) / np.arange(1, n + 1).sum()
    cache_node.set_read(job_class, DiscreteSampler(pop))
    if sizes is not None:
        cache_node.set_item_sizes(sizes)
    if caps is not None:
        cache_node.set_cost_caps(caps)
    cache_node.set_hit_class(job_class, hit_class)
    cache_node.set_miss_class(job_class, miss_class)
    P = model.init_routing_matrix()
    P.set(job_class, job_class, source, cache_node, 1.0)
    P.set(hit_class, hit_class, cache_node, sink, 1.0)
    P.set(miss_class, miss_class, cache_node, sink, 1.0)
    model.link(P)
    return model, cache_node


def _hit_ratio(node):
    return float(np.ravel(node.get_hit_ratio())[0])


def test_nc_rayint_is_an_alias_of_spm_without_sizes():
    # 'rayint', 'spm' and the default all name one SPM saddle point on a cache.
    from line_solver import SolverNC
    seen = {}
    for method_name in (None, 'spm', 'rayint'):
        mdl, node = _cache_model(14, 4, None, None)
        solver = SolverNC(mdl) if method_name is None else SolverNC(mdl, method=method_name)
        solver.get_avg_table()
        seen[method_name] = (solver.result.method, _hit_ratio(node))
    assert seen[None][0] == 'spm'
    assert seen['spm'] == seen[None]
    assert seen['rayint'] == seen[None]


def test_nc_spm_switches_to_the_size_tilted_kernel_with_item_sizes():
    # With storage costs the same three method names serve cache_spm_size instead, and
    # say so: the label names the kernel that produced the numbers.
    from line_solver import SolverNC
    sizes = [1] * 7 + [2] * 7
    seen = {}
    for method_name in (None, 'spm', 'rayint'):
        mdl, node = _cache_model(14, 4, sizes, None)
        solver = SolverNC(mdl) if method_name is None else SolverNC(mdl, method=method_name)
        solver.get_avg_table()
        seen[method_name] = (solver.result.method, _hit_ratio(node))
    assert seen[None][0] == 'spm.size'
    assert seen['spm'] == seen[None]
    assert seen['rayint'] == seen[None]


def test_nc_spm_size_tracks_exact_with_sizes_and_caps():
    from line_solver import SolverNC
    n, m, sizes, caps = 14, 4, [1] * 7 + [2] * 7, 6
    mdl, node = _cache_model(n, m, sizes, caps)
    SolverNC(mdl, method='exact').get_avg_table()
    exact = _hit_ratio(node)
    mdl2, node2 = _cache_model(n, m, sizes, caps)
    solver = SolverNC(mdl2, method='rayint')
    solver.get_avg_table()
    assert solver.result.method == 'spm.size'
    assert abs(_hit_ratio(node2) - exact) < 0.05


def test_nc_spm_size_degenerates_without_caps():
    from line_solver import SolverNC
    n, m, sizes = 14, 4, [1] * 7 + [2] * 7
    mdl, node = _cache_model(n, m, sizes, None)
    SolverNC(mdl, method='exact').get_avg_table()
    exact = _hit_ratio(node)
    mdl2, node2 = _cache_model(n, m, sizes, None)
    solver = SolverNC(mdl2, method='spm')
    solver.get_avg_table()
    assert solver.result.method == 'spm.size'
    assert abs(_hit_ratio(node2) - exact) < 0.05
