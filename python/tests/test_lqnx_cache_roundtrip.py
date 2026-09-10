"""Round-trip of the LINE .lqnx dialect: caches, item entries, setup times.

The stock LQN schema has no element for a cache, an item entry, or a task's
setup / delay-off time, so `writeXML` used to emit a valid file describing a
DIFFERENT model: a CacheTask became a plain task, an ItemEntry a plain entry,
and the hit/miss split a bare unweighted <post> that every reader takes as
50/50. Every engine then answered correctly for what it was handed and
disagreed with a golden computed from the real model.

COUNTS ARE NOT AN ACCEPTANCE TEST. The JAR's existing round-trip asserts only
nhosts/ntasks/nentries, which is exactly why it stayed green while the export
destroyed the model. Every assertion here is on a VALUE that the old writer
silently dropped.

Wire format: scratchpad LQNX_CACHE_SPEC.md, whose field set is the JSON
interchange's (linemodel_save.m).
"""
import os
import tempfile
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from line_solver import (Activity, ActivityPrecedence, CacheTask, DiscreteSampler, Entry,
                         Exp, GlobalConstants, ItemEntry, LayeredNetwork, Processor,
                         ReplacementStrategy, SchedStrategy, SetupTask, Task)


def _cache_model(caps=(1, 1), items=4, replacement=ReplacementStrategy.RR,
                 popularity=None):
    """Client -> CacheTask, the read forking into a hit and a miss activity."""
    m = LayeredNetwork('cacheRT')
    p1 = Processor(m, 'P1', 1, SchedStrategy.PS)
    t1 = Task(m, 'T1', 1, SchedStrategy.REF).on(p1)
    e1 = Entry(m, 'E1').on(t1)
    pc = Processor(m, 'Pc', 1, SchedStrategy.PS)
    ct = CacheTask(m, 'CT', items, list(caps), replacement, 1).on(pc)
    if popularity is None:
        popularity = DiscreteSampler(np.ones(items) / items)
    ie = ItemEntry(m, 'IE', items, popularity).on(ct)
    a1 = Activity(m, 'A1', Exp(1.0)).on(t1).bound_to(e1).synch_call(ie, 1)
    ac = Activity(m, 'Ac', Exp(1e8)).on(ct).bound_to(ie)
    hit = Activity(m, 'Ac_hit', Exp(1.0)).on(ct)
    miss = Activity(m, 'Ac_miss', Exp(0.5)).on(ct)
    ct.add_precedence(ActivityPrecedence.CacheAccess(ac, [hit, miss]))
    return m


def _setup_model(setup_mean=0.5, setup_scv=1.0, off_mean=2.0, off_scv=1.0):
    m = LayeredNetwork('setupRT')
    p1 = Processor(m, 'P1', 1, SchedStrategy.PS)
    t1 = Task(m, 'T1', 1, SchedStrategy.REF).on(p1)
    e1 = Entry(m, 'E1').on(t1)
    p2 = Processor(m, 'P2', 4, SchedStrategy.FCFS)
    f2 = SetupTask(m, 'F2', 6, SchedStrategy.FCFS).on(p2)
    f2.setSetupTime(Exp(1.0 / setup_mean))
    f2.setDelayOffTime(Exp(1.0 / off_mean))
    e2 = Entry(m, 'E2').on(f2)
    Activity(m, 'A1', Exp(1.0)).on(t1).bound_to(e1).synch_call(e2, 1)
    Activity(m, 'A2', Exp(3.0)).on(f2).bound_to(e2).replies_to(e2)
    return m


def _roundtrip(model):
    path = os.path.join(tempfile.mkdtemp(prefix='lqnx_rt_'), 'm.lqnx')
    model.writeXML(path)
    return path, LayeredNetwork.parse_xml(path)


def _task(model, name):
    return next(t for t in model.tasks if t.name == name)


def _entry(model, name):
    return next(e for e in model.entries if e.name == name)


# --- the cache -------------------------------------------------------------

def test_cache_task_survives_roundtrip():
    """A CacheTask must come back a CacheTask, not a plain Task."""
    _, back = _roundtrip(_cache_model())
    assert type(_task(back, 'CT')).__name__ == 'CacheTask'


def test_cache_items_capacity_and_replacement_survive():
    src = _cache_model(caps=(2, 1), items=7, replacement=ReplacementStrategy.LRU)
    _, back = _roundtrip(src)
    ct = _task(back, 'CT')
    assert int(ct.total_items) == 7
    # THE CAPACITY IS AN ARRAY. A multi-list cache is the normal case, and
    # collapsing it to a scalar silently changes the cache.
    assert [int(c) for c in np.atleast_1d(ct.cache_capacity)] == [2, 1]
    assert ct.replacement_strategy == ReplacementStrategy.LRU


def test_single_list_cache_capacity_survives():
    _, back = _roundtrip(_cache_model(caps=(3,), items=5))
    ct = _task(back, 'CT')
    assert [int(c) for c in np.atleast_1d(ct.cache_capacity)] == [3]


@pytest.mark.parametrize('strategy', [ReplacementStrategy.RR, ReplacementStrategy.FIFO,
                                      ReplacementStrategy.LRU, ReplacementStrategy.SFIFO])
def test_every_replacement_spelling_survives(strategy):
    _, back = _roundtrip(_cache_model(replacement=strategy))
    assert _task(back, 'CT').replacement_strategy == strategy


def test_cache_element_shape_matches_spec():
    """One <level capacity=> per list, in order, under a <cache> on the task."""
    path, _ = _roundtrip(_cache_model(caps=(2, 1), items=7))
    root = ET.parse(path).getroot()
    cache = root.find(".//task[@name='CT']/cache")
    assert cache is not None
    assert cache.get('items') == '7'
    assert cache.get('replacement') == 'RR'
    assert [lv.get('capacity') for lv in cache.findall('./level')] == ['2', '1']


# --- the item entry --------------------------------------------------------

def test_item_entry_and_cardinality_survive():
    _, back = _roundtrip(_cache_model(items=6))
    ie = _entry(back, 'IE')
    assert type(ie).__name__ == 'ItemEntry'
    assert int(ie.total_items) == 6


def test_access_popularity_values_survive():
    """The popularity is the access law; a uniform default would be a different
    model, so the probabilities themselves must come back."""
    probs = np.array([0.5, 0.25, 0.15, 0.10])
    _, back = _roundtrip(_cache_model(items=4, popularity=DiscreteSampler(probs)))
    ie = _entry(back, 'IE')
    assert np.allclose(np.asarray(ie.access_prob._probs).flatten(), probs)


def test_non_default_support_survives():
    """A DiscreteSampler whose support is not 1..n writes p then x; the reader
    tells the two apart by the cardinality."""
    probs = np.array([0.4, 0.6])
    values = np.array([3.0, 7.0])
    _, back = _roundtrip(_cache_model(items=2, popularity=DiscreteSampler(probs, values)))
    ie = _entry(back, 'IE')
    assert np.allclose(np.asarray(ie.access_prob._probs).flatten(), probs)
    assert np.allclose(np.asarray(ie.access_prob._values).flatten(), values)


# --- the hit / miss split --------------------------------------------------

def test_post_cache_precedence_survives_with_hit_miss_assignment():
    """The split must load as a CACHE_ACCESS whose outcomes are (hit, miss) in
    that order. A bare <post> loads as an unweighted branch, which is where the
    spurious 'exact 0.5' hit ratio came from."""
    _, back = _roundtrip(_cache_model())
    ct = _task(back, 'CT')
    precs = [p for p in ct.precedences
             if getattr(p.prec_type, 'name', '') == 'CACHE_ACCESS']
    assert len(precs) == 1, 'cache access precedence lost on round trip'
    prec = precs[0]
    assert [a.name for a in prec.pre_activities] == ['Ac']
    assert [a.name for a in prec.post_activities] == ['Ac_hit', 'Ac_miss']


def test_post_cache_result_attributes_written():
    path, _ = _roundtrip(_cache_model())
    root = ET.parse(path).getroot()
    post = root.find(".//post-CACHE")
    assert post is not None, 'cache access still written as a bare <post>'
    assert [(a.get('name'), a.get('cache-result'))
            for a in post.findall('./activity')] == [('Ac_hit', 'hit'), ('Ac_miss', 'miss')]


def test_cache_result_falls_back_to_document_order():
    """A file written before cache-result existed must still load, hit first."""
    path, _ = _roundtrip(_cache_model())
    tree = ET.parse(path)
    for act in tree.getroot().find('.//post-CACHE').findall('./activity'):
        del act.attrib['cache-result']
    tree.write(path)
    back = LayeredNetwork.parse_xml(path)
    prec = next(p for p in _task(back, 'CT').precedences
                if getattr(p.prec_type, 'name', '') == 'CACHE_ACCESS')
    assert [a.name for a in prec.post_activities] == ['Ac_hit', 'Ac_miss']


# --- setup and delay-off ---------------------------------------------------

def test_setup_and_delay_off_survive():
    _, back = _roundtrip(_setup_model(setup_mean=0.5, off_mean=2.0))
    f2 = _task(back, 'F2')
    assert f2.setup_time is not None, 'setup time dropped on round trip'
    assert f2.delay_off_time is not None, 'delay-off time dropped on round trip'
    assert f2.setup_time.getMean() == pytest.approx(0.5, rel=1e-9)
    assert f2.delay_off_time.getMean() == pytest.approx(2.0, rel=1e-9)


def test_setup_scv_survives():
    """A mean alone cannot rebuild the distribution family, so the SCV rides too."""
    m = _setup_model()
    from line_solver import Erlang
    _task(m, 'F2').setSetupTime(Erlang.fit_mean_and_scv(0.4, 0.25))
    _, back = _roundtrip(m)
    st = _task(back, 'F2').setup_time
    assert st.getMean() == pytest.approx(0.4, rel=1e-9)
    assert st.getSCV() == pytest.approx(0.25, rel=1e-6)


def test_setup_elements_absent_when_unset():
    """A plain task must not gain a <setup>; the elements are emitted only when
    the corresponding time is set."""
    m = LayeredNetwork('plain')
    p1 = Processor(m, 'P1', 1, SchedStrategy.PS)
    t1 = Task(m, 'T1', 1, SchedStrategy.REF).on(p1)
    e1 = Entry(m, 'E1').on(t1)
    Activity(m, 'A1', Exp(1.0)).on(t1).bound_to(e1)
    path, back = _roundtrip(m)
    root = ET.parse(path).getroot()
    assert root.find('.//setup') is None
    assert root.find('.//delay-off') is None
    assert root.find('.//cache') is None
    assert type(_task(back, 'T1')).__name__ == 'Task'


# --- tolerance -------------------------------------------------------------

def test_retrieval_flag_survives():
    m = _cache_model()
    _task(m, 'CT').retrieval = True
    path, back = _roundtrip(m)
    assert ET.parse(path).getroot().find(".//cache").get('retrieval') == 'true'
    assert _task(back, 'CT').retrieval is True


def test_zipf_popularity_survives():
    """Zipf is the second carried class: two parameters, s then n. Python holds
    `_s`/`_n` directly, so MATLAB's 1=p,2=x,3=s,4=n layout trap does not arise
    here -- but the wire form is the same two values, in that order."""
    from line_solver import Zipf
    path, back = _roundtrip(_cache_model(items=5, popularity=Zipf(1.2, 5)))
    ap = ET.parse(path).getroot().find('.//access-popularity')
    assert ap.get('name') == 'Zipf'
    assert [p.get('value') for p in ap.findall('./parameter')] == [repr(1.2), repr(5.0)]
    pop = _entry(back, 'IE').access_prob
    assert type(pop).__name__ == 'Zipf'
    assert float(pop.s) == pytest.approx(1.2)
    assert int(pop.n) == 5


def test_zipf_cardinality_split_is_not_applied():
    """The `name` attribute selects the rule, so a 2-item Zipf is never read as a
    DiscreteSampler and a 2-item DiscreteSampler is never read as a Zipf."""
    from line_solver import Zipf
    _, back = _roundtrip(_cache_model(items=2, popularity=Zipf(0.8, 2)))
    assert type(_entry(back, 'IE').access_prob).__name__ == 'Zipf'
    _, back2 = _roundtrip(_cache_model(items=2,
                                       popularity=DiscreteSampler(np.array([0.3, 0.7]))))
    assert type(_entry(back2, 'IE').access_prob).__name__ == 'DiscreteSampler'


def test_zipf_with_wrong_parameter_count_is_refused():
    """A reader meeting name='Zipf' with anything but two parameters refuses by
    name rather than guessing."""
    from line_solver import Zipf
    path, _ = _roundtrip(_cache_model(items=5, popularity=Zipf(1.2, 5)))
    tree = ET.parse(path)
    ap = tree.getroot().find('.//access-popularity')
    ET.SubElement(ap, 'parameter').set('value', '3.0')
    tree.write(path)
    with pytest.raises(RuntimeError, match=r"'Zipf' takes exactly two parameters"):
        LayeredNetwork.parse_xml(path)


def test_unknown_popularity_class_is_refused_on_read():
    path, _ = _roundtrip(_cache_model())
    tree = ET.parse(path)
    tree.getroot().find('.//access-popularity').set('name', 'Weibull')
    tree.write(path)
    with pytest.raises(RuntimeError, match=r"names the class 'Weibull'"):
        LayeredNetwork.parse_xml(path)


def test_unencodable_popularity_is_refused_by_name():
    """A popularity class this encoding cannot express is refused NAMING THE
    CLASS. Writing a nameless element without its parameters would reproduce the
    very bug the dialect exists to fix."""
    from line_solver import Exp as _Exp
    m = _cache_model()
    _entry(m, 'IE').access_prob = _Exp(1.0)
    with pytest.raises(RuntimeError, match=r"access-popularity.*'Exp'"):
        _roundtrip(m)


# --- SEMANTIC acceptance: the numbers, not just the fields ------------------
#
# The structural assertions above prove the FIELDS survive. They would still
# pass if a field landed somewhere the solver never reads. Solving the model on
# both sides of the round trip is what proves the dialect carries the MODEL.

def _solve(model):
    """Every metric of the layered table, keyed by element and metric name."""
    from line_solver import LN, MVA
    opts = MVA.default_options()
    opts.verbose = 0
    ln_opts = LN.default_options()
    ln_opts.verbose = 0
    df = LN(model, lambda m: MVA(m, opts), ln_opts).get_avg_table()
    metrics = [c for c in ('QLen', 'Util', 'RespT', 'ResidT', 'Tput') if c in df.columns]
    return {(r['Node'], c): float(r[c]) for _, r in df.iterrows() for c in metrics}, df


def _assert_same_solution(before, after):
    """Same elements, same numbers, and the SAME NaN MASK: which metrics are
    undefined for a Processor or a Task is part of the answer, so a NaN that
    turns into a number is a difference, not a tolerance question."""
    assert set(before) == set(after), 'round trip changed which elements exist'
    for key in sorted(before):
        b, a = before[key], after[key]
        if np.isnan(b) or np.isnan(a):
            assert np.isnan(b) and np.isnan(a), \
                'round trip changed the NaN mask at %s: %s -> %s' % (key, b, a)
            continue
        assert a == pytest.approx(b, rel=1e-9, abs=1e-12), \
            'round trip changed %s: %s -> %s' % (key, b, a)


def test_cache_model_solves_identically_after_roundtrip():
    src = _cache_model()
    _, back = _roundtrip(src)
    before, df_before = _solve(src)
    after, df_after = _solve(back)
    _assert_same_solution(before, after)
    # The hit and miss branches must not come back as an even split: an
    # unweighted <post> is exactly what the old exporter produced.
    tput = {r['Node']: r['Tput'] for _, r in df_after.iterrows()}
    assert tput['Ac'] == pytest.approx(0.0, abs=1e-9), \
        'cache access activity carries flow; the split degenerated'


def test_setup_model_solves_identically_after_roundtrip():
    src = _setup_model()
    _, back = _roundtrip(src)
    before, _ = _solve(src)
    after, _ = _solve(back)
    _assert_same_solution(before, after)


def _lcq_singlehost_model():
    """The real cache-LQN example, built from its own source so the test cannot
    drift from the model the golden was computed on."""
    here = os.path.dirname(os.path.abspath(__file__))
    src_path = os.path.join(here, '..', 'examples', 'advanced', 'layeredCQ',
                            'lcq_singlehost.py')
    if not os.path.exists(src_path):
        pytest.skip('lcq_singlehost example not present')
    src = open(src_path).read()
    cut = src.find('LN(')
    ns = {'__name__': 'lcq_probe'}
    exec(compile(src[:src.rfind('\n', 0, cut)], src_path, 'exec'), ns)
    models = [v for v in ns.values() if type(v).__name__ == 'LayeredNetwork']
    if not models:
        pytest.skip('could not build lcq_singlehost model')
    return models[0]


def test_real_cache_example_matches_golden_after_roundtrip():
    """The half that matters: solve the SHIPPED cache example after a round trip
    and land on the parity golden, not merely on the same fields."""
    here = os.path.dirname(os.path.abspath(__file__))
    # The shared goldens live IN THIS REPO, at line-dev.git/goldens (moved out of
    # line-test.git/parity-static on 2026-08-19). This test used to reach for a
    # SIBLING checkout and skip when it was absent -- a hole that reads exactly
    # like a pass, which is what the move removes. It is now unconditional.
    goldens = os.environ.get('LINE_GOLDENS') or os.path.join(
        os.path.dirname(os.path.dirname(here)), 'goldens')
    golden_path = os.path.join(goldens, 'baselines', 'lcq_singlehost.json')
    import json
    golden = json.load(open(golden_path))['solvers']['LN']

    _, back = _roundtrip(_lcq_singlehost_model())
    solved, _ = _solve(back)
    checked = 0
    for row in golden:
        node = row['Station']
        for metric in ('QLen', 'Util', 'RespT', 'ResidT', 'Tput'):
            expected = row.get(metric)
            if expected is None or expected == 'NaN':
                continue
            key = (node, metric)
            assert key in solved, 'element %s lost on round trip' % node
            assert solved[key] == pytest.approx(float(expected), rel=1e-4, abs=1e-6), \
                'round-tripped %s %s = %s, golden %s' % (node, metric, solved[key], expected)
            checked += 1
    assert checked > 10, 'golden comparison asserted almost nothing (%d cells)' % checked


def test_reader_tolerates_unknown_elements():
    """Rule 4 of the spec: a reader meeting an element it does not know must not
    crash, or the dialect cannot be extended again."""
    path, _ = _roundtrip(_cache_model())
    tree = ET.parse(path)
    task = tree.getroot().find(".//task[@name='CT']")
    ET.SubElement(task, 'future-extension').set('whatever', '1')
    tree.write(path)
    back = LayeredNetwork.parse_xml(path)
    assert type(_task(back, 'CT')).__name__ == 'CacheTask'
