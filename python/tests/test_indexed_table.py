"""
Regression tests for IndexedTable indexing semantics.

IndexedTable overloads __getitem__ so that a solver result table can be filtered
by model objects (``table[queue, jobclass]``). That overload must not capture
ordinary pandas indexers. It used to: ``_is_object_arg`` treated anything that is
not a builtin scalar or container as a model object, so a boolean ``pd.Series``
was routed to ``filterBy``, matched its ``hasattr(obj, 'name')`` fallback -- a
Series carries a ``.name`` attribute -- and came back as a SILENTLY EMPTY frame.
``table[table['Station'] == 'Q']`` therefore returned nothing while
``table.data[...]`` returned the row, which is what broke both shipped
``examples/advanced/passAndSwap`` scripts AFTER they had solved correctly.

Silent emptiness is the dangerous part: it reads as "no matching rows" rather
than as a failure, so these assert row counts against ``.data`` directly.
"""

import numpy as np

from line_solver import Network, Queue, Delay, ClosedClass, SchedStrategy, Exp, CTMC


def _table():
    model = Network('IndexedTableProbe')
    d = Delay(model, 'D')
    q = Queue(model, 'Q', SchedStrategy.FCFS)
    cl = ClosedClass(model, 'C', 2, d)
    d.setService(cl, Exp(1.0))
    q.setService(cl, Exp(2.0))
    model.link(Network.serialRouting(d, q))
    return CTMC(model).getAvgTable(), model, d, q, cl


def test_boolean_mask_filters_like_a_dataframe():
    t, _, _, _, _ = _table()
    mask = t['Station'] == 'Q'
    assert int(np.asarray(mask).sum()) == 1, 'probe model should have exactly one Q row'
    assert len(t[mask]) == 1, 'boolean mask must not be swallowed by object filtering'
    assert len(t[mask]) == len(t.data[np.asarray(mask)]), \
        'mask filtering must agree with the underlying DataFrame'


def test_compound_boolean_mask():
    t, _, _, _, _ = _table()
    sub = t[(t['Station'] == 'Q') & (t['JobClass'] == 'C')]
    assert len(sub) == 1
    # the value must be reachable the way the shipped examples read it
    assert np.isfinite(float(sub['QLen'].iloc[0]))


def test_empty_mask_is_empty_not_an_error():
    t, _, _, _, _ = _table()
    assert len(t[t['Station'] == 'NoSuchStation']) == 0


def test_object_filtering_still_works():
    t, _, _, q, cl = _table()
    assert len(t[q, cl]) == 1
    assert len(t.filterBy(q, cl)) == 1


def test_column_and_positional_access_unaffected():
    t, _, _, _, _ = _table()
    assert sorted(set(t['Station'])) == ['D', 'Q']
    assert list(t[['Station', 'QLen']].columns) == ['Station', 'QLen']
    assert len(t[0:1]) == 1


def test_station_names_are_not_quoted():
    """Guards the misreading that names carry embedded quotes: they do not."""
    t, _, _, _, _ = _table()
    for nm in t['Station']:
        assert not str(nm).startswith("'"), f'station name unexpectedly quoted: {nm!r}'
        assert str(nm) in ('D', 'Q')
