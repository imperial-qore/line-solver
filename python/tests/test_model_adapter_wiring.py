"""Regression tests for the two model transformations that had no entry point
and no caller until they were wired: Network.aggregate_chains (chain
aggregation, ModelAdapter.aggregate_chains) and Network.remove_class /
Network.without_class (class removal, ModelAdapter.remove_class).

Chain aggregation is exact on a product-form model, so the chain totals of the
aggregated model must reproduce the per-class totals of the original, and
sn_deaggregate_chain_results must recover the class-level metrics from them.
Class removal must leave a model identical to one built without that class,
and must leave the original model untouched."""

import numpy as np
import pytest

from line_solver import (
    Network, Queue, Delay, Source, Sink, ClassSwitch, OpenClass, ClosedClass,
    Exp, SchedStrategy, SolverMVA,
)
from line_solver.api.sn.demands import sn_get_demands_chain
from line_solver.api.sn.deaggregate import sn_deaggregate_chain_results

TOL = 1e-6


def _two_class_one_chain():
    """Closed model whose two classes switch into each other: one chain."""
    model = Network('cs2')
    delay = Delay(model, 'Delay')
    queue = Queue(model, 'Queue1', SchedStrategy.PS)
    c1 = ClosedClass(model, 'ClassA', 4, delay)
    c2 = ClosedClass(model, 'ClassB', 0, delay)
    delay.setService(c1, Exp(1.0))
    delay.setService(c2, Exp(2.0))
    queue.setService(c1, Exp(3.0))
    queue.setService(c2, Exp(4.0))
    cs = ClassSwitch(model, 'CS', np.array([[0.3, 0.7], [0.6, 0.4]]))
    P = model.initRoutingMatrix()
    P.set(c1, c1, Network.serialRouting(delay, queue, cs))
    P.set(c2, c2, Network.serialRouting(delay, queue, cs))
    P.set(c1, c1, cs, delay, 1.0)
    P.set(c2, c2, cs, delay, 1.0)
    model.link(P)
    return model, c1, c2


def _open_model(classes=('ClassA', 'ClassB')):
    """Open model with one PS queue and one open class per name given."""
    model = Network('open%d' % len(classes))
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue1', SchedStrategy.PS)
    sink = Sink(model, 'Sink')
    rate = {'ClassA': (0.5, 2.0), 'ClassB': (0.3, 1.5)}
    jobclasses = []
    P = model.initRoutingMatrix()
    for name in classes:
        jclass = OpenClass(model, name)
        source.setArrival(jclass, Exp(rate[name][0]))
        queue.setService(jclass, Exp(rate[name][1]))
        P.set(jclass, Network.serialRouting(source, queue, sink))
        jobclasses.append(jclass)
    model.link(P)
    return model, jobclasses


def _table_column(table, column):
    return np.asarray(table[column], dtype=float)


def test_aggregate_chains_collapses_the_chain():
    model, _, _ = _two_class_one_chain()
    sn = model.get_struct()
    assert sn.nclasses == 2 and sn.nchains == 1

    chain_model, alpha, deagg = model.aggregate_chains()
    assert chain_model.get_number_of_classes() == 1
    assert deagg.is_aggregated
    assert alpha.shape == (sn.nstations, sn.nclasses)
    # alpha weights the classes of a chain and sums to one per station
    np.testing.assert_allclose(alpha.sum(axis=1), np.ones(sn.nstations), atol=TOL)
    # the original model is untouched
    assert model.get_number_of_classes() == 2


def test_aggregate_chains_preserves_chain_totals():
    """Exact on this product-form model: chain metrics must equal the sum of
    the class metrics of the original."""
    model, _, _ = _two_class_one_chain()
    chain_model, _, _ = model.aggregate_chains()

    base = SolverMVA(model).getAvgTable()
    aggr = SolverMVA(chain_model).getAvgTable()

    M = model.get_number_of_stations()
    for column in ('QLen', 'Util', 'Tput'):
        base_col = _table_column(base, column).reshape(M, 2)
        aggr_col = _table_column(aggr, column).reshape(M, 1)
        np.testing.assert_allclose(base_col.sum(axis=1), aggr_col[:, 0], rtol=1e-5, atol=1e-8)


def test_aggregate_chains_deaggregation_recovers_class_metrics():
    model, _, _ = _two_class_one_chain()
    chain_model, _, deagg = model.aggregate_chains()

    sn = model.get_struct()
    M, K = sn.nstations, sn.nclasses
    demands = sn_get_demands_chain(sn)

    aggr = SolverMVA(chain_model).getAvgTable()
    Rchain = _table_column(aggr, 'RespT').reshape(M, 1)
    Tchain = _table_column(aggr, 'Tput').reshape(M, 1)
    refstat = int(np.asarray(sn.refstat).ravel()[0])
    Xchain = np.array([Tchain[refstat, 0]])

    # Qchain and Uchain are left empty on purpose, exactly as solver_mva does:
    # splitting the chain queue length by alpha alone ignores the per-class
    # service times, whereas the response-time branch rescales by ST/STchain
    result = sn_deaggregate_chain_results(
        sn, demands.Lchain, None, demands.STchain, demands.Vchain, deagg.alpha,
        None, None, Rchain, Tchain, None, Xchain)

    base = SolverMVA(model).getAvgTable()
    Qbase = _table_column(base, 'QLen').reshape(M, K)
    Tbase = _table_column(base, 'Tput').reshape(M, K)
    np.testing.assert_allclose(np.asarray(result.Q), Qbase, rtol=1e-4, atol=1e-6)
    np.testing.assert_allclose(np.asarray(result.T), Tbase, rtol=1e-4, atol=1e-6)


def test_aggregate_chains_is_identity_when_each_class_is_a_chain():
    model, _ = _open_model()
    chain_model, alpha, deagg = model.aggregate_chains()
    assert not deagg.is_aggregated
    assert chain_model.get_number_of_classes() == model.get_number_of_classes()

    base = SolverMVA(model).getAvgTable()
    copied = SolverMVA(chain_model).getAvgTable()
    for column in ('QLen', 'Util', 'Tput'):
        np.testing.assert_allclose(_table_column(copied, column),
                                   _table_column(base, column), rtol=1e-6, atol=1e-9)


def test_without_class_matches_a_model_built_without_it():
    model, (_, second) = _open_model()
    reduced = model.without_class(second)
    reference, _ = _open_model(classes=('ClassA',))

    assert reduced.get_number_of_classes() == 1
    assert model.get_number_of_classes() == 2  # the original stays intact

    got = SolverMVA(reduced).getAvgTable()
    want = SolverMVA(reference).getAvgTable()
    for column in ('QLen', 'Util', 'RespT', 'Tput'):
        np.testing.assert_allclose(_table_column(got, column),
                                   _table_column(want, column), rtol=1e-6, atol=1e-9)


def test_remove_class_mutates_in_place():
    model, (_, second) = _open_model()
    model.remove_class(second)
    assert model.get_number_of_classes() == 1
    assert model.get_class_names() == ['ClassA']
    reference, _ = _open_model(classes=('ClassA',))
    np.testing.assert_allclose(
        _table_column(SolverMVA(model).getAvgTable(), 'QLen'),
        _table_column(SolverMVA(reference).getAvgTable(), 'QLen'),
        rtol=1e-6, atol=1e-9)


def test_remove_class_rejects_the_last_class():
    model, (only,) = _open_model(classes=('ClassA',))
    with pytest.raises(RuntimeError):
        model.remove_class(only)


def test_remove_class_by_index():
    model, _ = _open_model()
    model.remove_class(1)
    assert model.get_class_names() == ['ClassA']


def test_remove_class_slices_the_class_switch_matrix():
    """A ClassSwitch matrix is indexed by class position, so removing a class
    must delete its row and column; the model must stay solvable."""
    model = Network('cs3')
    delay = Delay(model, 'Delay')
    queue = Queue(model, 'Queue1', SchedStrategy.PS)
    c1 = ClosedClass(model, 'ClassA', 3, delay)
    c2 = ClosedClass(model, 'ClassB', 0, delay)
    c3 = ClosedClass(model, 'ClassC', 2, delay)
    for jclass, rate in ((c1, 1.0), (c2, 2.0), (c3, 1.5)):
        delay.setService(jclass, Exp(rate))
        queue.setService(jclass, Exp(3.0))
    cs = ClassSwitch(model, 'CS', np.array([[0.3, 0.7, 0.0],
                                            [0.6, 0.4, 0.0],
                                            [0.0, 0.0, 1.0]]))
    P = model.initRoutingMatrix()
    for jclass in (c1, c2, c3):
        P.set(jclass, jclass, Network.serialRouting(delay, queue, cs))
        P.set(jclass, jclass, cs, delay, 1.0)
    model.link(P)
    SolverMVA(model).getAvgTable()

    reduced = model.without_class(c3)
    assert reduced.get_number_of_classes() == 2
    reduced_cs = reduced.get_node_by_name('CS')
    assert reduced_cs.get_class_switching_matrix().shape == (2, 2)
    np.testing.assert_allclose(reduced_cs.get_class_switching_matrix(),
                               np.array([[0.3, 0.7], [0.6, 0.4]]), atol=TOL)
    table = SolverMVA(reduced).getAvgTable()
    assert np.all(np.isfinite(_table_column(table, 'QLen')))


def test_copy_preserves_per_class_settings_beyond_service_and_routing():
    """Network.copy re-keys every per-class dictionary onto the copy's own
    class objects; a class capacity set on the original must survive."""
    model, (first, second) = _open_model()
    queue = model.get_node_by_name('Queue1')
    queue.setClassCapacity(second, 7)
    copied = model.copy()
    copied_queue = copied.get_node_by_name('Queue1')
    copied_second = copied.get_classes()[1]
    assert copied_queue.get_class_capacity(copied_second) == 7
