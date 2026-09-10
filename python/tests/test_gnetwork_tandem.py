"""Tandem G-network semantics under SolverAG (RCAT) and SolverCTMC.

Model: Source -> Queue1 -> Queue2 -> Sink, with a signal class routed along the
same chain. A removal signal fires ONCE, at the first station it reaches, and is
then annihilated (the state routine never re-emits it, and the MAM AG builder
folds the removal into the local rates of Queue1 only). Queue2 therefore loses
no jobs and its throughput equals its arrival rate.

The catastrophe cases are the regression for the INAP estimator: a catastrophe
makes the isolated process non birth-death, so the mean-of-ratios reversed rate
is no longer the RCAT fixed point and overestimates the departure rate (it
returned a Queue2 throughput of 1.5995 against an arrival rate of 0.8). Those
processes now use the rate-conservation form.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import os

import pytest

from line_solver import (Network, Source, Queue, Sink, OpenClass, Signal,
                         SignalType, SchedStrategy, Exp, SolverAG, SolverMAM, SolverCTMC,
                         SolverLDES)


def _ldes_engine_reason():
    """Reason to skip the LDES leg, or '' when the engine can be run.

    SolverLDES is a subprocess client of the JAR engine, so the leg is skipped
    in a checkout where no engine artifact has been built. The artifacts are
    probed directly rather than through SolverLDES._get_ldes_jar_path(), which
    downloads ldes.jar from SourceForge when it is absent and raises when that
    download fails -- neither belongs in a collection-time predicate.

    A built artifact that cannot be run is NOT a skip: with an engine present,
    a missing java is a misconfiguration, and dropping the coverage silently is
    what lets an engine regression go unnoticed.
    """
    for d in (SolverLDES._get_package_bin_dir(), SolverLDES._get_common_dir()):
        if os.path.isfile(os.path.join(d, 'ldes.jar')):
            return ''
    if SolverLDES._get_ldes_native_path():
        return ''
    return ("no LDES engine artifact; build one with "
            "'mvn clean package -P ldes' in jar/ (writes common/ldes.jar)")


_LDES_SKIP_REASON = _ldes_engine_reason()

MU1 = 2.0
MU2 = 3.0


def _tandem(lambda_pos, lambda_neg, signal_type):
    model = Network('GNetTandem')
    source = Source(model, 'Source')
    queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')

    pos = OpenClass(model, 'Positive')
    source.setArrival(pos, Exp(lambda_pos))
    queue1.setService(pos, Exp(MU1))
    queue2.setService(pos, Exp(MU2))

    neg = Signal(model, 'Negative', signal_type)
    source.setArrival(neg, Exp(lambda_neg))
    queue1.setService(neg, Exp(MU1))
    queue2.setService(neg, Exp(MU2))

    P = model.initRoutingMatrix()
    for c in (pos, neg):
        P.set(c, c, source, queue1, 1.0)
        P.set(c, c, queue1, queue2, 1.0)
        P.set(c, c, queue2, sink, 1.0)
    model.link(P)
    return model


def _pick(table, station, column):
    """Positive-class entry of a station row of an average table.

    Solvers return either an IndexedTable (which wraps the frame in .data) or a
    bare DataFrame, so both shapes are accepted.
    """
    data = getattr(table, 'data', table)
    row = data[(data['Station'] == station) & (data['JobClass'] == 'Positive')]
    return float(row[column].iloc[0])


# (name, signal type, lambda+, lambda-, Q1, Q2, T1)
CASES = [
    ('negative', SignalType.NEGATIVE, 1.0, 0.3, 0.769231, 0.408163, 0.869565),
    ('catastrophe', SignalType.CATASTROPHE, 1.0, 0.3, 0.666667, 0.363636, 0.8),
    ('catastrophe_strong', SignalType.CATASTROPHE, 1.0, 1.0, 0.414214, 0.242641, 0.585786),
]


@pytest.mark.parametrize('name,signal_type,lambda_pos,lambda_neg,q1,q2,t1', CASES,
                         ids=[c[0] for c in CASES])
@pytest.mark.parametrize('method', ['inap', 'inapinf'])
def test_mam_tandem(method, name, signal_type, lambda_pos, lambda_neg, q1, q2, t1):
    table = SolverAG(_tandem(lambda_pos, lambda_neg, signal_type), method).getAvgTable()
    assert _pick(table, 'Queue1', 'QLen') == pytest.approx(q1, abs=1e-4)
    assert _pick(table, 'Queue2', 'QLen') == pytest.approx(q2, abs=1e-4)
    assert _pick(table, 'Queue1', 'Tput') == pytest.approx(t1, abs=1e-4)
    # The signal fires at Queue1 only, so Queue2 must not lose jobs.
    assert _pick(table, 'Queue2', 'Tput') == pytest.approx(
        _pick(table, 'Queue1', 'Tput'), abs=1e-4)


@pytest.mark.parametrize('name,signal_type,lambda_pos,lambda_neg,q1,q2,t1', CASES,
                         ids=[c[0] for c in CASES])
@pytest.mark.skipif(bool(_LDES_SKIP_REASON), reason=_LDES_SKIP_REASON)
def test_ldes_tandem(name, signal_type, lambda_pos, lambda_neg, q1, q2, t1):
    """Sample-path check: the signal is annihilated at Queue1, and a
    catastrophe empties the station including the in-service job."""
    table = SolverLDES(_tandem(lambda_pos, lambda_neg, signal_type),
                       samples=500000, seed=23000).getAvgTable()
    assert _pick(table, 'Queue1', 'QLen') == pytest.approx(q1, abs=0.03)
    assert _pick(table, 'Queue2', 'QLen') == pytest.approx(q2, abs=0.03)
    assert _pick(table, 'Queue1', 'Tput') == pytest.approx(t1, abs=0.02)
    assert _pick(table, 'Queue2', 'Tput') == pytest.approx(
        _pick(table, 'Queue1', 'Tput'), abs=0.01)


@pytest.mark.parametrize('name,signal_type,lambda_pos,lambda_neg,q1,q2,t1', CASES,
                         ids=[c[0] for c in CASES])
def test_ctmc_tandem(name, signal_type, lambda_pos, lambda_neg, q1, q2, t1):
    table = SolverCTMC(_tandem(lambda_pos, lambda_neg, signal_type), cutoff=14).getAvgTable()
    assert _pick(table, 'Queue1', 'QLen') == pytest.approx(q1, abs=1e-3)
    assert _pick(table, 'Queue2', 'QLen') == pytest.approx(q2, abs=1e-3)
    assert _pick(table, 'Queue2', 'Tput') == pytest.approx(
        _pick(table, 'Queue1', 'Tput'), abs=1e-3)


@pytest.mark.parametrize('method', ['default', 'dec.source', 'dec.mmap', 'mna',
                                    'ldqbd', 'bgchain'])
def test_mam_methods_refuse_signals(method):
    """The RCAT builder, which moved to SolverAG, is the only code in LINE that
    reads sn.issignal. No MAM algorithm reads it, so every MAM method would solve
    the model with the signals turned into ordinary customers -- and the MAM
    envelope no longer declares the G-network names at all."""
    solver = SolverMAM(_tandem(1.0, 0.3, SignalType.NEGATIVE), method)
    with pytest.raises(RuntimeError, match='G-network signals'):
        solver.getAvgTable()
    feats = solver.getMethodFeatureSet(method)
    assert 'OpenSignal' not in feats
    assert 'ClosedSignal' not in feats
    assert 'SignalType_NEGATIVE' not in feats


@pytest.mark.parametrize('method', ['inap', 'inapplus', 'inapinf'])
def test_rcat_methods_keep_the_signal_features(method):
    solver = SolverAG(_tandem(1.0, 0.3, SignalType.NEGATIVE), method)
    feats = solver.getMethodFeatureSet(method)
    assert 'OpenSignal' in feats
    assert 'SignalType_NEGATIVE' in feats
    ok, reason = solver.supportsModelMethod(method)
    assert ok, reason


def test_signal_gate_leaves_ordinary_models_alone():
    """A model with no signal class is untouched by the gate."""
    model = Network('Tandem')
    source = Source(model, 'Source')
    queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    pos = OpenClass(model, 'Positive')
    source.setArrival(pos, Exp(1.0))
    queue1.setService(pos, Exp(MU1))
    queue2.setService(pos, Exp(MU2))
    P = model.initRoutingMatrix()
    P.set(pos, pos, source, queue1, 1.0)
    P.set(pos, pos, queue1, queue2, 1.0)
    P.set(pos, pos, queue2, sink, 1.0)
    model.link(P)

    # The claim under test is that the gate does not fire, so the model still
    # reaches the analyzer: Queue1 is the M/M/1 the gate must leave alone.
    table = SolverMAM(model, 'default').getAvgTable()
    assert _pick(table, 'Queue1', 'QLen') == pytest.approx(1.0, abs=1e-4)
