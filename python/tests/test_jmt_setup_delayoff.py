"""
Regression test: the JMT export must carry Queue setup/delay-off times.

The native Python JSIM writer previously emitted no delayOffTime/setUpTime
parameters at all, so a model built with setDelayOff simulated setup-free and
silently returned M/M/1 results. JMT casts each per-class entry of these
parameters to ServiceStrategy[] and reads element [0], so the strategy must be
nested inside a single-element array; a flat entry aborts the simulation with a
ClassCastException. The reference is the exact QBD solution.
"""
import numpy as np
import pytest

from line_solver import (Exp, Network, OpenClass, Queue, Sink, SchedStrategy,
                         SolverJMT, SolverJMTOptions, Source, VerboseLevel)
from line_solver.api.mam.qbd import qbd_setupdelayoff

LAMBDA, MU = 0.5, 1.0
SETUP_RATE, DELAYOFF_RATE = 2.0, 4.0


def _build(with_delayoff):
    model = Network('setup_delayoff')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1', 0)
    source.setArrival(oclass, Exp(LAMBDA))
    queue.setService(oclass, Exp(MU))
    if with_delayoff:
        queue.setDelayOff(oclass, Exp(SETUP_RATE), Exp(DELAYOFF_RATE))
    P = model.initRoutingMatrix()
    P.set(oclass, oclass, source, queue, 1.0)
    P.set(oclass, oclass, queue, sink, 1.0)
    model.link(P)
    return model


def _jmt_qlen(model):
    options = SolverJMTOptions()
    options.seed = 23000
    options.samples = 500000
    options.verbose = VerboseLevel.SILENT
    return float(SolverJMT(model, options).getAvgTable().QLen[1])


def test_export_nests_strategies_in_class_array(tmp_path):
    import xml.etree.ElementTree as ET

    from line_solver.api.solvers.jmt import handler

    model = _build(True)
    path = str(tmp_path / 'setup_delayoff.jsimg')
    handler._write_jsim_file(model.getStruct(), path, SolverJMTOptions(), model)
    root = ET.parse(path).getroot()

    for name in ('delayOffTime', 'setUpTime'):
        params = root.findall(".//parameter[@name='%s']" % name)
        assert len(params) == 1, '%s must be exported exactly once' % name
        rows = params[0].findall("subParameter[@name='%s']" % name)
        assert len(rows) == 1, '%s must carry one row per class' % name
        # JMT reads setUpStrategies[class][0], so the row itself is an array.
        assert rows[0].get('array') == 'true'
        assert rows[0].get('classPath') == 'jmt.engine.NetStrategies.ServiceStrategy'
        assert rows[0].find("subParameter[@name='ServiceTimeStrategy']") is not None


def test_delayoff_matches_qbd():
    reference = qbd_setupdelayoff(LAMBDA, MU, SETUP_RATE, 1.0, DELAYOFF_RATE, 1.0)
    qlen = _jmt_qlen(_build(True))
    assert np.abs(qlen - reference) / reference < 0.05


def test_delayoff_is_not_ignored():
    # Without setup the station is a plain M/M/1 with QLen = rho/(1-rho) = 1.
    plain = _jmt_qlen(_build(False))
    assert np.abs(plain - 1.0) < 0.05
    assert _jmt_qlen(_build(True)) > plain * 1.1
