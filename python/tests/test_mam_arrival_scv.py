"""
Regression test: the MAM dec.source arrival MMAP must preserve the
arrival variability of open chains. A defect previously replaced the
source arrival process with a Poisson MMAP, making PH/PH/1 results
arrival-scv-independent (all equal to the M/PH/1 value). References are
the exact CTMC values.
"""
import warnings

import numpy as np
import pytest

from line_solver import APH, Network, OpenClass, Queue, Sink, SolverMAM, Source


def _build(ca):
    model = Network('g')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue')
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'C')
    source.setArrival(oclass, APH.fitMeanAndSCV(2.0, ca))
    queue.setService(oclass, APH.fitMeanAndSCV(1.0, 4.0))
    model.link(Network.serialRouting(source, queue, sink))
    return model


def _qlen(model):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        q = np.asarray(SolverMAM(model).getAvgQLen())
    return q[1, 0] if q.ndim > 1 else q[1]


class TestMamArrivalScv:
    def test_phph1_arrival_scv_dependence(self):
        assert _qlen(_build(2.0)) == pytest.approx(2.08423, abs=5e-2)
        assert _qlen(_build(16.0)) == pytest.approx(4.62730, abs=2e-1)
