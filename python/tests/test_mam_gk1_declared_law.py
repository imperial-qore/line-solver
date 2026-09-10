"""
The MMAP[K]/G[K]/1 arm must read the law the user DECLARED, not the phase-type
fit that replaces it.

``sn_nonmarkov_toph`` runs inside ``SolverMAM.runAnalyzer``, before the dispatch,
and retags ``sn.procid`` PH or ME; DET is the only tag it preserves. A gate
reading the live ``procid`` therefore held for a station with a Det beside the
other law and fell to the fit everywhere else: an M/Gamma/1 at rho 0.5 with SCV 3
read 1.0785538590 against 1.5.

Pollaczek-Khinchine needs only the first two moments of an M/G/1 service, so it
is an oracle independent of every path in the solver.

Mirrors jar SolverMamGeneralServiceTest.declaredLawSelectsTheTransformAnalysis-
WithoutADetBesideIt and the cpp test_mam_gk1_solver case of the same name.
"""

import numpy as np
import pytest

from line_solver import (Network, Source, Queue, Sink, OpenClass, Exp, Det,
                         Uniform, Gamma, Lognormal, SchedStrategy, SolverMAM)


def _qlen(service, arr_rate=0.5):
    model = Network('gk1_declared')
    src = Source(model, 'Source')
    qu = Queue(model, 'Queue', SchedStrategy.FCFS)
    snk = Sink(model, 'Sink')
    oc = OpenClass(model, 'Class1')
    src.setArrival(oc, Exp(arr_rate))
    qu.setService(oc, service)
    model.link(Network.serialRouting(src, qu, snk))
    return float(SolverMAM(model).getAvgTable().QLen[1])


def _pk(lam, m1, scv):
    rho = lam * m1
    return rho + lam * lam * (m1 * m1 * (1.0 + scv)) / (2.0 * (1.0 - rho))


@pytest.mark.filterwarnings('ignore::UserWarning')
@pytest.mark.parametrize('law', [
    Det(1.0),
    Uniform(0.2, 1.8),
    Gamma.fitMeanAndSCV(1.0, 3.0),
    Lognormal.fitMeanAndSCV(1.0, 2.0),
])
def test_declared_law_reaches_the_transform_analysis(law):
    expected = _pk(0.5, float(law.getMean()), float(law.getSCV()))
    assert abs(_qlen(law) - expected) < 1e-6 * max(1.0, abs(expected))


@pytest.mark.filterwarnings('ignore::UserWarning')
def test_declared_tags_are_snapshotted_before_the_conversion():
    # The snapshot is what the gate reads; without it the live tag is PH here.
    from line_solver.constants import ProcessType
    model = Network('snapshot')
    src = Source(model, 'Source')
    qu = Queue(model, 'Queue', SchedStrategy.FCFS)
    snk = Sink(model, 'Sink')
    oc = OpenClass(model, 'Class1')
    src.setArrival(oc, Exp(0.5))
    qu.setService(oc, Gamma.fitMeanAndSCV(1.0, 3.0))
    model.link(Network.serialRouting(src, qu, snk))
    solver = SolverMAM(model)
    solver.getAvgTable()
    declared = getattr(solver.sn, 'procid_declared', None)
    assert declared is not None
    assert declared[1, 0] == ProcessType.GAMMA
    assert solver.sn.procid[1, 0] == ProcessType.PH
