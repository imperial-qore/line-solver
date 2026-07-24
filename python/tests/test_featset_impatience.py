"""
Regression tests: the Reneging and Balking feature gates.

Both names were declared in every feature set but never SET by
getUsedLangFeatures, in any codebase, so they gated nothing: MVA, NC, MAM and
FLD accepted a setPatience/setBalking model and solved it impatience-free,
returning the plain rho/(1-rho) = 9.0 against an exact 1.1565 (CTMC).

Support was established by measurement, not by the declarations: JMT honours
reneging (1.2041) yet declared neither, and the Python JSIM writer emitted no
Balking element at all, so a balking model exported from Python simulated
balking-free. Both are fixed; JMT now declares Reneging and Balking.
"""
import warnings

import pytest

from line_solver import (BalkingStrategy, Exp, Network, OpenClass, Queue,
                         SchedStrategy, Sink, SolverCTMC, SolverJMT,
                         SolverMVA, SolverNC, SolverSSA, Source)

LAMBDA, MU = 0.9, 1.0


def _build(kind):
    model = Network(kind)
    source = Source(model, 'S')
    queue = Queue(model, 'Q', SchedStrategy.FCFS)
    sink = Sink(model, 'K')
    oclass = OpenClass(model, 'C', 0)
    source.setArrival(oclass, Exp(LAMBDA))
    queue.setService(oclass, Exp(MU))
    if kind == 'reneging':
        queue.setPatience(oclass, Exp(0.5))
    elif kind == 'balking':
        queue.setCapacity(10)
        queue.setBalking(oclass, BalkingStrategy.QUEUE_LENGTH, [(3, 10, 1.0)])
    model.link(Network.serialRouting(source, queue, sink))
    return model


def _supports(solver, model):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return solver.supports(model)


def test_features_are_detected():
    assert _build('reneging').getUsedLangFeatures().list['Reneging'] is True
    assert _build('balking').getUsedLangFeatures().list['Balking'] is True
    plain = _build('plain').getUsedLangFeatures()
    assert plain.list['Reneging'] is False
    assert plain.list['Balking'] is False


@pytest.mark.parametrize('kind', ['reneging', 'balking'])
@pytest.mark.parametrize('solver', [SolverMVA, SolverNC])
def test_solvers_that_ignore_impatience_reject(solver, kind):
    assert _supports(solver, _build(kind)) is False


@pytest.mark.parametrize('kind', ['reneging', 'balking'])
@pytest.mark.parametrize('solver', [SolverCTMC, SolverSSA, SolverJMT])
def test_solvers_that_honour_impatience_accept(solver, kind):
    assert _supports(solver, _build(kind)) is True


@pytest.mark.parametrize('solver', [SolverMVA, SolverNC, SolverCTMC, SolverSSA, SolverJMT])
def test_plain_model_is_unaffected(solver):
    # The gate must reject the feature, not the model class.
    assert _supports(solver, _build('plain')) is True


def test_jmt_export_honours_balking():
    # The JSIM writer emitted no Balking element, so balking was silently
    # dropped. CTMC is exact here: 3.9694 without balking, 1.3687 with.
    from line_solver import SolverJMTOptions, VerboseLevel
    options = SolverJMTOptions()
    options.seed = 11
    options.samples = 500000
    options.verbose = VerboseLevel.SILENT
    balking = float(SolverJMT(_build('balking'), options).getAvgTable().QLen[1])
    assert abs(balking - 1.3687) / 1.3687 < 0.05


def test_jmt_export_honours_reneging():
    from line_solver import SolverJMTOptions, VerboseLevel
    options = SolverJMTOptions()
    options.seed = 11
    options.samples = 500000
    options.verbose = VerboseLevel.SILENT
    reneging = float(SolverJMT(_build('reneging'), options).getAvgTable().QLen[1])
    # CTMC exact is 1.1565; simulation noise dominates the tolerance.
    assert abs(reneging - 1.1565) / 1.1565 < 0.10
