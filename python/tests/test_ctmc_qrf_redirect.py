"""The QRF reduction bounds moved out of SolverCTMC into SolverBA, and the
refusal has to carry that forwarding address.

Dropping the names from ``listValidMethods`` is what makes SolverCTMC refuse
them, and it is also what loses the address, so the two have to be declared
together. They were not: ``runAnalyzerChecks`` sits above every dispatcher and
reported the flat "the 'qrf.bas' method is unsupported by this solver". The gate
asks ``unsupportedMethodReason`` first now.

The python entry point had additionally OUTLIVED the move: it kept dispatching
to ``solver_ctmc_qrf_analyzer``, the very analyzer SolverBA itself calls
(``solver_ba_analyzer:435``), so it was a second front door onto one
computation. ``test_solver_ba_gives_what_the_ctmc_path_used_to`` pins that the
address answers with the same numbers, which is what makes retiring the door
free.
"""

import numpy as np
import pytest

from line_solver import (Network, Queue, ClosedClass, Exp, SchedStrategy,
                         SolverBA, SolverCTMC)

QRF_METHODS = ['qrf', 'qrf.bas', 'qrf.mmi', 'qrf.rsrd']


def _cqn():
    """A two-queue closed cycle: no delay station, which QRF refuses."""
    m = Network('cqn')
    q1 = Queue(m, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Q2', SchedStrategy.FCFS)
    c = ClosedClass(m, 'C', 2, q1)
    q1.setService(c, Exp(1.0))
    q2.setService(c, Exp(2.0))
    m.link(Network.serialRouting(q1, q2))
    return m


@pytest.mark.parametrize('moved', QRF_METHODS)
def test_qrf_methods_redirect_to_solver_ba(moved):
    with pytest.raises(Exception) as exc:
        SolverCTMC(_cqn(), method=moved).getAvgTable()
    msg = str(exc.value)
    assert 'moved out of SolverCTMC into the dedicated SolverBA solver' in msg
    assert "SolverBA(model, '%s')" % moved in msg


def test_the_redirect_also_fires_with_checks_disabled():
    """enableChecks=False skips the name gate, so runAnalyzer must refuse too --
    otherwise the retired entry point is still reachable by that route."""
    s = SolverCTMC(_cqn(), method='qrf.bas')
    s.enableChecks = False
    with pytest.raises(Exception) as exc:
        s.getAvgTable()
    assert 'moved out of SolverCTMC' in str(exc.value)


def test_unknown_method_keeps_the_flat_refusal():
    with pytest.raises(Exception) as exc:
        SolverCTMC(_cqn(), method='nosuchmethod').getAvgTable()
    msg = str(exc.value)
    assert 'unsupported by this solver' in msg
    assert 'moved out of SolverCTMC' not in msg


def test_solver_ba_gives_what_the_ctmc_path_used_to():
    """The address answers, and with the same numbers the retired CTMC entry
    point produced: measured QLen [1.7778, 0.22222] on this model either way."""
    table = SolverBA(_cqn(), method='qrf.bas').getAvgTable()
    qlen = np.asarray(table['QLen'], dtype=float).ravel()
    np.testing.assert_allclose(qlen, [1.77778, 0.22222], rtol=1e-4)
