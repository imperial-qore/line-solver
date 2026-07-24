"""
A multiclass station whose service is matrix-exponential falls outside the exact
RAP/RAP/1 QBD and is answered by the phase-type approximation MMAPPH1FCFS, which
is not exact for a non-phase-type service. The user is told so, and must be told
so for every model concerned: this warning reports a correctness limitation, so
suppressing it as a repeat would let the second and later models in a session
read as clean.

Mirrors jar/src/test/java/jline/solvers/mam/MamMEWarningTest.java.
"""

import io
import os
from contextlib import redirect_stdout

import numpy as np
import pytest

from line_solver import (Network, Source, Queue, Sink, OpenClass, Exp, ME,
                         SchedStrategy, SolverMAM)
from line_solver.api.io.logging import line_warning_always


# A matrix-exponential with a negative entry in alpha, whose density has an
# interior zero, so it admits no phase-type representation of any order.
NON_PH_ALPHA = np.array([0.61058991931158258, -0.15547146730086722,
                         0.54488154798928464])
NON_PH_A = np.array([[-1.0, 0.0, 0.0],
                     [0.0, -2.0, 2.0],
                     [0.0, -2.0, -2.0]])

MARKER = "matrix-exponential or rational service process"


def _two_class_me_queue():
    model = Network("M/ME/1 two classes")
    source = Source(model, "Source")
    queue = Queue(model, "Queue", SchedStrategy.FCFS)
    sink = Sink(model, "Sink")
    class1 = OpenClass(model, "Class1", 0)
    class2 = OpenClass(model, "Class2", 0)
    source.setArrival(class1, Exp(0.2))
    source.setArrival(class2, Exp(0.2))
    queue.setService(class1, ME(NON_PH_ALPHA, NON_PH_A))
    queue.setService(class2, ME(NON_PH_ALPHA, NON_PH_A))
    P = model.initRoutingMatrix()
    P.set(class1, Network.serialRouting(source, queue, sink))
    P.set(class2, Network.serialRouting(source, queue, sink))
    model.link(P)
    return model


def _solve_capturing(num_solves):
    """Return the lines written to stdout while solving num_solves models."""
    buf = io.StringIO()
    from line_solver.api.io.logging import LineLogger
    logger = LineLogger.get_instance()
    saved_stdout = logger.stdout
    logger.stdout = buf
    try:
        with redirect_stdout(buf):
            for _ in range(num_solves):
                SolverMAM(_two_class_me_queue()).getAvgTable()
    finally:
        logger.stdout = saved_stdout
    return buf.getvalue()


def _count(output):
    return output.count(MARKER)


# The three tests below capture the native LineLogger stream while solving.
# Under lang='java' the solve is delegated to jline.jar over JSON, so the
# warning is emitted by the JAR (covered by MamMEWarningTest.java) and never
# reaches the python logger being captured here.
_skip_java_dispatch = pytest.mark.skipif(
    os.environ.get('LINE_SOLVER_LANG') == 'java',
    reason="captures the native LineLogger; under java dispatch the warning is "
           "emitted inside jline.jar, not by the python solver")


@_skip_java_dispatch
def test_me_multiclass_warning_is_raised():
    output = _solve_capturing(1)
    assert _count(output) == 1
    line = [ln for ln in output.splitlines() if MARKER in ln][0]
    assert line.startswith("Warning [solver_mam_basic]: ")
    # The station is named, not indexed: station indices are 1-based in MATLAB
    # and 0-based here and in the JAR, so an index would make the same message
    # read differently in each codebase.
    assert "Station Queue has a matrix-exponential or rational service process" in line
    assert "here 2 classes, 1 servers" in line
    assert "not exact for this service process" in line


@_skip_java_dispatch
def test_me_warning_raised_once_per_solve_not_per_iteration():
    # The fixed-point iteration visits the station several times per solve; the
    # user needs the limitation stated once, not once per sweep.
    assert _count(_solve_capturing(1)) == 1


@_skip_java_dispatch
def test_me_warning_repeats_across_models():
    # warnings.warn under the default filter would have shown this only for the
    # first model, and line_warning would have hidden it as a 60-second repeat,
    # leaving the second and third models looking clean.
    assert _count(_solve_capturing(3)) == 3


def test_line_warning_always_does_not_suppress_repeats():
    from line_solver.api.io.logging import LineLogger
    logger = LineLogger.get_instance()
    buf = io.StringIO()
    saved_stdout = logger.stdout
    logger.stdout = buf
    try:
        for _ in range(3):
            line_warning_always("test_mam_me_warning", "identical message %d", 7)
    finally:
        logger.stdout = saved_stdout
    assert buf.getvalue().count(
        "Warning [test_mam_me_warning]: identical message 7") == 3
