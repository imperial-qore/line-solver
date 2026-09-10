"""The peak rate that normalizes utilization is a MODEL INPUT, not a default.

Utilization at a rate-dependent station is reported as ``U = T*E[S]/peak``, the
same fraction-of-capacity as the ``T*S/c`` of an ordinary multiserver station.
Load dependence always knows its peak, ``max(c, max alpha)``, because the user
supplied every alpha.  CLASS and JOINT dependence do not: the scaling is a
CALLABLE, so recovering ``max_n beta(n)`` would mean sweeping the population
lattice -- which needs a bound the callable does not carry, and an open class
has no bound at all.  So the peak must be declared, and a declaration without
one is refused rather than guessed at.

The refusal lands on the SETTER here, as it does in MATLAB, so a peak-less
model can never reach a solver.  The C++ port instead accepts an empty peak and
refuses at every reader; both contracts end in an error, which is the property
these tests pin.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, Network, OpenClass, Queue,
                         SchedStrategy, Sink, Source)


def _beta(n):
    """min(n_1, 2): a callable whose lattice peak is 2."""
    return np.array([min(float(np.atleast_1d(n)[0]), 2.0)])


def _closed_model():
    model = Network('cdpeak')
    delay = Delay(model, 'Delay')
    queue = Queue(model, 'Queue', SchedStrategy.PS)
    cclass = ClosedClass(model, 'Class1', 3, delay)
    delay.setService(cclass, Exp(1.0))
    queue.setService(cclass, Exp(1.0))
    model.link(model.serialRouting(delay, queue))
    return model, queue


def test_class_dependence_without_peak_is_refused():
    _, queue = _closed_model()
    with pytest.raises(ValueError, match='peak'):
        queue.set_class_dependence(_beta)


def test_joint_dependence_without_peak_is_refused():
    _, queue = _closed_model()
    with pytest.raises(ValueError, match='peak'):
        queue.set_joint_dependence(_beta)


@pytest.mark.parametrize('bad', [0.0, -2.0])
def test_non_positive_peak_is_refused(bad):
    # a peak of zero would divide, and a negative one would flip the sign
    _, queue = _closed_model()
    with pytest.raises(ValueError, match='positive'):
        queue.set_class_dependence(_beta, bad)
    with pytest.raises(ValueError, match='positive'):
        queue.set_joint_dependence(_beta, bad)


def test_declared_peak_reaches_the_struct_unchanged():
    # beta maxes at 2 over the lattice but is DECLARED as 5, and the struct must
    # carry 5: the declaration is the normalizer, never a sweep of the callable.
    model, queue = _closed_model()
    queue.set_class_dependence(_beta, 5.0)
    sn = model.getStruct()
    peak = np.asarray(sn.cdscalingpeak)
    assert np.nanmax(peak) == pytest.approx(5.0)


def test_open_class_callable_has_no_lattice_to_sweep():
    # THE CASE THAT MOTIVATES THE ERROR: beta(n) = 1 + n is unbounded and the
    # population is infinite, so there is nothing to sweep. A derived peak here
    # would just be beta(0) = 1.
    model = Network('openCd')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.PS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1')
    source.setArrival(oclass, Exp(0.5))
    queue.setService(oclass, Exp(1.0))
    model.link(model.serialRouting(source, queue, sink))

    def growing(n):
        return np.array([1.0 + float(np.atleast_1d(n)[0])])

    with pytest.raises(ValueError, match='peak'):
        queue.set_class_dependence(growing)
    queue.set_class_dependence(growing, 11.0)
    assert np.nanmax(np.asarray(model.getStruct().cdscalingpeak)) == pytest.approx(11.0)


def test_the_peak_reaches_the_ldes_engine_over_the_wire():
    # MATLAB and Python reach the LDES engine only through model.json, so a peak
    # the writer drops is a peak the engine cannot enforce -- it would refuse a
    # model the user declared correctly, or normalize by the wrong number. The
    # engine's own numbers are pinned in cpp/tests/test_dependent_peak_required.
    import json
    import os
    import tempfile

    from line_solver.io.linemodel_io import load_model, save_model

    model, queue = _closed_model()
    queue.set_class_dependence(_beta, 2.0)

    fd, path = tempfile.mkstemp(suffix='.json')
    os.close(fd)
    try:
        save_model(model, path)
        with open(path) as fh:
            doc = json.load(fh)
        nodes = {n['name']: n for n in doc['model']['nodes']}
        blk = nodes['Queue']['classDependence']
        assert 'peak' in blk, 'the declared peak did not reach the wire'
        assert max(blk['peak']) == pytest.approx(2.0)
        # and it survives the round trip rather than being re-swept from the table
        back = load_model(path)
        assert np.nanmax(np.asarray(back.getStruct().cdscalingpeak)) == pytest.approx(2.0)
    finally:
        os.unlink(path)


def test_global_dependence_also_requires_its_peak():
    # the Whittle primitive reads the whole population matrix, so its peak is
    # even less recoverable; every codebase makes it mandatory
    model, _ = _closed_model()
    with pytest.raises((ValueError, TypeError)):
        model.set_global_dependence(lambda n: np.ones(1))
