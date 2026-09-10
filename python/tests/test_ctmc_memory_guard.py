"""Two-ended coverage for the CTMC state-space estimator behind the memory gate.

The estimator had no test in any codebase, which is how two defects survived in
it. Both were OVER-estimates, and an over-estimate is the failure mode that hides:
it refuses a model that solves in seconds, and the refusal reads as a legitimate
capacity limit rather than as a bug.

Every test therefore pins BOTH ends. A test that only asserts "this is admitted"
would pass over an estimator that admits everything, and the gate exists precisely
to refuse. The refusal cases are as load-bearing as the admissions.
"""

import math

import pytest

from line_solver import (ClosedClass, Delay, Disabled, Exp, Network, OpenClass,
                         Queue, SchedStrategy, Sink, Source)
from line_solver.api.solvers.ctmc.memory_guard import (ctmc_memory_gate,
                                                       state_space_log_size)


class _Options:
    def __init__(self, cutoff=None):
        self.cutoff = cutoff


def _nstates(model, cutoff=None):
    return math.exp(state_space_log_size(model.getStruct(), _Options(cutoff)))


def _admits(model, cutoff=None):
    ok, _ = ctmc_memory_gate(state_space_log_size(model.getStruct(), _Options(cutoff)))
    return ok


def _closed_ps(populations):
    """Delay + PS queue, one class per entry of `populations`, served everywhere."""
    model = Network('closed_ps')
    think = Delay(model, 'Think')
    queue = Queue(model, 'Queue', SchedStrategy.PS)
    for r, n in enumerate(populations):
        jobclass = ClosedClass(model, 'Class%d' % (r + 1), n, think, 0)
        think.setService(jobclass, Exp(1.0))
        queue.setService(jobclass, Exp(2.0))
    model.link(Network.serialRouting(think, queue))
    return model


def _closed_fcfs(populations):
    """As `_closed_ps` but the queue keeps a class SEQUENCE, so states multiply."""
    model = Network('closed_fcfs')
    think = Delay(model, 'Think')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    for r, n in enumerate(populations):
        jobclass = ClosedClass(model, 'Class%d' % (r + 1), n, think, 0)
        think.setService(jobclass, Exp(1.0))
        queue.setService(jobclass, Exp(2.0))
    model.link(Network.serialRouting(think, queue))
    return model


def _disjoint_routes(nroutes, cutoff):
    """One PS route per open class, each DISABLED for every other class.

    This is the shape of examples/advanced/loadDependent/ld_whittle_bandwidth:
    a class can only ever be at its own route, so the true space is
    (cutoff+1)**nroutes rather than anything spread over all the stations.
    """
    model = Network('disjoint')
    source = Source(model, 'Source')
    sink = Sink(model, 'Sink')
    routes = [Queue(model, 'Route%d' % (s + 1), SchedStrategy.PS) for s in range(nroutes)]
    classes = [OpenClass(model, 'Flow%d' % (s + 1), 0) for s in range(nroutes)]
    for s in range(nroutes):
        source.setArrival(classes[s], Exp(0.2))
        for t in range(nroutes):
            if t == s:
                routes[t].setService(classes[s], Exp(1.0))
            else:
                routes[t].setService(classes[s], Disabled())
    routing = model.initRoutingMatrix()
    for s in range(nroutes):
        routing.set(classes[s], classes[s], source, routes[s], 1.0)
        routing.set(classes[s], classes[s], routes[s], sink, 1.0)
    model.link(routing)
    return model


@pytest.mark.parametrize('populations,expected', [
    ([4], 5),
    ([3, 3], 16),
    ([2, 2, 2], 27),
])
def test_share_stations_are_counted_exactly(populations, expected):
    """With no ordered buffer the count is stars and bars, and it is EXACT.

    Pinning the value rather than an inequality is deliberate: a bound would be
    satisfied by an estimator that has drifted by orders of magnitude, which is
    the only drift that matters here.
    """
    assert _nstates(_closed_ps(populations)) == pytest.approx(expected, rel=1e-9)


def test_a_class_disabled_at_a_station_is_not_placed_there():
    """The placement term must spread a class over the stations that ADMIT it.

    Counting all M stations priced ld_whittle_bandwidth at C(9,6)**3 = 592704
    states -- 7852 GB under the quadratic byte model -- against a true 7**3 = 343,
    and the gate refused a model that solves in under a minute.
    """
    cutoff = 6
    model = _disjoint_routes(3, cutoff)
    assert _nstates(model, cutoff) == pytest.approx((cutoff + 1) ** 3, rel=1e-9)
    assert _admits(model, cutoff)


def test_the_exact_dp_is_reached_when_its_own_cost_is_small():
    """The DP grid is prod_k (n_k+1), NOT the per-station box over every station.

    The second proxy is exponential in the station count: on five ordered stations
    holding two classes of three it read 16**5 and tripped the 1e6 guard, sending a
    DP costing some 1280 steps to the geometric fallback, which priced
    mqn_multiserver_fcfs at 2.4e13 GB. The three other codebases never used it.
    """
    model = Network('five_fcfs')
    queues = [Queue(model, 'Queue%d' % (i + 1), SchedStrategy.FCFS) for i in range(5)]
    for i, queue in enumerate(queues):
        queue.setNumberOfServers(i + 1)
    closed = ClosedClass(model, 'ClosedClass', 3, queues[0], 0)
    for queue in queues:
        queue.setService(closed, Exp(1.0))
    model.link(Network.serialRouting(*queues))
    # The geometric fallback would price this past any budget; the exact DP does not.
    assert _admits(model, 3)


def test_an_ordered_buffer_still_multiplies_by_its_sequences():
    """FCFS must cost strictly more than PS at the same population.

    Omitting the sequence factor UNDER-estimates, and an under-estimate is the
    direction that reaches an out-of-memory rather than a refusal: one omission
    priced a model at 726 states against 225840 enumerated and the solve then
    attempted a dense (225840, 225840) array.
    """
    assert _nstates(_closed_fcfs([3, 3])) > _nstates(_closed_ps([3, 3]))


def test_a_genuinely_large_model_is_still_refused():
    """The admitting end is only meaningful if the refusing end still holds."""
    model = _closed_fcfs([10, 10, 10])
    assert _nstates(model) > 1e12
    assert not _admits(model)
