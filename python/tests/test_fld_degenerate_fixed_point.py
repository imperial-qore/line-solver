"""The first-order fluid drift on a model whose fixed point is NOT ISOLATED.

`min(E[n], c)` is flat above the server count: a saturated station delivers its
full capacity whatever its queue holds, so the drift cannot tell one split of
the mass between two saturated stations from another. Every such split is an
equilibrium, and the raw first-order method returns whichever one the
integrator stopped at -- on two identical saturated stations in a closed cycle
it came back with [9 1] where the exact chain says [5 5], and [8 2] on the same
pair with two servers each. Population was conserved; the answer was simply on
the wrong point of the surface.

The repair replaces `min(E[n], c)` by `E[min(n, c)]` under the station's
equilibrium geometric marginal, which is strictly increasing and therefore
isolates one fixed point, and it RE-INTEGRATES the drift rather than adjusting
the point that came back. It fires only where the probe finds a null direction,
so a well-posed model is untouched -- which is what the second half of this
file asserts, and it is the half that would catch an over-eager repair.

Why a closure and not a smoothed min: any smoothing sharp enough to stay
faithful to min away from the kink is numerically flat far from it. The
Boltzmann softmin at alpha = 20 carries a restoring force of exp(-160) at the
[9 1] point, and the p-norm trades faithfulness against selection directly --
pstar = 2 recovers [5 5], pstar = 8 gives [7.64 2.36], pstar = 128 gives
[8.94 1.06]. The closure has no such trade-off: its slope comes from the
VARIANCE of the marginal, not from a smoothing width.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, FLD, CTMC, Network, OpenClass,
                         Queue, SchedStrategy, Sink, Source, GlobalConstants,
                         VerboseLevel)

GlobalConstants.setVerbose(VerboseLevel.SILENT)


def _identical_pair(sched, N, servers=1):
    """Two identical stations on a closed cycle, saturated by construction."""
    m = Network('degenerate')
    q1 = Queue(m, 'Q1', sched)
    q2 = Queue(m, 'Q2', sched)
    q1.setNumberOfServers(servers)
    q2.setNumberOfServers(servers)
    c = ClosedClass(m, 'C', N, q1)
    q1.setService(c, Exp(1.0))
    q2.setService(c, Exp(1.0))
    m.link(Network.serialRouting(q1, q2))
    return m


@pytest.mark.parametrize('sched,N,servers', [
    (SchedStrategy.PS, 10, 1),
    (SchedStrategy.FCFS, 6, 1),
    (SchedStrategy.PS, 10, 2),
])
def test_a_degenerate_pair_lands_on_the_symmetric_point(sched, N, servers):
    """The continuum is resolved to the point the exact chain agrees with."""
    QN = FLD(_identical_pair(sched, N, servers), method='matrix').avg_table()['QLen'].values
    assert np.allclose(QN, [N / 2.0, N / 2.0], atol=1e-4), \
        'expected the symmetric split, got %s' % (list(QN),)


@pytest.mark.parametrize('sched,N,servers', [
    (SchedStrategy.PS, 10, 1),
    (SchedStrategy.FCFS, 6, 1),
])
def test_the_repaired_point_matches_the_exact_chain(sched, N, servers):
    """Not merely symmetric: the same answer the CTMC gives."""
    QN = FLD(_identical_pair(sched, N, servers), method='matrix').avg_table()['QLen'].values
    QX = CTMC(_identical_pair(sched, N, servers)).avg_table()['QLen'].values
    assert np.allclose(QN, QX, atol=1e-4)


def _asymmetric_pair(sched, N, d1, d2):
    m = Network('asym')
    q1 = Queue(m, 'Q1', sched)
    q2 = Queue(m, 'Q2', sched)
    c = ClosedClass(m, 'C', N, q1)
    q1.setService(c, Exp(1.0 / d1))
    q2.setService(c, Exp(1.0 / d2))
    m.link(Network.serialRouting(q1, q2))
    return m


def _delay_and_queue(sched, N, Z, d):
    m = Network('dz')
    de = Delay(m, 'Delay')
    q = Queue(m, 'Q', sched)
    c = ClosedClass(m, 'C', N, de)
    de.setService(c, Exp(1.0 / Z))
    q.setService(c, Exp(1.0 / d))
    m.link(Network.serialRouting(de, q))
    return m


def _open_tandem(sched, lam, d1, d2):
    m = Network('oqn')
    s = Source(m, 'Source')
    q1 = Queue(m, 'Q1', sched)
    q2 = Queue(m, 'Q2', sched)
    k = Sink(m, 'Sink')
    c = OpenClass(m, 'C')
    s.setArrival(c, Exp(lam))
    q1.setService(c, Exp(1.0 / d1))
    q2.setService(c, Exp(1.0 / d2))
    m.link(Network.serialRouting(s, q1, q2, k))
    return m


# The values on the right are what `matrix` returned BEFORE the repair existed.
# They are here to assert the repair is inert on a well-posed model, so they are
# not a golden of the fluid approximation itself -- only that it did not move.
@pytest.mark.parametrize('build,expected', [
    (lambda: _asymmetric_pair(SchedStrategy.PS, 10, 1.0, 2.0), [0.5, 9.5]),
    (lambda: _asymmetric_pair(SchedStrategy.FCFS, 6, 1.0, 2.0), [0.5, 5.5]),
    (lambda: _delay_and_queue(SchedStrategy.PS, 10, 1.0, 0.5), [2.0, 8.0]),
    (lambda: _delay_and_queue(SchedStrategy.FCFS, 6, 1.0, 0.5), [2.0, 4.0]),
    (lambda: _open_tandem(SchedStrategy.PS, 0.5, 1.0, 1.0), [0.0, 0.5, 0.5]),
    (lambda: _open_tandem(SchedStrategy.FCFS, 0.5, 1.0, 1.0), [0.0, 0.5, 0.5]),
])
def test_a_well_posed_model_is_untouched(build, expected):
    """An isolated fixed point must not be re-integrated, let alone moved."""
    QN = FLD(build(), method='matrix').avg_table()['QLen'].values
    assert np.allclose(QN, expected, atol=1e-4), \
        'the repair moved a non-degenerate model: %s' % (list(QN),)


def test_a_delay_alone_does_not_look_degenerate():
    """An INF station carries theta = x and no min() at all.

    It is excluded from the probe: it can never be the pinned coordinate, and
    including it would make every Delay+Queue model read as degenerate.
    """
    QN = FLD(_delay_and_queue(SchedStrategy.PS, 10, 1.0, 0.5),
             method='matrix').avg_table()['QLen'].values
    assert np.allclose(QN, [2.0, 8.0], atol=1e-4)
