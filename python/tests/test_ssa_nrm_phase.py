"""
Native-Python tests for the NRM phase expansion at INF/PS-family stations.

Mirror of the JAR ``SolverSSANrmPhaseTest``. Non-exponential service at an
INF/PS-family station is expanded over service phases rather than collapsed
onto its mean rate; these cases assert the expansion reproduces the exact
insensitive means.

The fixture is a closed two-job network whose station under test has mean
service 0.5 and whose Delay has mean service 1.0. PS and INF are both
insensitive stations under BCMP, so the exact means depend on the service
distribution only through its mean:

    PS  station: Q = 0.8,   U = 0.6,   X = 1.2
    INF station: Q = 2/3,   U = 2/3,   X = 4/3

Insensitivity is what makes this a usable regression: every distribution below
has mean 0.5, so all of them must return the SAME exact value, and no CTMC
solve is needed to obtain it. The Exp case is the control -- with nph == 1 the
slot map degenerates to the pre-expansion flat class index, so Exp must
reproduce its pre-expansion value and any Exp deviation indicts the layout
rather than the phase logic.

These tests exist because the expansion was previously covered only by
all-exponential models, in which the dependency graph's arithmetic decode of a
state index is still correct and its defects are therefore invisible: the
failure mode is a quiet bias that GROWS with distance from exponential (see
_kb/log.md [2026-07-17]). The Erlang/HyperExp/Coxian cases carry the signal; a
green all-exponential suite proves nothing about them.

Tolerance: at SAMPLES the single-seed noise floor on Q is sd ~= 0.002 (0.2%),
calibrated by an 8-seed spread on the Exp control. REL_TOL is set to 1% -- ~5
sd above that floor, so the fixed-seed assertion is not flaky, yet below the
signature of the decode defect this covers (Erlang-3 +1.3%, HyperExp +4.0%),
so a recurrence fails rather than passing inside a loose band.
"""

import numpy as np
import pytest

from line_solver import (
    Network, Delay, Queue, Exp, Erlang, HyperExp, Coxian,
    ClosedClass, SchedStrategy, SolverSSA,
)

SAMPLES = 500000
SEED = 23000
REL_TOL = 0.01
MEAN = 0.5


def _dists():
    return [
        ('Exp', Exp(1.0 / MEAN)),
        ('Erlang-2', Erlang.fitMeanAndOrder(MEAN, 2)),
        ('Erlang-3', Erlang.fitMeanAndOrder(MEAN, 3)),
        ('HyperExp-4', HyperExp.fitMeanAndSCV(MEAN, 4.0)),
        ('Coxian-3', Coxian.fitMeanAndSCV(MEAN, 3.0)),
    ]


def _build(sched, dist):
    model = Network('phase')
    queue = Queue(model, 'Q1', sched)
    delay = Delay(model, 'Delay')
    cc = ClosedClass(model, 'C1', 2, queue)
    queue.setService(cc, dist)
    delay.setService(cc, Exp(1.0))
    model.link(Network.serialRouting(queue, delay))
    return model, queue


def _nrm_qlen(model, queue):
    solver = SolverSSA(model, 'nrm', seed=SEED, samples=SAMPLES, verbose=False)
    return float(np.atleast_2d(solver.getAvgQLen())[queue.getStationIndex() - 1, 0])


@pytest.mark.parametrize('name,dist', _dists(), ids=[d[0] for d in _dists()])
@pytest.mark.parametrize('sched,exact,label', [
    (SchedStrategy.PS, 0.8, 'PS'),
    (SchedStrategy.INF, 2.0 / 3.0, 'INF'),
], ids=['PS', 'INF'])
def test_phase_expansion_insensitive(sched, exact, label, name, dist):
    """Every mean-0.5 distribution must return the same exact insensitive mean."""
    model, queue = _build(sched, dist)
    got = _nrm_qlen(model, queue)
    err = abs(got - exact) / exact
    assert err < REL_TOL, (
        f'{label}/{name}: NRM QLen {got:.4f} deviates from the exact '
        f'insensitive value {exact:.4f} by {100.0 * err:.2f}%')
