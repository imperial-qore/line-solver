"""
Native-Python tests for the NRM's DROP-rule finite capacity regions, with the
region on an INTERIOR station.

Mirror of the JAR ``SolverSSANrmFcrDropTest``.

The interior placement is the whole point of this fixture. Under DROP the
refused job is DESTROYED: the departure fires at its full rate and the job never
reaches the destination. A previous implementation instead CENSORED the refused
transition by scaling the propensity with the admitted share, which holds the
job at its SOURCE. The two are observationally identical when the region is fed
by a Source -- a Source's population is fictitious, so "destroy the arriving
job" and "censor the arrival" cannot be told apart -- and every pre-existing FCR
fixture placed the region exactly there. They diverge qualitatively as soon as
the blocked job would be leaving a real queue: the upstream queue then grows
without bound and nothing is ever lost.

Fixture: Source(1.0) -> Q1/PS(3.0) -> Q2/PS(1.5) -> Sink, region over {Q2},
globalMaxJobs = 1, DropStrategy.DROP. The exact CTMC (cutoff-invariant from 16
upward) gives Q1 = 0.5 and T_Q2 = 0.6: 40% of jobs are lost at the boundary and
Q1 is stable. The censoring implementation returned a non-stationary Q1 (108.7 /
112.1 / 303.6 at 100k / 400k / 1.6M samples) and T_Q2 -> 0.995.

The exact values are hard-coded rather than obtained from a CTMC solve so the
assertion cannot drift with the reference solver, and because Q1 = 0.5 is
analytic: with the boundary loss Q1 is an M/M/1-PS at utilisation 1/3.
"""

import numpy as np
import pytest

from line_solver import (
    Network, Source, Sink, Queue, Exp, OpenClass,
    SchedStrategy, DropStrategy, SolverSSA,
)

SAMPLES = 200000
SEEDS = [23000, 24000, 25000, 26000, 27000]
REL_TOL = 0.02


def _build(interior, cap):
    model = Network('fcr_drop')
    source = Source(model, 'Source')
    q1 = Queue(model, 'Q1', SchedStrategy.PS)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    sink = Sink(model, 'Sink')
    c1 = OpenClass(model, 'C1')
    source.setArrival(c1, Exp(1.0))
    q1.setService(c1, Exp(3.0))
    q2.setService(c1, Exp(1.5))
    model.link(Network.serialRouting(source, q1, q2, sink))
    region = model.addRegion([q2 if interior else q1])
    region.setGlobalMaxJobs(cap)
    region.setDropRule(c1, DropStrategy.DROP)
    return model


def _mean_q1_and_t2(interior, cap):
    """Mean over SEEDS of (QLen at Q1, Tput at Q2). Station rows: 1 = Q1, 2 = Q2."""
    q1s, t2s = [], []
    for seed in SEEDS:
        solver = SolverSSA(_build(interior, cap), 'nrm',
                           seed=seed, samples=SAMPLES, verbose=False)
        q = np.atleast_2d(solver.getAvgQLen())
        t = np.atleast_2d(solver.getAvgTput())
        # The NRM must actually have run. The handler downgrades an ineligible
        # model to the serial engine SILENTLY, and the serial engine is CORRECT
        # here -- so a fallback would make this test pass while never exercising
        # the code it exists to cover.
        ran = getattr(solver._result, 'method', None)
        assert ran == 'nrm', f'NRM did not run: method was {ran!r}'
        q1s.append(float(q[1, 0]))
        t2s.append(float(t[2, 0]))
    return float(np.mean(q1s)), float(np.mean(t2s))


@pytest.mark.parametrize('interior,cap,exact_q1,exact_t2,label', [
    (True, 1, 0.5, 0.6, 'interior'),
    (False, 2, 0.38462, 0.92307, 'source-fed'),
], ids=['interior', 'source-fed'])
def test_fcr_drop(interior, cap, exact_q1, exact_t2, label):
    """DROP must destroy the refused job.

    ``interior``: the case the censoring build failed with a divergent Q1.
    ``source-fed``: the pre-existing configuration, which must not regress.
    """
    got_q1, got_t2 = _mean_q1_and_t2(interior, cap)
    for got, exact, what in ((got_q1, exact_q1, 'Q1'), (got_t2, exact_t2, 'T_Q2')):
        err = abs(got - exact) / exact
        assert err < REL_TOL, (
            f'{label} region {what}: NRM returned {got:.5f}, exact is {exact} '
            f'({100.0 * err:.2f}% off)')
