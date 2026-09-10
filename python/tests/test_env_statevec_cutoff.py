"""
Regression: the ENV state-vector analyzer must be CUTOFF-INVARIANT when one
environment stage is individually unstable.

`pre` seeded every stage's entry distribution with THAT STAGE'S OWN stationary
law, `ctmc_solve_reducible(Q_e)`. A stage whose arrival rate is at or above its
own service rate has no stationary law, so the reducible solve returned the
stationary law of the TRUNCATED generator -- mass piled against the truncation
wall, mean growing linearly with the cutoff (uniform, mean N/2, for a critical
M/M/1 truncated at N).

The fixed point itself was never wrong: chaining the resolvent is exactly the
stationary equation of the joint (queue, stage) chain, so it contracts to the
right answer. What grew with the cutoff was the number of sweeps needed to drain
the seed, which made the reported number drift AWAY from the truth as the cutoff
was RAISED -- the natural response to a suspect answer made it worse, and a small
cutoff looking right was a coincidence.

Here the Slow stage is critical (lambda = mu = 0.8) while the system is stable on
average (mu_bar = 1.4), the regime that exposed it. Both cutoffs sit far past the
real mass, so anything that moves between them is the seed and not the model.
MATLAB, before the fix and on the same sweep budget, gave 1.478346 at cutoff 40
and 28.608647 at cutoff 300 against an exact 1.478318.

THE SWEEP BUDGET IS PART OF THE TEST. The fixed seed reaches the fixed point in
~220 sweeps at EVERY cutoff; the old per-stage seed needs ~350 at cutoff 40 but
~650 at cutoff 150, which is why the LARGER cutoff is what discriminates here.
Raising iter_max far enough lets the old seed converge too and SILENTLY DEFANGS
this test -- it would still pass while testing nothing.

Also covers the alias: 'blend' is 'statevec' in MATLAB, the JAR and C++, and
python listed it without wiring it, so a 'blend' run silently fell through to the
MEAN-FIELD analyzer -- a different approximation, which returned 0.948301 here.
"""

import warnings

import numpy as np
import pytest

from line_solver import (Network, Source, Queue, Sink, OpenClass, Exp,
                         SchedStrategy, Environment, SolverCTMC, SolverENV)

#: Exact joint (queue, stage) CTMC answer for this environment.
EXACT_Q = 1.478318
LAMBDA = 0.8


def _mm1(lam, mu):
    """Open M/M/1: Source -> FCFS Queue -> Sink, one class."""
    m = Network('mm1')
    src = Source(m, 'Source')
    q = Queue(m, 'Q', SchedStrategy.FCFS)
    snk = Sink(m, 'Sink')
    c = OpenClass(m, 'C')
    src.setArrival(c, Exp(lam))
    q.setService(c, Exp(mu))
    m.link(Network.serialRouting(src, q, snk))
    return m


def _env():
    e = Environment('mode')
    e.addStage('Slow', 'degraded', _mm1(LAMBDA, 0.8))       # critical on its own
    e.addStage('Fast', 'operational', _mm1(LAMBDA, 2.0))
    e.addTransition('Slow', 'Fast', Exp(1.0))
    e.addTransition('Fast', 'Slow', Exp(1.0))
    e.init()
    return e


def _solve(cutoff, method='statevec', iter_max=400):
    """Queue-station QLen and throughput, plus any non-convergence warnings."""
    factory = lambda m: SolverCTMC(m, 'exact', timespan=[0, 200],
                                   cutoff=cutoff, verbose=False)
    options = {'iter_max': iter_max, 'iter_tol': 1e-10,
               'verbose': False, 'method': method}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        QN, _, _, TN = SolverENV(_env(), factory, options).getAvg()[:4]
        noconv = [str(w.message) for w in caught
                  if 'did not converge' in str(w.message)]
    # Station 1 is the Queue; station 0 is the Source, whose QLen is inf by
    # convention and must not be summed into the total.
    return float(np.asarray(QN).ravel()[1]), float(np.asarray(TN).ravel()[1]), noconv


@pytest.mark.parametrize('cutoff', [40, 150])
def test_critical_stage_is_cutoff_invariant(cutoff):
    qlen, tput, noconv = _solve(cutoff)
    # A run that stopped on iter_max says nothing about the fixed point, so
    # assert convergence first -- that is the half that used to be silent.
    assert not noconv, 'the fixed point must converge within iter_max'
    # Flow balance at the queue is the cheap tell the defect violated by up to
    # 33%: in steady state the queue must clear exactly what arrives.
    assert tput == pytest.approx(LAMBDA, abs=1e-6)
    assert qlen == pytest.approx(EXACT_Q, abs=1e-4)


def test_raising_the_cutoff_does_not_move_the_answer():
    small, _, _ = _solve(40)
    large, _, _ = _solve(150)
    assert small == pytest.approx(large, abs=1e-6)


def test_blend_is_an_alias_for_statevec():
    statevec = _solve(40, method='statevec')
    blend = _solve(40, method='blend')
    assert blend[0] == pytest.approx(statevec[0], abs=1e-9)
    assert blend[1] == pytest.approx(statevec[1], abs=1e-9)


def test_non_convergence_is_reported():
    """Exhausting iter_max must warn, not return the last iterate silently."""
    _, _, noconv = _solve(60, iter_max=2)
    assert noconv, 'a run that stopped on iter_max must warn'
