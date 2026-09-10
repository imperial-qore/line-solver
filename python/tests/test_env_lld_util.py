"""
Regression: the ENV state-vector analyzer must report the same utilization as the
plain CTMC for a load-dependent station.

_avg_from_pi (statevec.py), which mirrors MATLAB solver_ctmc_avg_from_pi.m, had
drifted from solver_ctmc_analyzer.m: in its load/class-dependent branch it
accumulated the bare per-class capacity share nir[k]*schedparam/(nir @ schedparam)
without weighting by the current scaling lldnow or normalizing by the effective
capacity ceff = max(nservers, max(lldscaling)). That is a P(busy)-style value, not
a busy-server fraction, and it overstated a load-dependent station's utilization
(0.709091 against the exact 0.490909 on the fixture below).

The oracle needs no external reference: an environment whose stages are IDENTICAL
cannot change anything, so it must reproduce the plain CTMC solution of that one
model exactly. Queue lengths already matched before the fix, which isolates the
defect to the utilization branch.
"""

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, ClosedClass, Exp, SchedStrategy,
                         Environment, SolverCTMC, SolverENV, SolverOptions)
from line_solver.constants import SolverType

TOL = 1e-9
N = 3
C = 2


def _build():
    """Closed model with a delay and a load-dependent PS queue, alpha(n)=min(n,c)."""
    m = Network('base')
    d = Delay(m, 'D')
    q = Queue(m, 'Q', SchedStrategy.PS)
    jc = ClosedClass(m, 'C', N, d)
    d.setService(jc, Exp(1.0))
    q.setService(jc, Exp(2.0))
    q.setLoadDependence(np.minimum(np.arange(1, N + 1), C))
    m.link(Network.serialRouting(d, q))
    return m


def test_env_statevec_matches_ctmc_on_identical_lld_stages():
    want = np.asarray(SolverCTMC(_build(), verbose=False).getAvgTable()['Util'], dtype=float)

    env = Environment('IdenticalStages')
    env.addStage(0, 'S1', 'operational', _build())
    env.addStage(1, 'S2', 'operational', _build())
    env.addTransition(0, 1, Exp(0.5))
    env.addTransition(1, 0, Exp(0.5))
    env.init()

    options = SolverOptions(SolverType.ENV)
    options.method = 'statevec'
    options.verbose = False
    # The state-vector analyzer under test is native: the delegating ENV engines
    # solve every stage by the fluid transient and refuse a CTMC-staged ensemble
    # outright, so the RECOMBINATION is pinned here while the stage solvers keep
    # the ambient lang and the CTMC oracle above is read the same way.
    options.lang = 'python'

    solver = SolverENV(env, lambda m: SolverCTMC(m, timespan=[0, 1e6], verbose=False),
                       options=options)
    got = np.asarray(solver.getAvgTable()['Util'], dtype=float)

    # the environment is a no-op, so every station's utilization must agree exactly
    np.testing.assert_allclose(got, want, atol=TOL)


def test_env_statevec_matches_ctmc_on_queue_lengths_too():
    """Queue lengths were already correct; pin them so the fix cannot disturb them."""
    want = np.asarray(SolverCTMC(_build(), verbose=False).getAvgTable()['QLen'], dtype=float)

    env = Environment('IdenticalStages')
    env.addStage(0, 'S1', 'operational', _build())
    env.addStage(1, 'S2', 'operational', _build())
    env.addTransition(0, 1, Exp(0.5))
    env.addTransition(1, 0, Exp(0.5))
    env.init()

    options = SolverOptions(SolverType.ENV)
    options.method = 'statevec'
    options.verbose = False
    # The state-vector analyzer under test is native: the delegating ENV engines
    # solve every stage by the fluid transient and refuse a CTMC-staged ensemble
    # outright, so the RECOMBINATION is pinned here while the stage solvers keep
    # the ambient lang and the CTMC oracle above is read the same way.
    options.lang = 'python'

    solver = SolverENV(env, lambda m: SolverCTMC(m, timespan=[0, 1e6], verbose=False),
                       options=options)
    got = np.asarray(solver.getAvgTable()['QLen'], dtype=float)
    np.testing.assert_allclose(got, want, atol=TOL)
