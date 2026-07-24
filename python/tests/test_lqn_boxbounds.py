"""
Validates the LQN Majumdar-Woodside robust box bounds (SolverLN mwba.upper /
mwba.lower) on a two-tier client-server LQN.

Reference task T1 (mult 10, think 5.0) on processor P1 runs a single activity
(host demand 1.0) that makes 2.5 synchronous calls to entry E2 on an infinite
DB task on processor P2 (host demand 0.8). Per-chain processor demands:
D_P1 = 1.0, D_P2 = 2.5*0.8 = 2.0, Z = 5, N = 10. The processor-utilization
upper bound gives X_ref+ = 1/D_P2 = 0.5 (P2 bottleneck); the Majumdar-Woodside
lower bound gives X_ref- = N/(Z + N*(D_P1+D_P2)) = 10/35 = 0.2857. The exact LN
throughput (~0.499) lies within [0.2857, 0.5].
"""

import numpy as np
import pytest

from line_solver import (LayeredNetwork, Processor, Task, Entry, Activity, Exp,
                          SchedStrategy, SolverLN, SolverMVA)

TOL = 1e-4


def _build():
    m = LayeredNetwork('cd')
    P1 = Processor(m, 'P1', 1, SchedStrategy.PS)
    P2 = Processor(m, 'P2', 1, SchedStrategy.PS)
    T1 = Task(m, 'T1', 10, SchedStrategy.REF).on(P1)
    T1.set_think_time(Exp.fit_mean(5.0))
    T2 = Task(m, 'T2', float('inf'), SchedStrategy.INF).on(P2)
    E1 = Entry(m, 'E1').on(T1)
    E2 = Entry(m, 'E2').on(T2)
    Activity(m, 'A1', Exp.fit_mean(1.0)).on(T1).bound_to(E1).synch_call(E2, 2.5)
    Activity(m, 'A2', Exp.fit_mean(0.8)).on(T2).bound_to(E2).replies_to(E2)
    return m


def _tput(df, node):
    names = list(np.asarray(df['Node']))
    tput = list(np.asarray(df['Tput']))
    return float(tput[names.index(node)])


def test_lqn_boxbounds_bracket_exact():
    up = SolverLN(_build(), lambda mm: SolverMVA(mm), method='mwba.upper').getAvgTable()
    lo = SolverLN(_build(), lambda mm: SolverMVA(mm), method='mwba.lower').getAvgTable()
    ex = SolverLN(_build(), lambda mm: SolverMVA(mm)).getAvgTable()

    x_up = _tput(up, 'T1')
    x_lo = _tput(lo, 'T1')
    x_ex = _tput(ex, 'T1')

    assert x_up == pytest.approx(0.5, abs=TOL)
    assert x_lo == pytest.approx(10.0 / 35.0, abs=TOL)
    assert x_lo <= x_ex + TOL <= x_up + 2 * TOL

    # called-entry throughput scales with the 2.5 call multiplicity
    assert _tput(up, 'E2') == pytest.approx(2.5 * 0.5, abs=TOL)
    assert _tput(lo, 'E2') == pytest.approx(2.5 * 10.0 / 35.0, abs=TOL)
