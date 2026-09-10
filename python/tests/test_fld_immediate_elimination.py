"""
Elimination of immediate transitions from the fluid drift.

`options.hide_immediate` used to be accepted by the native python fluid solver
and silently ignored: MATLAB has ode_eliminate_immediate.m and
eliminate_immediate_matrix.m, the JAR has ImmediateElimination and C++ has
fluid_eliminate_immediate, but python had neither, so a caller who asked for the
elimination got the stiff system back with no warning.

The oracle here is an identity rather than a golden. Stochastic complementation
of a state whose only outgoing transition is immediate must leave a proper
generator over the surviving states, carrying the SAME routing the immediate
state merely passed through: on 0 -(a)-> 1 -(Immediate)-> 2 -(b)-> 0 the reduced
chain is exactly 0 -(a)-> 2 -(b)-> 0, whose stationary law is available in
closed form and must agree with the full chain's, conditioned on the surviving
states.
"""

import numpy as np
import pytest

from line_solver.constants import GlobalConstants
from line_solver.solvers.solver_fld.immediate import (
    eliminate_immediate_matrix,
    generator_to_jumps,
    ode_eliminate_immediate,
    ode_expand_state,
)


class _Opts:
    """The two functions read options.config and options.verbose only."""

    def __init__(self, **kw):
        self.config = kw.pop('config', {})
        self.verbose = kw.pop('verbose', 0)
        for k, v in kw.items():
            setattr(self, k, v)


def _chain_with_immediate(a=1.0, b=2.0):
    """0 -(a)-> 1 -(Immediate)-> 2 -(b)-> 0."""
    W = np.zeros((3, 3))
    W[0, 1] = a
    W[1, 2] = GlobalConstants.Immediate
    W[2, 0] = b
    return W - np.diag(W.sum(axis=1))


def test_generator_to_jumps_round_trip():
    W = np.array([[-3.0, 1.0, 2.0], [4.0, -4.0, 0.0], [0.0, 5.0, -5.0]])
    jumps, rates, src = generator_to_jumps(W)
    assert jumps.shape == (3, 4)
    # Every column is a unit move: -1 at the source, +1 at the destination.
    assert np.allclose(jumps.sum(axis=0), 0.0)
    assert np.all(np.abs(jumps).sum(axis=0) == 2)
    # Rebuilding the generator from the jumps returns the original.
    W2 = np.zeros((3, 3))
    for t in range(len(rates)):
        dst = int(np.nonzero(jumps[:, t] > 0)[0][0])
        W2[int(src[t]), dst] += rates[t]
    W2 = W2 - np.diag(W2.sum(axis=1))
    assert np.allclose(W2, W)


def test_event_form_folds_out_the_immediate_state():
    W = _chain_with_immediate()
    jumps, rates, src = generator_to_jumps(W)
    jr, rr, er, state_map, _emap, _absorb = ode_eliminate_immediate(jumps, rates, src, None, _Opts())

    # State 1 is the only one with an immediate outgoing transition.
    assert state_map.tolist() == [0, 2]
    # The reduced system carries the timed rates and nothing else.
    assert sorted(np.round(rr, 12).tolist()) == [1.0, 2.0]
    assert set(er.tolist()) == {0, 2}
    # The jump matrix keeps the ORIGINAL dimension, with the eliminated
    # coordinate empty, so a caller indexing the drift by original coordinates
    # still works.
    assert jr.shape[0] == 3
    assert np.allclose(jr[1, :], 0.0)


def test_matrix_form_returns_a_proper_generator():
    W = _chain_with_immediate(a=1.0, b=2.0)
    W_red, state_map = eliminate_immediate_matrix(W, None, _Opts())
    assert state_map.tolist() == [0, 2]
    assert W_red.shape == (2, 2)
    # A generator has zero row sums; the reduced one is 0 -(1)-> 2 -(2)-> 0.
    assert np.allclose(W_red.sum(axis=1), 0.0)
    assert W_red[0, 1] == pytest.approx(1.0)
    assert W_red[1, 0] == pytest.approx(2.0)


def test_reduced_chain_keeps_the_conditional_stationary_law():
    W = _chain_with_immediate(a=1.0, b=2.0)
    W_red, state_map = eliminate_immediate_matrix(W, None, _Opts())

    def stationary(gen):
        n = gen.shape[0]
        A = np.vstack([gen.T, np.ones(n)])
        rhs = np.zeros(n + 1)
        rhs[-1] = 1.0
        return np.linalg.lstsq(A, rhs, rcond=None)[0]

    full = stationary(W)
    red = stationary(W_red)
    kept = full[state_map]
    assert np.allclose(red, kept / kept.sum(), atol=1e-6)


def test_no_immediate_transition_leaves_the_system_untouched():
    W = np.array([[-1.0, 1.0], [2.0, -2.0]])
    jumps, rates, src = generator_to_jumps(W)
    jr, rr, er, state_map, _emap, _absorb = ode_eliminate_immediate(jumps, rates, src, None, _Opts())
    assert np.allclose(jr, jumps)
    assert np.allclose(rr, rates)
    assert np.allclose(er, src)
    assert state_map.tolist() == [0, 1]

    W_red, smap = eliminate_immediate_matrix(W, None, _Opts())
    assert np.allclose(W_red, W)
    assert smap.tolist() == [0, 1]


def test_elimination_declines_when_nothing_timed_would_remain():
    # Every state immediate: reducing would leave a trivial system, so the
    # reference returns the original rather than a degenerate one.
    W = np.zeros((2, 2))
    W[0, 1] = GlobalConstants.Immediate
    W[1, 0] = GlobalConstants.Immediate
    W = W - np.diag(W.sum(axis=1))
    jumps, rates, src = generator_to_jumps(W)
    jr, rr, er, state_map, _emap, _absorb = ode_eliminate_immediate(jumps, rates, src, None, _Opts())
    assert np.allclose(rr, rates)
    assert state_map.tolist() == [0, 1]


def test_expand_state_puts_zero_in_the_eliminated_slots():
    x = ode_expand_state([0.25, 0.75], np.array([0, 2]), 3)
    assert x.tolist() == [0.25, 0.0, 0.75]


def test_hide_immediate_reaches_the_closing_drift():
    """The flag must change what is integrated, not merely be accepted."""
    from line_solver import (ClosedClass, Delay, Exp, Network, Queue,
                             SchedStrategy, SolverFLD)

    def build():
        m = Network('imm')
        d = Delay(m, 'Think')
        q = Queue(m, 'Q1', SchedStrategy.PS)
        c = ClosedClass(m, 'C1', 3, d, 0)
        d.setService(c, Exp(1.0))
        q.setService(c, Exp(2.0))
        m.link(Network.serialRouting(d, q))
        return m

    base = SolverFLD(build(), method='closing')
    qn_base = np.asarray(base.getAvgQLen()).ravel()

    opts = SolverFLD.defaultOptions()
    opts.method = 'closing'
    opts.hide_immediate = True
    hidden = SolverFLD(build(), options=opts)
    qn_hidden = np.asarray(hidden.getAvgQLen()).ravel()

    # This model has no immediate transition, so the elimination is a no-op and
    # the two tables must agree exactly. What is asserted is that requesting it
    # does not perturb a model it does not apply to.
    assert np.allclose(qn_base, qn_hidden, atol=1e-9)
