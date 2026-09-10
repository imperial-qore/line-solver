"""
Tests SolverNC.getProbSysMarg, the joint law of the per-station TOTAL queue
lengths evaluated through the Calame permanent identity.

The assertions run from weakest to strongest: the law must normalize, its first
moment must reproduce the exact CTMC queue lengths, and it must keep both of
those properties on the zero-bearing demand matrices (a class that skips a
station) and on several infinite servers, each of which contributes its own
1/n_j!. The last tests fix the refusals: an approximate permanent engine on a
structurally zero demand matrix, a multiserver station, and the getter on a
solver that does not implement it.
"""

import itertools

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, Network, Queue, SchedStrategy,
                         SolverCTMC, SolverMVA, SolverNC)

TOL = 1e-9


def _dense_model():
    """Delay + two PS queues, 2 classes, N = (2,1). Dense demand matrix."""
    model = Network("marg_dense")
    d = Delay(model, "Delay")
    q1 = Queue(model, "Queue1", SchedStrategy.PS)
    q2 = Queue(model, "Queue2", SchedStrategy.PS)
    c1 = ClosedClass(model, "Class1", 2, d, 0)
    c2 = ClosedClass(model, "Class2", 1, d, 0)
    d.setService(c1, Exp.fitMean(0.7))
    d.setService(c2, Exp.fitMean(1.3))
    q1.setService(c1, Exp.fitMean(1.5))
    q1.setService(c2, Exp.fitMean(0.8))
    q2.setService(c1, Exp.fitMean(0.9))
    q2.setService(c2, Exp.fitMean(1.1))
    P = model.initRoutingMatrix()
    P.set(c1, c1, d, q1, 1.0)
    P.set(c1, c1, q1, q2, 1.0)
    P.set(c1, c1, q2, d, 1.0)
    P.set(c2, c2, d, q1, 1.0)
    P.set(c2, c2, q1, q2, 1.0)
    P.set(c2, c2, q2, d, 1.0)
    model.link(P)
    return model


def _zero_demand_model():
    """Same three stations, but class 2 never visits Queue2."""
    model = Network("marg_zero")
    d = Delay(model, "Delay")
    q1 = Queue(model, "Queue1", SchedStrategy.PS)
    q2 = Queue(model, "Queue2", SchedStrategy.PS)
    c1 = ClosedClass(model, "Class1", 2, d, 0)
    c2 = ClosedClass(model, "Class2", 1, d, 0)
    d.setService(c1, Exp.fitMean(0.7))
    d.setService(c2, Exp.fitMean(1.3))
    q1.setService(c1, Exp.fitMean(1.5))
    q1.setService(c2, Exp.fitMean(0.8))
    q2.setService(c1, Exp.fitMean(0.9))
    q2.setService(c2, Exp.fitMean(1.1))
    P = model.initRoutingMatrix()
    P.set(c1, c1, d, q1, 1.0)
    P.set(c1, c1, q1, q2, 1.0)
    P.set(c1, c1, q2, d, 1.0)
    # Class 2 bypasses Queue2 entirely.
    P.set(c2, c2, d, q1, 1.0)
    P.set(c2, c2, q1, d, 1.0)
    model.link(P)
    return model


def _two_delay_model():
    """TWO infinite servers plus one queue: each delay carries its own 1/n_j!."""
    model = Network("marg_twodelay")
    d1 = Delay(model, "Delay1")
    d2 = Delay(model, "Delay2")
    q = Queue(model, "Queue1", SchedStrategy.PS)
    c1 = ClosedClass(model, "Class1", 2, d1, 0)
    c2 = ClosedClass(model, "Class2", 1, d1, 0)
    d1.setService(c1, Exp.fitMean(0.7))
    d1.setService(c2, Exp.fitMean(1.3))
    d2.setService(c1, Exp.fitMean(1.1))
    d2.setService(c2, Exp.fitMean(0.6))
    q.setService(c1, Exp.fitMean(1.5))
    q.setService(c2, Exp.fitMean(0.8))
    P = model.initRoutingMatrix()
    P.set(c1, c1, d1, d2, 1.0)
    P.set(c1, c1, d2, q, 1.0)
    P.set(c1, c1, q, d1, 1.0)
    P.set(c2, c2, d1, d2, 1.0)
    P.set(c2, c2, d2, q, 1.0)
    P.set(c2, c2, q, d1, 1.0)
    model.link(P)
    return model


def _compositions(M, total):
    """Every composition of total into M nonnegative parts."""
    for cut in itertools.combinations(range(total + M - 1), M - 1):
        prev = -1
        out = []
        for c in cut:
            out.append(c - prev - 1)
            prev = c
        out.append(total + M - 2 - prev)
        yield out


def _law_moments(model, M, total):
    """Total mass and the induced mean queue length per station."""
    solver = SolverNC(model)
    mass = 0.0
    mean = np.zeros(M)
    for n in _compositions(M, total):
        p, _ = solver.getProbSysMarg(n)
        assert p >= -TOL, "a probability must not be negative, was %g at %s" % (p, n)
        mass += p
        mean += p * np.asarray(n, dtype=float)
    return mass, mean


def _ctmc_totals(model, M):
    """Station totals from the exact CTMC solution, classes summed out."""
    qlen = np.atleast_2d(np.asarray(SolverCTMC(model).getAvgQLen(), dtype=float))
    return np.sum(qlen, axis=1)[:M]


def test_dense_model_normalizes_and_matches_ctmc():
    """A normalizing bug and a per-state bug are distinguishable only if the
    mean is tested alongside the mass."""
    mass, mean = _law_moments(_dense_model(), 3, 3)
    assert mass == pytest.approx(1.0, abs=1e-12)
    np.testing.assert_allclose(mean, _ctmc_totals(_dense_model(), 3), atol=1e-9)


def test_zero_demand_normalizes_and_matches_ctmc():
    """A structural zero must not disturb either property under 'exact'."""
    mass, mean = _law_moments(_zero_demand_model(), 3, 3)
    assert mass == pytest.approx(1.0, abs=1e-12)
    np.testing.assert_allclose(mean, _ctmc_totals(_zero_demand_model(), 3), atol=1e-9)


def test_two_delays_normalize_and_match_ctmc():
    """Dividing once, as a single aggregated delay row would, leaves the law
    unnormalized."""
    mass, mean = _law_moments(_two_delay_model(), 3, 3)
    assert mass == pytest.approx(1.0, abs=1e-12)
    np.testing.assert_allclose(mean, _ctmc_totals(_two_delay_model(), 3), atol=1e-9)


def test_zero_population_class_and_infeasible_state():
    model = Network("marg_empty")
    d = Delay(model, "Delay")
    q = Queue(model, "Queue1", SchedStrategy.PS)
    c1 = ClosedClass(model, "Class1", 2, d, 0)
    c2 = ClosedClass(model, "Class2", 0, d, 0)
    d.setService(c1, Exp.fitMean(0.7))
    d.setService(c2, Exp.fitMean(1.3))
    q.setService(c1, Exp.fitMean(1.5))
    q.setService(c2, Exp.fitMean(0.8))
    P = model.initRoutingMatrix()
    P.set(c1, c1, d, q, 1.0)
    P.set(c1, c1, q, d, 1.0)
    P.set(c2, c2, d, q, 1.0)
    P.set(c2, c2, q, d, 1.0)
    model.link(P)

    solver = SolverNC(model)
    mass = sum(solver.getProbSysMarg(n)[0] for n in _compositions(2, 2))
    assert mass == pytest.approx(1.0, abs=1e-12)

    # A state whose total does not match the population is impossible.
    p, lp = solver.getProbSysMarg([1, 0])
    assert p == 0.0
    assert lp == -np.inf


@pytest.mark.parametrize("engine", ["spm", "bethe", "heur", "huberlaw", "adapart"])
def test_approximate_engines_refuse_a_structural_zero(engine):
    """Sinkhorn stalls without full support and the Bethe gap is a
    state-dependent lower bound that does not cancel under normalization."""
    solver = SolverNC(_zero_demand_model())
    with pytest.raises(ValueError, match="zero"):
        solver.getProbSysMarg([1, 1, 1], engine)
    # The exact engine is unaffected on the same state.
    assert solver.getProbSysMarg([1, 1, 1])[0] > 0


# The 3-station 2-class demand matrix the saddle-point engine is calibrated on.
SPM_L = np.array([[0.286, 0.437], [1.001, 0.782], [0.294, 0.633]])


def test_spm_engine_error_falls_as_the_populations_grow():
    """The saddle point is asymptotic in the COLUMN multiplicities, which here
    are the class populations, so its relative error must fall like 1/min(N)
    rather than staying flat as the sampling and mean-field engines do."""
    from line_solver.api.pfqn import pfqn_jointmarg

    previous = np.inf
    for k in (1, 2, 3, 4):
        N = np.array([k, k])
        worst, mass_exact, mass_spm = 0.0, 0.0, 0.0
        for n in _compositions(3, 2 * k):
            n = np.array(n)
            pex, _ = pfqn_jointmarg(n, SPM_L, N, None, None, 'exact')
            psp, lsp = pfqn_jointmarg(n, SPM_L, N, None, None, 'spm')
            mass_exact += pex
            mass_spm += psp
            assert abs(lsp - np.log(psp)) < 1e-12
            if pex > 1e-12:
                # the expansion overestimates the permanent, hence the probability
                assert psp >= pex * (1.0 - 1e-9)
                worst = max(worst, abs(psp - pex) / pex)
        assert mass_exact == pytest.approx(1.0, abs=1e-9)
        assert worst < 0.2 / k          # tracks 1/(8 min N_r), here R = 2
        assert worst < previous         # and it must actually improve
        previous = worst
        assert mass_spm > 1.0           # a uniform overestimate, so the mass exceeds one


def test_spm_engine_bias_cancels_under_renormalization():
    """The bias is nearly the same at every state, so a caller sweeping the
    lattice and renormalizing to sum to one is left with an order of magnitude
    less error than the raw relative error."""
    from line_solver.api.pfqn import pfqn_jointmarg

    previous = np.inf
    for k in (1, 2, 3):
        N = np.array([k, k])
        exact, spm = [], []
        for n in _compositions(3, 2 * k):
            exact.append(pfqn_jointmarg(np.array(n), SPM_L, N, None, None, 'exact')[0])
            spm.append(pfqn_jointmarg(np.array(n), SPM_L, N, None, None, 'spm')[0])
        exact = np.array(exact)
        spm = np.array(spm)
        raw = np.max(np.abs(spm - exact)[exact > 1e-12] / exact[exact > 1e-12])
        tvd = 0.5 * float(np.abs(spm / spm.sum() - exact).sum())
        assert tvd < 0.1 * raw
        assert tvd < previous
        previous = tvd


def test_spm_engine_is_reachable_from_the_solver():
    solver = SolverNC(_dense_model())
    exact, _ = solver.getProbSysMarg([1, 1, 1])
    approx, _ = solver.getProbSysMarg([1, 1, 1], 'spm')
    assert approx > 0
    assert abs(approx - exact) < 0.5 * exact


def test_multiserver_is_refused():
    model = Network("marg_multiserver")
    d = Delay(model, "Delay")
    q = Queue(model, "Queue1", SchedStrategy.FCFS)
    q.setNumberOfServers(2)
    c1 = ClosedClass(model, "Class1", 3, d, 0)
    d.setService(c1, Exp.fitMean(0.7))
    q.setService(c1, Exp.fitMean(1.5))
    P = model.initRoutingMatrix()
    P.set(c1, c1, d, q, 1.0)
    P.set(c1, c1, q, d, 1.0)
    model.link(P)

    with pytest.raises(ValueError, match="multiserver"):
        SolverNC(model).getProbSysMarg([1, 2])


def test_other_solvers_refuse_the_getter():
    with pytest.raises(NotImplementedError, match="getProbSysMarg is not supported"):
        SolverMVA(_dense_model()).getProbSysMarg([1, 1, 1])


def _started_model():
    """Delay + PS queue + FCFS queue, 2 classes: one shared server, one buffered."""
    model = Network("started")
    d = Delay(model, "Delay")
    q1 = Queue(model, "Q1", SchedStrategy.PS)
    q2 = Queue(model, "Q2", SchedStrategy.FCFS)
    c1 = ClosedClass(model, "C1", 2, d, 0)
    c2 = ClosedClass(model, "C2", 1, d, 0)
    d.setService(c1, Exp(1 / 0.7))
    d.setService(c2, Exp(1 / 1.3))
    q1.setService(c1, Exp(1 / 1.5))
    q1.setService(c2, Exp(1 / 0.8))
    q2.setService(c1, Exp(1 / 0.9))
    q2.setService(c2, Exp(1 / 1.1))
    P = model.initRoutingMatrix()
    P.set(c1, c1, Network.serialRouting(d, q1, q2))
    P.set(c2, c2, Network.serialRouting(d, q1, q2))
    model.link(P)
    return model


def _rows(space):
    return np.atleast_2d(np.asarray(space)).astype(int).tolist()


def test_from_marginal_and_started_pinned_against_matlab():
    """A SHARED SERVER holds every job present, so the state carries n and not
    s; writing s there lost the queued jobs until 2026-08-09. An ORDERED buffer
    carries the waiting class tags with the started counts in phase one."""
    from line_solver.api.state.marginal import fromMarginalAndStarted

    sn = _started_model().get_struct()
    assert _rows(fromMarginalAndStarted(sn, 1, [2, 0], [1, 0])) == [[2, 0]]
    assert _rows(fromMarginalAndStarted(sn, 1, [1, 1], [0, 1])) == [[1, 1]]
    assert _rows(fromMarginalAndStarted(sn, 2, [2, 0], [1, 0])) == [[1, 1, 0]]
    assert _rows(fromMarginalAndStarted(sn, 2, [1, 1], [0, 1])) == [[1, 0, 1]]
    # The empty station keeps the decoder's width: one buffer column here.
    assert _rows(fromMarginalAndStarted(sn, 2, [0, 0], [0, 0])) == [[0, 0, 0]]


def test_from_marg_and_started_pinned_against_matlab():
    """The union over every (n,s) pair consistent with the two totals. Class 2
    has population 1, so the split [0,2] is never enumerated."""
    from line_solver.lang.state import State

    model = _started_model()
    assert _rows(State.fromMargAndStarted(model, 1, 2, 1)) == [[2, 0], [1, 1]]
    assert _rows(State.fromMargAndStarted(model, 2, 2, 1)) == [[2, 1, 0], [1, 1, 0], [1, 0, 1]]
    assert _rows(State.fromMargAndStarted(model, 2, 0, 0)) == [[0, 0, 0]]


def test_from_marg_excludes_a_disabled_class_split():
    """Class 2 has population 1, so the station cannot hold two of them and the
    split must be dropped BEFORE fromMarginal is asked for it: an empty local
    space is absorbed by the cartesian product rather than annihilating it."""
    from line_solver.lang.state import State

    model = _started_model()
    assert _rows(State.fromMarg(model, 1, 2)) == [[2, 0], [1, 1]]
    assert _rows(State.fromMarg(model, 1, 0)) == [[0, 0]]


def test_citations_report_the_permanent():
    solver = SolverNC(_dense_model())
    solver.getProbSysMarg([1, 1, 1])
    refs = " ".join(e["ref"] for e in solver.citations())
    assert "Ryser" in refs
