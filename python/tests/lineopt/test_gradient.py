"""
Tests for the gradient-based optimizer path and analytic sensitivities.

Covers:
  - analytic performance sensitivities for open product-form networks match
    the M/M/1 closed form and finite differences;
  - the analytic objective gradient matches a finite-difference gradient;
  - the gradient optimizer path converges to the known optimum, using analytic
    sensitivities (no per-dimension solver calls in the gradient).
"""

import numpy as np
import pytest

line_solver = pytest.importorskip("line_solver")

from line_solver import (Network, Queue, Source, Sink, OpenClass, Exp,
                         SchedStrategy)
from line_solver import (OptimizationProblem, ServiceRate, MinimizeCost,
                         ResponseTimeConstraint)
from line_solver.opt.solver import LineOptSolver


def _build_mm1():
    m = Network("MM1")
    s = Source(m, "Arrivals")
    q = Queue(m, "Server", SchedStrategy.FCFS)
    k = Sink(m, "Departures")
    j = OpenClass(m, "Jobs")
    s.setArrival(j, Exp(3.0))
    q.setService(j, Exp(4.0))
    m.addLink(s, q)
    m.addLink(q, k)
    return m, q, j


def _problem():
    m, q, j = _build_mm1()
    p = OptimizationProblem(m)
    p.add_variable(ServiceRate(q, j, bounds=(3.5, 8.0)))
    p.set_objective(MinimizeCost(rate_cost={q: 20.0}))
    p.add_constraint(ResponseTimeConstraint(q, j, max_value=0.5))
    return p


def test_open_sensitivity_matches_mm1_closed_form():
    from line_solver.opt.sensitivity import compute_model_sensitivities
    m, q, j = _build_mm1()
    q.setService(j, Exp(5.0))  # mu = 5, lambda = 3
    sens = compute_model_sensitivities(m)
    assert sens is not None
    # d RespT / d mu = -1/(mu - lambda)^2 = -0.25 at mu=5
    d = sens['RespT'][('Server', 'Jobs')][('rate', 'Server', 'Jobs')]
    assert d == pytest.approx(-0.25, rel=1e-9)
    # d QLen / d mu with QLen = rho/(1-rho): -lambda/(mu-lambda)^2 = -0.75
    dq = sens['QLen'][('Server', 'Jobs')][('rate', 'Server', 'Jobs')]
    assert dq == pytest.approx(-3.0 / (2.0 ** 2), rel=1e-9)


def test_analytic_gradient_matches_finite_difference():
    p = _problem()
    sol = LineOptSolver(p, seed=42, optimizer='gradient')
    sol._start_time = 0.0
    sol._deadline = float('inf')
    sol._caches = [{}]
    sol._cache = sol._caches[0]
    sol._best_value = float('inf')
    sol._best_x = None
    for xv in (0.4, 0.5, 0.7, 0.9):
        x = np.array([xv])
        ga = sol._analyticGradient(x)
        assert ga is not None
        orig = sol._analyticGradient
        sol._analyticGradient = lambda _x: None
        gfd = sol._objectiveGradient(x)
        sol._analyticGradient = orig
        assert ga[0] == pytest.approx(gfd[0], rel=1e-4, abs=1e-4)


def test_gradient_path_converges_to_optimum():
    p = _problem()
    sol = LineOptSolver(p, seed=42, optimizer='gradient', gradient_restarts=2)
    result = sol.solve()
    rate = result.variable_values['Server_Jobs_rate']
    # optimum is the smallest rate meeting RT <= 0.5, i.e. mu* = 5
    assert rate == pytest.approx(5.0, abs=0.1)
    assert result.feasible


def test_auto_selects_gradient_for_continuous():
    p = _problem()
    sol = LineOptSolver(p, seed=42, optimizer='auto')
    assert sol._shouldUseGradient() is True


def _closed_problem():
    from line_solver import Network, Delay, Queue, ClosedClass
    m = Network("csens")
    d0 = Delay(m, "Think")
    q1 = Queue(m, "Q1", SchedStrategy.PS)
    q2 = Queue(m, "Q2", SchedStrategy.PS)
    c1 = ClosedClass(m, "C1", 5, d0)
    d0.setService(c1, Exp(1.0))
    q1.setService(c1, Exp(2.0))
    q2.setService(c1, Exp(1.8))
    m.link(Network.serialRouting(d0, q1, q2))
    p = OptimizationProblem(m)
    p.add_variable(ServiceRate(q1, c1, bounds=(1.2, 4.0)))
    p.add_variable(ServiceRate(q2, c1, bounds=(1.2, 4.0)))
    p.set_objective(MinimizeCost(rate_cost={q1: 10.0, q2: 10.0}))
    p.add_constraint(ResponseTimeConstraint(q1, c1, max_value=1.0))
    p.add_constraint(ResponseTimeConstraint(q2, c1, max_value=1.0))
    return p, m


def test_closed_network_analytic_sensitivities_available():
    from line_solver.opt.sensitivity import compute_model_sensitivities
    _, m = _closed_problem()
    assert compute_model_sensitivities(m) is not None


def test_closed_analytic_gradient_matches_fd():
    p, _ = _closed_problem()
    sol = LineOptSolver(p, seed=1, optimizer='gradient')
    sol._start_time = 0.0
    sol._deadline = float('inf')
    sol._caches = [{}]
    sol._cache = sol._caches[0]
    sol._best_value = float('inf')
    sol._best_x = None
    # probe a constraint-active point so the metric sensitivity term is exercised
    for xv in [(0.02, 0.05), (0.15, 0.15)]:
        x = np.array(xv)
        ga = sol._analyticGradient(x)
        assert ga is not None
        orig = sol._analyticGradient
        sol._analyticGradient = lambda _x: None
        gfd = sol._objectiveGradient(x)
        sol._analyticGradient = orig
        rel = np.max(np.abs(ga - gfd) / (np.abs(gfd) + 1e-3))
        assert rel < 1e-4
