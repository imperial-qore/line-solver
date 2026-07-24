"""Exact symbolic sensitivity (SolverCTMC.getSensitivity) and fluid Jacobian
(SolverFLD.getJacobian / getSymbolicDrift).

The drift is checked against the field the solver actually integrates,
recomputed here in double precision, so a divergence in the exported
expression shows up as a number rather than as a differently spelled formula.
The Jacobian is checked against a central difference of that same field, and
the CTMC sensitivity against a brute-force difference of the reward itself.

Counterpart of the JAR SolverFluidSymbolicSageTest and of the MATLAB
@SolverFLD/getJacobian and @SolverCTMC/getSensitivity. The native engine is
sympy, so nothing here needs the line-sage-rest service.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from line_solver import (Network, Queue, Delay, ClosedClass, Exp, Erlang,
                         SchedStrategy)
from line_solver.solvers.solver_ctmc import SolverCTMC
from line_solver.solvers.solver_fld import SolverFLD
from line_solver.solvers.solver_fld.ode.pnorm import fluid_ode_pnorm

sympy = pytest.importorskip('sympy')
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                        convert_xor)

TRANSFORMATIONS = standard_transformations + (convert_xor,)


def _evaluate(exprs, variables, x):
    """Numeric value of expression strings written with '^' for powers."""
    syms = [sympy.Symbol(v, real=True) for v in variables]
    local = dict(zip(variables, syms))
    subs = dict(zip(syms, x))
    return np.array([float(parse_expr(e, local_dict=local,
                                      transformations=TRANSFORMATIONS).subs(subs))
                     for e in exprs])


def fluid_model():
    """Closed model: Delay Erlang(1, 2) -> PS Queue Exp(2) with 2 servers, N = 4."""
    model = Network('FluidSymbolic')
    delay = Delay(model, 'Think')
    queue = Queue(model, 'Queue1', SchedStrategy.PS)
    queue.setNumberOfServers(2)
    jobs = ClosedClass(model, 'Class1', 4, delay)
    delay.setService(jobs, Erlang.fitMeanAndOrder(1.0, 2))
    queue.setService(jobs, Exp(2.0))
    model.link(Network.serialRouting(delay, queue))
    return model


def integrated_field(sys, x):
    """The field the p-norm method integrates, evaluated at x."""
    n = sys.nstates
    SQ = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            SQ[a, b] = 1.0 if sys.stateStation[a] == sys.stateStation[b] else 0.0
    Sa = np.array([sys.S[sys.stateStation[a]] for a in range(n)])
    pstar = np.array([sys.pstar[sys.stateStation[a]] for a in range(n)])
    return fluid_ode_pnorm(0.0, np.asarray(x, dtype=float), sys.W.T, SQ, Sa,
                           sys.Alambda, pstar)


def test_min_scaled_drift_is_refused_not_one_sidedly_differentiated():
    solver = SolverFLD(fluid_model(), method='matrix')
    with pytest.raises(ValueError) as excinfo:
        solver.getSymbolicDrift()
    assert 'min(n_i, S_i)' in str(excinfo.value)


def test_smooth_drift_matches_the_integrated_field():
    solver = SolverFLD(fluid_model(), method='matrix', pstar=8.0)
    rhs, variables, sys = solver.getSymbolicDrift()
    n = sys.nstates
    assert len(rhs) == n
    x = np.array([0.5 + 0.37 * (s + 1) for s in range(n)])
    assert np.allclose(_evaluate(rhs, variables, x), integrated_field(sys, x),
                       atol=1e-9, rtol=0)


def test_jacobian_matches_a_central_difference_of_the_field():
    solver = SolverFLD(fluid_model(), method='matrix', pstar=8.0)
    J, rhs, variables, equilibria = solver.getJacobian()
    sys = solver.getSymbolicDrift()[2]
    n = sys.nstates
    assert len(J) == n
    assert equilibria is None
    x = np.array([0.5 + 0.37 * (s + 1) for s in range(n)])
    h = 1e-6
    numeric = np.zeros((n, n))
    for j in range(n):
        xp = x.copy()
        xm = x.copy()
        xp[j] += h
        xm[j] -= h
        numeric[:, j] = (integrated_field(sys, xp) - integrated_field(sys, xm)) / (2 * h)
    exact = np.array([_evaluate(row, variables, x) for row in J])
    assert np.allclose(exact, numeric, atol=1e-6, rtol=0)


def ctmc_model(mu=2.0):
    """Closed model: Delay Exp(1) -> PS Queue Exp(mu), N = 2."""
    model = Network('CtmcSensitivity')
    delay = Delay(model, 'Think')
    queue = Queue(model, 'Queue1', SchedStrategy.PS)
    jobs = ClosedClass(model, 'Class1', 2, delay)
    delay.setService(jobs, Exp(1.0))
    queue.setService(jobs, Exp(mu))
    model.link(Network.serialRouting(delay, queue))
    return model


def _set_service_rate(model, value):
    model.getStationByName('Queue1').setService(
        model.getClassByName('Class1'), Exp(value))


PARAM = {'name': 'mu', 'value': 2.0, 'set': _set_service_rate}


def mean_queue_length(mu):
    solver = SolverCTMC(ctmc_model(mu))
    pi = np.asarray(solver.getSteadyState(), dtype=float).flatten()
    return float(pi @ solver.getStateSpaceAggr()[:, 1])


@pytest.mark.parametrize('method', ['fd', 'symbolic'])
def test_sensitivity_matches_a_brute_force_difference(method):
    solver = SolverCTMC(ctmc_model())
    reward = solver.getStateSpaceAggr()[:, 1]
    S, SS, dpi, pi = solver.getSensitivity(PARAM, reward, method)

    # E[Q] = (2/mu + 4/mu^2) / (1 + 2/mu + 2/mu^2), so dE[Q]/dmu is -0.28 at
    # mu = 2 in closed form; the brute-force difference confirms the model
    # builds the chain that formula describes.
    h = 1e-4
    numeric = (mean_queue_length(2.0 + h) - mean_queue_length(2.0 - h)) / (2 * h)
    assert S == pytest.approx(-0.28, abs=1e-8)
    assert S == pytest.approx(numeric, abs=1e-6)
    assert SS == pytest.approx((2.0 / 0.8) * S, rel=1e-12)
    assert np.isclose(pi.sum(), 1.0)
    assert np.isclose(dpi.sum(), 0.0, atol=1e-9)


def test_symbolic_sensitivity_is_more_accurate_than_finite_differences():
    """The exact d(pi)/d(x_e) removes the O(step^2) term of the 'fd' method.

    The remaining error of 'symbolic' is in the rate map alone, which is
    affine here, so it is at round-off; 'fd' differences through the solve and
    is not.
    """
    solver = SolverCTMC(ctmc_model())
    reward = solver.getStateSpaceAggr()[:, 1]
    S_fd = solver.getSensitivity(PARAM, reward, 'fd')[0]
    S_sym = solver.getSensitivity(PARAM, reward, 'symbolic')[0]
    assert abs(S_sym + 0.28) < abs(S_fd + 0.28)
    assert abs(S_sym + 0.28) < 1e-12


def _scale_all_rates(model, value):
    model.getStationByName('Think').setService(
        model.getClassByName('Class1'), Exp(1.0 * value))
    model.getStationByName('Queue1').setService(
        model.getClassByName('Class1'), Exp(2.0 * value))


@pytest.mark.parametrize('method', ['fd', 'symbolic'])
def test_uniform_time_rescaling_has_zero_sensitivity(method):
    """Scaling every rate by theta rescales time only, so pi is invariant.

    Every event depends on theta here, so the chain rule has to sum several
    non-zero terms to exactly zero: a term dropped or mis-scaled cannot cancel.
    """
    solver = SolverCTMC(ctmc_model())
    reward = solver.getStateSpaceAggr()[:, 1]
    param = {'name': 'scale', 'value': 1.0, 'set': _scale_all_rates}
    S, SS, dpi, pi = solver.getSensitivity(param, reward, method)
    assert S == pytest.approx(0.0, abs=1e-9)
    assert np.allclose(dpi, 0.0, atol=1e-9)


def test_unknown_method_is_rejected():
    solver = SolverCTMC(ctmc_model())
    with pytest.raises(ValueError):
        solver.getSensitivity(PARAM, None, 'central')
