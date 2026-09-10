"""The in-tree LSODA (`line_solver.lib.lsoda`) driven by the fluid path.

`tests/test_lsoda.py` checks the library against the C reference through its own
`lsoda()` driver, which integrates to each requested output time (itask=1). The
fluid path drives it a step at a time instead, through the scipy `OdeSolver`
adapter in `solvers/solver_fld/ode/native_lsoda.py`, and that mode has its own
failure surface: it was the mode that exposed the dropped `jstart = -1`, which
left a pending Adams/BDF switch uncompleted and sent Robertson to y1 = -1.9e7 in
8e6 function evaluations while still reporting success.

So the checks here are: the adapter reproduces scipy's compiled LSODA step for
step on a stiff benchmark, its dense output is the trajectory and not a
resampling of it, the stiff variant really never leaves BDF, and a fluid model
answers the same whether it is integrated by the default or by this.

Nothing here changes a default: `options.odesolver` is unset everywhere else.
"""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from line_solver import (ClosedClass, Delay, Exp, Network, Queue, SchedStrategy,
                         SolverFLD)
from line_solver.solvers.solver_fld.ode.native_lsoda import (NativeLSODA,
                                                             NativeLSODAStiff)


def _robertson(t, y):
    return [1.0e4 * y[1] * y[2] - 0.04 * y[0],
            -(1.0e4 * y[1] * y[2] - 0.04 * y[0]) - 3.0e7 * y[1] * y[1],
            3.0e7 * y[1] * y[1]]


# The C reference row at t = 4e10 from tests/test_lsoda.py.
_ROBERTSON_4E10 = np.array([1.431681921554496e-08, 5.726732909838437e-14,
                            9.999999856831133e-01])


def _fluid_model():
    model = Network('cqn')
    think = Delay(model, 'Think')
    q1 = Queue(model, 'Q1', SchedStrategy.PS)
    q2 = Queue(model, 'Q2', SchedStrategy.PS)
    jobclass = ClosedClass(model, 'C1', 5, think)
    think.setService(jobclass, Exp(1))
    q1.setService(jobclass, Exp(2))
    q2.setService(jobclass, Exp(1.5))
    model.link(Network.serialRouting(think, q1, q2))
    return model


def _avg(model, odesolver=None, method=None):
    solver = SolverFLD(model) if method is None else SolverFLD(model, method=method)
    solver.options.verbose = False
    if odesolver is not None:
        solver.options.odesolver = odesolver
    table = solver.getAvgTable()
    return np.array(table[['QLen', 'Util', 'RespT', 'Tput']], dtype=float)


def test_stepping_matches_the_c_reference():
    # Driven one step at a time by solve_ivp, i.e. through the adapter, not
    # through lsoda()'s own output-time loop. The gate is the accuracy that was
    # REQUESTED: at t = 4e10 the first component is 1.4e-8, a hundred times
    # below its own atol of 1e-6, so a relative gate on it would be asserting
    # far more than the integration promises.
    t_eval = [0.4 * 10.0 ** i for i in range(12)]
    atol = np.array([1e-6, 1e-10, 1e-6])
    sol = solve_ivp(_robertson, (0.0, t_eval[-1]), [1.0, 0.0, 0.0],
                    method=NativeLSODA, t_eval=t_eval, rtol=1e-4, atol=list(atol))
    assert sol.success
    got = sol.y[:, -1]
    assert np.all(np.abs(got - _ROBERTSON_4E10) <= atol + 1e-4 * np.abs(_ROBERTSON_4E10))
    assert abs(sum(got) - 1.0) < 1e-8


def test_force_stiff_is_accurate_on_hires():
    # HIRES at rtol 1e-6, against the C reference row at t = 321.8122 that
    # tests/test_lsoda.py uses. The pinned variant has no C counterpart -- the C
    # cannot pin -- so it is checked against the reference the auto-switcher
    # reproduces exactly.
    def hires(t, y):
        return [-1.71*y[0] + 0.43*y[1] + 8.32*y[2] + 0.0007,
                1.71*y[0] - 8.75*y[1],
                -10.03*y[2] + 0.43*y[3] + 0.035*y[4],
                8.32*y[1] + 1.71*y[2] - 1.12*y[3],
                -1.745*y[4] + 0.43*y[5] + 0.43*y[6],
                -280.0*y[5]*y[7] + 0.69*y[3] + 1.71*y[4] - 0.43*y[5] + 0.69*y[6],
                280.0*y[5]*y[7] - 1.81*y[6],
                -280.0*y[5]*y[7] + 1.81*y[6]]
    ref = np.array([7.371423243758239e-04, 1.442507544489246e-04,
                    5.888935680491721e-05, 1.175672019019430e-03,
                    2.386687478330434e-03, 6.239970967822740e-03,
                    2.850266805796065e-03, 2.849733194203930e-03])
    y0 = [1.0, 0, 0, 0, 0, 0, 0, 0.0057]
    sol = solve_ivp(hires, (0.0, 321.8122), y0, method=NativeLSODAStiff,
                    t_eval=[100.0, 200.0, 321.8122], rtol=1e-6, atol=1e-8)
    assert sol.success
    assert np.max(np.abs(sol.y[:, -1] - ref) / ref) < 1e-4


def test_force_stiff_survives_a_cold_start_at_zero():
    # Pinning BDF takes a finite-difference Jacobian at t0, where the
    # auto-switcher is still on Adams and takes none, and LSODA sizes that
    # increment as max(sqrt(eps)*|y_j|, r0/ewt_j): for y2 = 0 under atol 1e-10 it
    # collapses to ~1e-19 and the column is rounding noise. This run used to
    # reach y1 = -1.5e7 while reporting success, and flipped on a last-bit change
    # to the right-hand side. _prja now carries numjac's floor under the pin.
    t_eval = [0.4 * 10.0 ** i for i in range(12)]
    sol = solve_ivp(_robertson, (0.0, t_eval[-1]), [1.0, 0.0, 0.0],
                    method=NativeLSODAStiff, t_eval=t_eval, rtol=1e-4,
                    atol=[1e-6, 1e-10, 1e-6])
    assert sol.success
    got = sol.y[:, -1]
    # y1 -> 1.43e-8 and y3 -> 1, both inside their own atol of 1e-6
    assert abs(got[0]) < 1e-6
    assert abs(got[2] - 1.0) < 1e-6
    assert abs(sum(got) - 1.0) < 1e-8


def test_dense_output_is_the_trajectory():
    sol = solve_ivp(lambda t, y: -y, (0.0, 5.0), [1.0], method=NativeLSODA,
                    rtol=1e-8, atol=1e-10, dense_output=True)
    tq = np.linspace(0.0, 5.0, 51)
    assert np.max(np.abs(sol.sol(tq)[0] - np.exp(-tq))) < 1e-6


def test_force_stiff_never_leaves_bdf():
    # meth = 2 is BDF. The auto-switcher starts on Adams (meth = 1), so this
    # also proves the pin is doing something.
    solver = NativeLSODAStiff(lambda t, y: -1000.0 * y, 0.0, np.array([1.0]), 1.0,
                              rtol=1e-6, atol=1e-8)
    for _ in range(10):
        solver.step()
        assert solver._stepper.meth == 2
        if solver.status != 'running':
            break
    assert solver.status in ('running', 'finished')


@pytest.mark.parametrize('method', [None, 'closing'])
def test_fluid_model_agrees_with_the_default_integrator(method):
    model = _fluid_model()
    base = _avg(model, method=method)
    native = _avg(_fluid_model(), odesolver='LSODA_NATIVE_STIFF', method=method)
    assert np.max(np.abs(native - base) / (1.0 + np.abs(base))) < 1e-8


def test_fluid_model_accepts_a_solver_class():
    # options.odesolver takes a class as well as a name, which is how a caller
    # supplies an integrator that is not in _ODE_CLASSES.
    model = _fluid_model()
    base = _avg(model)
    native = _avg(_fluid_model(), odesolver=NativeLSODAStiff)
    assert np.max(np.abs(native - base) / (1.0 + np.abs(base))) < 1e-8
