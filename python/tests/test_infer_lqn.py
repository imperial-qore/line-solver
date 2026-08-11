# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.
"""Tests for native LQN parameter identification (infer_lqn, EKF).

Method: Zheng, Yang, Woodside, Litoiu, Iszlai, "Tracking Time-Varying
Parameters in Software Systems with Extended Kalman Filters", CASCON 2005.
Two hidden parameters are estimated: the reference-task think time (the paper's
Z) and the P2 host demand of activity AS3 (the paper's S_d), from the
observation vector [R(E1), U(P1), U(P2)].
"""
import numpy as np

from line_solver.layered import LayeredNetwork, Processor, Task, Entry, Activity
from line_solver.constants import SchedStrategy
from line_solver.distributions import Exp
from line_solver import SolverLN
from line_solver.inference import (
    infer_lqn, infer_lqn_jacobian, infer_lqn_setparams, infer_lqn_getobs,
)


def build_lqn():
    # deliberately light populations: native SolverLN is ~5x slower than MATLAB,
    # and this size recovers with margin while keeping the test CI-friendly.
    model = LayeredNetwork('paramident_LQN')
    P1 = Processor(model, 'P1', 1, SchedStrategy.PS)
    P2 = Processor(model, 'P2', 1, SchedStrategy.PS)
    T1 = Task(model, 'T1', 5, SchedStrategy.REF).on(P1)
    T1.set_think_time(Exp.fit_mean(0.5))
    T2 = Task(model, 'T2', 5, SchedStrategy.FCFS).on(P1)
    T2.set_think_time(Exp.fit_mean(1.0 / 3.0))
    T3 = Task(model, 'T3', 3, SchedStrategy.FCFS).on(P2)
    T3.set_think_time(Exp.fit_mean(0.25))
    E1 = Entry(model, 'E1').on(T1)
    E2 = Entry(model, 'E2').on(T2)
    E3 = Entry(model, 'E3').on(T3)
    A1 = Activity(model, 'AS1', Exp.fit_mean(0.1)).on(T1)
    A1.bound_to(E1).synch_call(E2, 1)
    A2 = Activity(model, 'AS2', Exp.fit_mean(0.05)).on(T2)
    A2.bound_to(E2).synch_call(E3, 5).replies_to(E2)
    A3 = Activity(model, 'AS3', Exp.fit_mean(0.02)).on(T3)
    A3.bound_to(E3).replies_to(E3)
    return model


def solve_obs(model, param_spec, obs_spec, a):
    infer_lqn_setparams(model, param_spec, a)
    solver = SolverLN(model, verbose=False)
    QN, UN, RN, TN, _AN, _WN = solver.get_ensemble_avg()
    names = solver.lqn.names
    metrics = {'QLen': QN, 'Util': UN, 'RespT': RN, 'Tput': TN}
    return infer_lqn_getobs(names, metrics, obs_spec)


def test_noisefree_recovery():
    model = build_lqn()
    param_spec = [{'type': 'think', 'name': 'T1'},
                  {'type': 'hostdem', 'name': 'AS3'}]
    obs_spec = [{'metric': 'RespT', 'name': 'E1'},
                {'metric': 'Util', 'name': 'P1'},
                {'metric': 'Util', 'name': 'P2'}]
    a_true = np.array([0.5, 1.0 / 50.0])

    z = solve_obs(model, param_spec, obs_spec, a_true)
    # Each EKF step costs 1 + len(param_spec) SolverLN solves (base + one
    # forward-difference column per parameter). On a noise-free constant Z the
    # filter is converged by step 4; the tail only refines the third digit, so
    # 8 steps keep >2x margin on every assertion below at a quarter of the cost.
    Z = np.tile(z.reshape(-1, 1), (1, 8))

    opt = {'a0': np.array([0.8, 1.0 / 35.0]), 'QFac': 1e-3, 'RFac': 0.2,
           'gammaT': 1.0, 'aTrue': a_true}
    _model, info = infer_lqn(model, param_spec, obs_spec, Z, opt)

    relerr = np.abs(info['ahat'][:, -1] - a_true) / a_true
    assert relerr[0] < 1e-2, 'think-time relerr too large: %g' % relerr[0]
    assert relerr[1] < 5e-3, 'host-demand relerr too large: %g' % relerr[1]
    assert info['Er'] < 1e-1


def test_util_monotone_in_demand():
    model = build_lqn()
    param_spec = [{'type': 'think', 'name': 'T1'},
                  {'type': 'hostdem', 'name': 'AS3'}]
    obs_spec = [{'metric': 'Util', 'name': 'P2'}]
    z_lo = solve_obs(model, param_spec, obs_spec, np.array([0.5, 1.0 / 50.0]))
    z_hi = solve_obs(model, param_spec, obs_spec, np.array([0.5, 1.0 / 25.0]))
    assert z_hi[0] > z_lo[0]


def test_jacobian_analytic():
    hfun = lambda a: np.array([a[0] ** 2, 2 * a[1], a[0] * a[1]])
    a = np.array([3.0, 5.0])
    H, h0 = infer_lqn_jacobian(hfun, a, 1e-6, 1e-9)
    Jexact = np.array([[2 * a[0], 0.0], [0.0, 2.0], [a[1], a[0]]])
    assert np.allclose(h0, hfun(a), atol=1e-12)
    assert np.allclose(H, Jexact, atol=1e-3)
