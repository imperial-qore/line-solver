#!/usr/bin/env python3
"""Identify hidden LQN parameters from measured performance (native Python).

Demonstrates infer_lqn, an Extended Kalman Filter that tracks hidden LQN
parameters (host demands, think times) from measurable performance data,
following Zheng, Yang, Woodside, Litoiu, Iszlai, "Tracking Time-Varying
Parameters in Software Systems with Extended Kalman Filters", CASCON 2005.

Two parameters are hidden: the reference-task think time (the paper's Z) and the
P2 host demand of activity AS3 (the paper's service demand S_d). The measurable
vector is [R(E1), U(P1), U(P2)]. A measurement sequence with a step change plus
noise is synthesised, then the parameter trajectory is recovered.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

import numpy as np

from line_solver.layered import LayeredNetwork, Processor, Task, Entry, Activity
from line_solver.constants import SchedStrategy
from line_solver.distributions import Exp
from line_solver import SolverLN
from line_solver.inference import infer_lqn, infer_lqn_setparams, infer_lqn_getobs


def build_lqn():
    # light populations keep the native SolverLN loop fast
    model = LayeredNetwork('paramident_LQN')
    P1 = Processor(model, 'P1', 1, SchedStrategy.PS)
    P2 = Processor(model, 'P2', 1, SchedStrategy.PS)
    T1 = Task(model, 'T1', 5, SchedStrategy.REF).on(P1)
    T1.set_think_time(Exp.fit_mean(0.5))
    T2 = Task(model, 'T2', 5, SchedStrategy.FCFS).on(P1)
    T2.set_think_time(Exp.fit_mean(1.0 / 3.0))
    T3 = Task(model, 'T3', 3, SchedStrategy.FCFS).on(P2)
    T3.set_think_time(Exp.fit_mean(0.25))
    Entry(model, 'E1').on(T1)
    Entry(model, 'E2').on(T2)
    Entry(model, 'E3').on(T3)
    E1, E2, E3 = model.entries
    Activity(model, 'AS1', Exp.fit_mean(0.1)).on(T1).bound_to(E1).synch_call(E2, 1)
    Activity(model, 'AS2', Exp.fit_mean(0.05)).on(T2).bound_to(E2).synch_call(E3, 5).replies_to(E2)
    Activity(model, 'AS3', Exp.fit_mean(0.02)).on(T3).bound_to(E3).replies_to(E3)
    return model


def main():
    np.random.seed(12345)
    model = build_lqn()
    param_spec = [{'type': 'think', 'name': 'T1'},
                  {'type': 'hostdem', 'name': 'AS3'}]
    obs_spec = [{'metric': 'RespT', 'name': 'E1'},
                {'metric': 'Util', 'name': 'P1'},
                {'metric': 'Util', 'name': 'P2'}]

    nsteps = 16
    a_true_seq = np.zeros((2, nsteps))
    a_true_seq[0, :] = 0.5
    a_true_seq[1, :] = 1.0 / 50.0
    a_true_seq[0, nsteps // 2:] = 1.0            # think time doubles mid-run
    a_true_seq[1, 4:12] = 1.0 / 25.0            # S_d pulse

    no = len(obs_spec)
    Z = np.zeros((no, nsteps))
    for k in range(nsteps):
        infer_lqn_setparams(model, param_spec, a_true_seq[:, k])
        solver = SolverLN(model, verbose=False)
        QN, UN, RN, TN, _AN, _WN = solver.get_ensemble_avg()
        z = infer_lqn_getobs(solver.lqn.names,
                             {'QLen': QN, 'Util': UN, 'RespT': RN, 'Tput': TN}, obs_spec)
        Z[:, k] = z * (1 + 0.02 * np.random.randn(no))   # ~2% measurement noise

    opt = {'a0': np.array([0.7, 1.0 / 40.0]), 'QFac': 0.1, 'RFac': 0.2,
           'gammaT': 1.0, 'aTrue': a_true_seq[:, -1]}
    _model, info = infer_lqn(model, param_spec, obs_spec, Z, opt)

    print('\nStep |  Z_true  Z_hat |  Sd_true  Sd_hat | ||e||')
    for k in range(nsteps):
        print('%4d | %6.3f  %6.3f | %7.4f  %7.4f | %.3g' % (
            k + 1, a_true_seq[0, k], info['ahat'][0, k],
            a_true_seq[1, k], info['ahat'][1, k],
            float(np.linalg.norm(info['e'][:, k]))))
    print('\nFinal: Z = %.4f (true %.4f), Sd = %.5f (true %.5f)' % (
        info['ahat'][0, -1], a_true_seq[0, -1], info['ahat'][1, -1], a_true_seq[1, -1]))
    print('RMS tracking Ea = %.4g, prediction Er = %.4g' % (info['Ea'], info['Er']))


if __name__ == '__main__':
    main()
