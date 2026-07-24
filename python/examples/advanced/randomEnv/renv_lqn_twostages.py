"""LayeredNetwork operating in a two-stage random environment (Python native).

Mirrors matlab/examples/advanced/randomEnv/renv_lqn_twostages.m. SolverENV runs
over a LayeredNetwork (LQN) base model through the uniform model/solver
interface (no branching in ENV; meanfield analyzer only). The environment
alternates between an UP stage (fast database) and a DOWN stage (slow database);
the environment-averaged total throughput must lie between the two single-stage
LQN solutions.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
from line_solver.layered import LayeredNetwork, Processor, Task, Entry, Activity
from line_solver.constants import SchedStrategy
from line_solver.distributions import Exp
from line_solver import Environment, SolverENV, SolverLN, SolverFluid, SolverCTMC


def build_lqn(name, db_mean):
    m = LayeredNetwork(name)
    P1 = Processor(m, 'ClientProcessor', 1, SchedStrategy.PS)
    P2 = Processor(m, 'DBProcessor', 1, SchedStrategy.PS)
    T1 = Task(m, 'ClientTask', 5, SchedStrategy.REF).on(P1)
    T1.set_think_time(Exp.fit_mean(5.0))
    T2 = Task(m, 'DBTask', float('inf'), SchedStrategy.INF).on(P2)
    E1 = Entry(m, 'ClientEntry').on(T1)
    E2 = Entry(m, 'DBEntry').on(T2)
    A1 = Activity(m, 'ClientActivity', Exp.fit_mean(1.0)).on(T1)
    A1.bound_to(E1).synch_call(E2, 2.5)
    A2 = Activity(m, 'DBActivity', Exp.fit_mean(db_mean)).on(T2)
    A2.bound_to(E2).replies_to(E2)
    return m


def _tail(cell):
    # Disabled / off-block cells are None -> NaN (excluded from the totals).
    if cell is None:
        return np.nan
    if isinstance(cell, dict):
        metric = cell.get('metric')
    else:
        metric = getattr(cell, 'metric', None)
    if metric is not None and len(metric) > 0:
        return float(metric[-1])
    return np.nan


def stage_aggregate(model, T, layer_factory):
    """Steady transient aggregate (block-diagonal) Q and T for one LQN stage,
    using the same SolverLN layout SolverENV consumes."""
    s = SolverLN(model, layer_factory, timespan=[0, T], verbose=False)
    Qt, Ut, Tt = s.getTranAvg()
    M = len(Qt)
    K = len(Qt[0]) if M > 0 else 0
    Q = np.full((M, K), np.nan)
    Tm = np.full((M, K), np.nan)
    for i in range(M):
        for k in range(K):
            Q[i, k] = _tail(Qt[i][k])
            Tm[i, k] = _tail(Tt[i][k])
    return Q, Tm


def _env_tput_sum(env, layer_factory, T):
    ln_factory = lambda m: SolverLN(m, layer_factory, timespan=[0, T],
                                    iter_max=8, verbose=False)
    s = SolverENV(env, ln_factory, {'iter_max': 3, 'iter_tol': 0.05, 'verbose': False})
    _, _, TN = s.avg()
    TN = np.atleast_2d(TN)
    return float(np.nansum(TN[np.isfinite(TN)]))


def env2_stage_tput(a, b, T, layer_factory):
    """Env-averaged throughput for the two-stage UP/DOWN environment at switch
    rates a (UP->DOWN) and b (DOWN->UP); probEnv(UP)=b/(a+b)."""
    up = build_lqn('UP', 0.8)
    down = build_lqn('DOWN', 3.0)
    env = Environment('R', 2)
    env.addStage(0, 'UP', 'operational', up)
    env.addStage(1, 'DOWN', 'degraded', down)
    env.addTransition(0, 1, Exp(a))
    env.addTransition(1, 0, Exp(b))
    return _env_tput_sum(env, layer_factory, T)


def build_lqn_small(name, db_mean):
    """Compact LQN (unit multiplicity, a single DB call) whose layer submodels
    have tiny CTMC state spaces, so the three-stage case can use exact and
    numerically stable CTMC layer solvers in native Python."""
    model = LayeredNetwork(name)
    P1 = Processor(model, 'CP', 1, SchedStrategy.PS)
    P2 = Processor(model, 'DP', 1, SchedStrategy.PS)
    T1 = Task(model, 'CT', 2, SchedStrategy.REF).on(P1)
    T1.set_think_time(Exp.fit_mean(4.0))
    T2 = Task(model, 'DT', float('inf'), SchedStrategy.INF).on(P2)
    E1 = Entry(model, 'CE').on(T1)
    E2 = Entry(model, 'DE').on(T2)
    A1 = Activity(model, 'CA', Exp.fit_mean(1.0)).on(T1)
    A1.bound_to(E1).synch_call(E2, 1.0)
    A2 = Activity(model, 'DA', Exp.fit_mean(db_mean)).on(T2)
    A2.bound_to(E2).replies_to(E2)
    return model


def stage_ctmc_tput(model, T):
    """Single-stage aggregate steady throughput via SolverLN with CTMC layers."""
    Qt, Ut, Tt = SolverLN(model, lambda mm: SolverCTMC(mm, verbose=False),
                          timespan=[0, T], iter_max=6, verbose=False).getTranAvg()
    tot = 0.0
    for i in range(len(Tt)):
        for k in range(len(Tt[i])):
            v = _tail(Tt[i][k])
            if np.isfinite(v):
                tot += v
    return tot


def env3_stage_ctmc(T):
    """Three-stage UP/MID/DOWN environment over compact LQN stages with CTMC
    layers. Returns (env throughput, fast-stage throughput, slow-stage
    throughput)."""
    up = build_lqn_small('UP', 0.5)
    mid = build_lqn_small('MID', 1.0)
    down = build_lqn_small('DOWN', 2.0)
    env = Environment('R3', 3)
    env.addStage(0, 'UP', 'operational', up)
    env.addStage(1, 'MID', 'degraded', mid)
    env.addStage(2, 'DOWN', 'failed', down)
    env.addTransition(0, 1, Exp(0.3))
    env.addTransition(1, 2, Exp(0.3))
    env.addTransition(2, 1, Exp(0.6))
    env.addTransition(1, 0, Exp(0.6))
    ctmc = lambda mm: SolverCTMC(mm, verbose=False)
    ln_factory = lambda m: SolverLN(m, ctmc, timespan=[0, T], iter_max=6, verbose=False)
    s = SolverENV(env, ln_factory, {'iter_max': 2, 'iter_tol': 0.1, 'verbose': False})
    _, _, TN = s.avg()
    TN = np.atleast_2d(TN)
    x_env = float(np.nansum(TN[np.isfinite(TN)]))
    return x_env, stage_ctmc_tput(up, T), stage_ctmc_tput(down, T)


def main():
    # Fluid layer solvers (matching the MATLAB example); a short transient
    # window keeps the native fluid ODE well-conditioned, and the SolverLN
    # inner iteration is capped so the meanfield loop stays fast in Python.
    T = 5
    layer_factory = lambda mm: SolverFluid(mm, verbose=False)
    ln_factory = lambda m: SolverLN(m, layer_factory, timespan=[0, T],
                                    iter_max=8, verbose=False)

    up = build_lqn('LQN_UP', 0.8)
    down = build_lqn('LQN_DOWN', 3.0)

    env = Environment('DBReliability', 2)
    env.addStage(0, 'UP', 'operational', up)
    env.addStage(1, 'DOWN', 'degraded', down)
    env.addTransition(0, 1, Exp(0.2))   # mean UP time = 5
    env.addTransition(1, 0, Exp(1.0))   # mean DOWN time = 1

    opts = {'iter_max': 3, 'iter_tol': 0.05, 'verbose': False}
    env_solver = SolverENV(env, ln_factory, opts)
    QN, UN, TN = env_solver.avg()
    QN = np.atleast_2d(QN)
    TN = np.atleast_2d(TN)
    print("ENV over LQN ran: aggregate size = %d x %d" % (QN.shape[0], QN.shape[1]))
    print(env_solver.getAvgTable())

    upQ, upT = stage_aggregate(up, T, layer_factory)
    downQ, downT = stage_aggregate(down, T, layer_factory)

    Xup = float(np.nansum(upT))
    Xdown = float(np.nansum(downT))
    Xenv = float(np.nansum(TN[np.isfinite(TN)]))
    QsumEnv = float(np.nansum(QN))
    QsumUp = float(np.nansum(upQ))
    QsumDown = float(np.nansum(downQ))
    print("Total throughput  UP=%.4f  DOWN=%.4f  ENV=%.4f" % (Xup, Xdown, Xenv))
    print("Aggregate Q (sum) UP=%.4f  DOWN=%.4f  ENV=%.4f" % (QsumUp, QsumDown, QsumEnv))

    # (1) The solver runs and returns finite metrics.
    assert np.all(np.isfinite(QN)), "ENV aggregate Q contains non-finite values"
    # (2) Population conservation (loose): the env-averaged aggregate is a
    #     CDF-weighted transient average, so it tracks the steady single-stage
    #     population to within a small transient tolerance.
    assert abs(QsumEnv - QsumUp) <= 0.05 * QsumUp and \
        abs(QsumEnv - QsumDown) <= 0.05 * QsumDown, \
        "ENV aggregate does not conserve the closed population of the stages"
    # (3) A physically monotone scalar (total throughput) lies between the two
    #     single-stage solutions: a slower database slows the whole system.
    lo, hi = min(Xup, Xdown), max(Xup, Xdown)
    tol = 1e-2 * max(1.0, hi)
    assert lo - tol <= Xenv <= hi + tol, \
        "ENV-averaged throughput is not bracketed by the single-stage solutions"

    # (4) Quantitative coupling: env-averaged throughput increases with the
    #     stationary probability of the fast UP stage, P(UP)=b/(a+b).
    x_low = env2_stage_tput(1.0, 0.2, T, layer_factory)   # P(UP)=0.167
    x_high = env2_stage_tput(0.2, 1.0, T, layer_factory)  # P(UP)=0.833
    print("Monotonicity   P(UP)=0.167 -> %.4f   P(UP)=0.833 -> %.4f" % (x_low, x_high))
    assert x_high > x_low + 1e-3, "ENV throughput is not monotone in P(UP)"
    assert lo - tol <= x_low and x_high <= hi + tol, \
        "ENV throughputs escape the single-stage bracket"

    # (5) Three-stage environment (UP/MID/DOWN), exercising the E>2 coupling.
    #     Uses CTMC layers on a compact LQN: exact and numerically stable, unlike
    #     the native fluid ODE which is fragile for the deeper coupling. The
    #     environment-averaged throughput stays bracketed by the fastest and
    #     slowest single-stage solutions.
    x3, t3fast, t3slow = env3_stage_ctmc(4.0)
    lo3, hi3 = min(t3fast, t3slow), max(t3fast, t3slow)
    tol3 = 1e-2 * max(1.0, hi3)
    print("Three-stage (CTMC) ENV throughput = %.4f (bracket [%.4f, %.4f])" % (x3, lo3, hi3))
    assert np.isfinite(x3), "Three-stage ENV throughput is not finite"
    assert lo3 - tol3 <= x3 <= hi3 + tol3, \
        "Three-stage ENV-averaged throughput is not bracketed"

    print("PASS: ENV-over-LQN meanfield ran; throughput bracketed, monotone in P(UP), 3-stage bracketed.")


if __name__ == '__main__':
    main()
