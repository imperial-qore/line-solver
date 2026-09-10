"""
LDES warm start from an auxiliary solver

An auxiliary solver is passed to SolverLDES as an argument, its steady-state
distribution is computed, and that distribution decides the initial state of
the simulation. Since the simulation then starts (approximately) in steady
state, the initialization bias vanishes and no warmup samples are discarded,
so a target accuracy is reached with fewer simulated events (and less time)
than a cold-started run.

Auxiliary-solver dispatch inside SolverLDES.initFromSolver:
 - SolverCTMC: the exact stationary distribution over the aggregate state
   space is computed and the initial state is its mode;
 - any other solver (e.g. SolverMVA): the steady-state mean queue lengths
   are rounded to a placement conserving the closed populations.

The benchmark model is a closed NEAR-BALANCED tandem, whose job split mixes
slowly: the bias of the default cold start (all jobs at the reference
station) persists beyond what the MSER-5 transient filter can remove.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import time
import numpy as np
from line_solver import (Network, Delay, Queue, ClosedClass, Exp,
                         SchedStrategy, SolverLDES, SolverMVA, SolverCTMC)


def build_model(n):
    """Closed near-balanced tandem: Think(Exp,1) -> Queue1(Exp,1.0) -> Queue2(Exp,0.98)."""
    model = Network('ldesWarmStart')
    think = Delay(model, 'Think')
    queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    jobs = ClosedClass(model, 'Jobs', n, think)
    think.setService(jobs, Exp(1.0))
    queue1.setService(jobs, Exp(1.0))
    queue2.setService(jobs, Exp(0.98))  # near-balanced bottleneck
    P = model.initRoutingMatrix()
    P.set(jobs, jobs, think, queue1, 1.0)
    P.set(jobs, jobs, queue1, queue2, 1.0)
    P.set(jobs, jobs, queue2, think, 1.0)
    model.link(P)
    return model


def run_set(n, exact, samples, seeds, init_sol):
    """Average L1 relative error and total wall-clock over the seed set."""
    err = 0.0
    rtime = 0.0
    for s in seeds:
        model = build_model(n)
        t0 = time.time()
        solver = SolverLDES(model, samples=samples, seed=s)
        if init_sol is not None:
            solver.options.init_sol = init_sol
            solver.options.tranfilter = 'fixed'
            solver.options.warmupfrac = 0.0
        qn = np.asarray(solver.getAvgQLen()).flatten()
        rtime += time.time() - t0
        err += np.abs(qn - exact).sum() / exact.sum()
    return err / len(seeds), rtime


if __name__ == '__main__':
    N = 100
    seeds = [23000, 23001, 23002]
    budgets = [10000, 50000, 200000]

    # Exact reference solution (exact MVA; the model is product-form)
    exact = np.asarray(SolverMVA(build_model(N), method='exact').getAvgQLen()).flatten()
    print('Exact mean queue lengths:', np.round(exact, 2))

    # Warm-start placement, computed ONCE by passing the auxiliary solver to
    # SolverLDES; the derived options.init_sol is then reused across runs.
    t0 = time.time()
    m = build_model(N)
    proto = SolverLDES(m, SolverMVA(m, method='exact'))
    init_sol = proto.options.init_sol
    init_time = time.time() - t0
    print('Warm placement from SolverMVA (%.3fs): %s' % (init_time, init_sol))

    # On a small instance, SolverCTMC yields the mode of the exact stationary
    # distribution instead (feasible when the state space is small):
    ms = build_model(20)
    proto_ctmc = SolverLDES(ms, SolverCTMC(ms, 'exact', force=True))
    print('CTMC distribution-mode placement (N=20):', proto_ctmc.options.init_sol)
    print()

    print('%10s  %16s  %16s' % ('samples', 'COLD err/time', 'WARM err/time'))
    cold_time, warm_time = 0.0, init_time
    cold_at, warm_at = None, None
    for b in budgets:
        ce, ct = run_set(N, exact, b, seeds, None)
        we, wt = run_set(N, exact, b, seeds, init_sol)
        cold_time += ct
        warm_time += wt
        if cold_at is None and ce < 0.10:
            cold_at = b
        if warm_at is None and we < 0.10:
            warm_at = b
        print('%10d  %7.2f%% %6.2fs  %7.2f%% %6.2fs' % (b, 100 * ce, ct, 100 * we, wt))
        if cold_at is not None and warm_at is not None:
            break
    print()
    print('Samples to reach 10%% error: COLD=%s, WARM=%s' % (cold_at, warm_at))
