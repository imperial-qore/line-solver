"""
renv_container_terminal - Daily-cycle container terminal in a random environment

Advanced example for ENV's STATE-VECTOR analyzer (options['method'] =
'statevec'), which carries the full per-stage state distribution across
environment switches instead of collapsing it to marginal mean queue lengths.

Case study (after the fyp26 SOQN blending port): the Rotterdam container
terminal of Dhingra et al. Container handling demand varies over the 24 hours of
a day; each hour is one environment "stage" with its own demand intensity and
(random, exponential) duration. The environment visits the 24 hourly stages in a
fixed daily cycle 1->2->...->24->1.

Here the semi-open SOQN of the original study is rendered as the equivalent
finite-token CLOSED network required by ENV: N straddle carriers cycle between a
yard staging Delay and a multi-server quay-crane Queue. The hourly demand
modulates the yard staging rate (the closed-network analog of the time-varying
external arrival intensity): busier hours stage containers to the cranes faster,
so the cranes congest during the daily peaks.

The script compares, for the day-averaged metrics:
  (a) statevec  - full state-vector blending (this feature),
  (b) meanfield - the default mean-field (marginal mean-queue-length) coupling,
  (c) exact     - the stationary solution of the full joint
                  (hour x network-state) random-environment CTMC.
The state-vector blend tracks the exact joint solution far more closely than the
mean-field collapse.

A second section repeats the comparison with the internal terminal handling
(quay cranes -> stacking cranes) collapsed by Norton's theorem into a single
closed, load-dependent Flow-Equivalent Server (FES). This shows that a closed
FES (load-dependent station, sn.lldscaling) is fully supported by the
state-vector analyzer's CTMC backend.

A third section analyses the OPEN counterpart: the hourly demand becomes an
external Markov-modulated Poisson arrival stream (MMPP(24)) into a multi-server
queue with an unbounded buffer, solved exactly by MAM as a QBD. This is the open
/ infinite-buffer regime that the closed-CTMC state-vector path cannot
represent, and corresponds to the joint-MMPP baseline of the fyp26 study.
"""

import numpy as np

from line_solver import (ClosedClass, CTMC, Delay, ENV, Exp, GlobalConstants,
                         MAM, MAP, MVA, Network, OpenClass, Queue,
                         SchedStrategy, Sink, Source, VerboseLevel)
from line_solver.environment import Environment


def terminal_model(stage_rate, crane_rate, n_cranes, N):
    """One hour-stage network: N carriers cycling yard-Delay -> quay-crane Queue."""
    qn = Network('Terminal')
    yard = Delay(qn, 'Yard')
    cranes = Queue(qn, 'QuayCranes', SchedStrategy.FCFS)
    cranes.setNumberOfServers(n_cranes)
    containers = ClosedClass(qn, 'Containers', N, yard)
    yard.setService(containers, Exp(stage_rate))
    cranes.setService(containers, Exp(crane_rate))
    qn.link(Network.serialRouting(yard, cranes))
    return qn


def terminal_fes_model(stage_rate, fes_rate):
    """One hour-stage with the internal terminal handling as a closed
    load-dependent FES: N containers cycling yard-Delay -> FES, where the FES
    serves at rate mu(n) = fes_rate[n-1] when n containers are inside it."""
    N = len(fes_rate)
    qn = Network('TerminalFES')
    yard = Delay(qn, 'Yard')
    fesq = Queue(qn, 'TerminalFES', SchedStrategy.PS)
    containers = ClosedClass(qn, 'Containers', N, yard)
    yard.setService(containers, Exp(stage_rate))
    fesq.setService(containers, Exp(fes_rate[0]))            # base rate mu(1)
    fesq.setLoadDependence(np.asarray(fes_rate) / fes_rate[0])  # alpha(n) = mu(n)/mu(1)
    qn.link(Network.serialRouting(yard, fesq))
    return qn


def fes_rate_curve(mu_quay, mu_stack, N):
    """Norton flow-equivalent rates: throughput of the isolated internal
    subnetwork (quay cranes -> stacking cranes, processor sharing) at n = 1..N
    containers, with the rest of the terminal short-circuited by a
    near-instantaneous delay."""
    mu = np.zeros(N)
    for n in range(1, N + 1):
        sub = Network('TerminalInternals')
        ref = Delay(sub, 'ShortCircuit')
        quay = Queue(sub, 'QuayCranes', SchedStrategy.PS)
        stack = Queue(sub, 'StackingCranes', SchedStrategy.PS)
        cls = ClosedClass(sub, 'Containers', n, ref)
        ref.setService(cls, Exp(1e6))                 # ~instantaneous short-circuit
        quay.setService(cls, Exp(mu_quay))
        stack.setService(cls, Exp(mu_stack))
        sub.link(Network.serialRouting(ref, quay, stack))
        T = MVA(sub, 'exact', verbose=False).getAvg()[3]
        mu[n - 1] = float(np.sum(T[0, :]))            # subnetwork throughput at n
    return mu


def exact_joint_metrics(env, ctmc_factory):
    """Day-averaged metrics from the exact joint (hour x network-state) CTMC."""
    from line_solver.api.mc import ctmc_solve
    from line_solver.api.solvers.ctmc.handler import solver_ctmc
    from line_solver.solvers.solver_env.statevec import _avg_from_pi

    renvQ, _ = ENV(env, ctmc_factory).getGenerator()
    pi_joint = np.asarray(ctmc_solve(renvQ), dtype=float).ravel()

    mdl = env.getEnsemble()
    E = len(mdl)
    M = mdl[0].getNumberOfStations()
    K = mdl[0].getNumberOfClasses()
    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    TN = np.zeros((M, K))
    off = 0
    for e in range(E):
        opts = ctmc_factory(mdl[e]).getOptions   # an attribute here, not a call
        ret = solver_ctmc(mdl[e].getStruct(), opts)
        ns = np.asarray(ret.infgen).shape[0]
        blk = pi_joint[off:off + ns]
        off += ns
        prob_env = float(np.sum(blk))
        if prob_env > 0:
            blk = blk / prob_env
        QNe, UNe, _, TNe = _avg_from_pi(ret.sn, blk, ret.space_aggr, ret.space_hashed,
                                        ret.arvRates, ret.depRates)
        QN += prob_env * np.asarray(QNe)
        UN += prob_env * np.asarray(UNe)
        TN += prob_env * np.asarray(TNe)
    return QN, UN, TN


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    # Daily demand schedule (Dhingra Rotterdam terminal, 24 hours)
    daily_rate_hr = np.array([6, 30, 40, 62, 76, 79, 119, 164, 152, 130, 79, 70,
                              57, 57, 113, 130, 162, 202, 148, 118, 92, 62, 36, 8],
                             dtype=float)
    daily_duration_hr = np.array([0.51, 0.84, 0.76, 0.98, 0.67, 2.33, 1.80, 0.62,
                                  0.36, 1.30, 1.02, 0.98, 0.86, 1.56, 1.30, 1.57,
                                  0.34, 0.22, 0.47, 1.33, 1.56, 0.93, 0.71, 1.11],
                                 dtype=float)
    E = daily_rate_hr.size                # 24 hourly environment stages

    # Terminal physical parameters (kept small so the CTMC is exact and fast)
    N = 6            # straddle carriers / containers circulating (closed tokens)
    n_cranes = 2     # quay cranes (multi-server FCFS queue)
    crane_rate = 16  # moves/hr served per crane
    stage_rate = daily_rate_hr / N        # per-token yard completion rate per hour

    # Build the random environment: one stage per hour, cyclic daily transitions
    env = Environment('RotterdamDailyCycle')
    for h in range(E):
        env.addStage('Hour%02d' % h, 'operational',
                     terminal_model(stage_rate[h], crane_rate, n_cranes, N))
    # Deterministic-on-average daily cycle h -> h+1 with exponential hour durations
    for h in range(E):
        env.addTransition('Hour%02d' % h, 'Hour%02d' % ((h + 1) % E),
                          Exp(1.0 / daily_duration_hr[h]))
    env.init()

    print('Rotterdam container terminal: %d hourly stages, %d carriers, %d cranes.'
          % (E, N, n_cranes))

    # Inner solver: CTMC with a finite transient horizon (required by statevec)
    T_HORIZON = 12   # hours; covers the hour-duration CDF tails (>99%)

    def ctmc_factory(m):
        return CTMC(m, 'exact', timespan=[0, T_HORIZON], verbose=False)

    base_opt = {'iter_max': 100, 'iter_tol': 1e-5, 'verbose': False}

    # (a) State-vector analyzer
    opt_statevec = dict(base_opt, method='statevec')
    Qsv, Usv, _, Tsv = ENV(env, ctmc_factory, opt_statevec).getAvg()[:4]

    # (b) Mean-field analyzer (default coupling)
    opt_meanfield = dict(base_opt, method='meanfield')
    Qmf, Umf, _, Tmf = ENV(env, ctmc_factory, opt_meanfield).getAvg()[:4]

    # (c) Exact joint random-environment CTMC (ground truth)
    Qex, Uex, Tex = exact_joint_metrics(env, ctmc_factory)

    crane = 1   # 0-based station index of the QuayCranes queue
    print('\n=== Day-averaged quay-crane metrics (container terminal) ===')
    print('%-10s %12s %12s %12s' % ('analyzer', 'QLen', 'Util', 'Tput'))
    for label, Q, U, T in (('exact', Qex, Uex, Tex), ('statevec', Qsv, Usv, Tsv),
                           ('meanfield', Qmf, Umf, Tmf)):
        print('%-10s %12.5f %12.5f %12.5f'
              % (label, np.sum(Q[crane, :]), np.sum(U[crane, :]), np.sum(T[crane, :])))

    print('\nQLen error vs exact:  statevec = %.3e , meanfield = %.3e'
          % (abs(np.sum(Qsv[crane, :] - Qex[crane, :])),
             abs(np.sum(Qmf[crane, :] - Qex[crane, :]))))
    print('Util error vs exact:  statevec = %.3e , meanfield = %.3e'
          % (abs(np.sum(Usv[crane, :] - Uex[crane, :])),
             abs(np.sum(Umf[crane, :] - Uex[crane, :]))))

    print('\nDay-averaged crane-queue table (state-vector analyzer):')
    print(ENV(env, ctmc_factory, opt_statevec).getAvgTable())

    # ===== Closed flow-equivalent-server (FES) variant =====================
    print('\n=== Closed FES variant (Norton-aggregated terminal internals) ===')
    mu_quay, mu_stack = 16.0, 20.0
    fes_rate = fes_rate_curve(mu_quay, mu_stack, N)
    print('Norton FES rate curve mu(n) = %s' % np.array2string(fes_rate, precision=5))

    env_fes = Environment('RotterdamDailyCycleFES')
    for h in range(E):
        env_fes.addStage('Hour%02d' % h, 'operational',
                         terminal_fes_model(stage_rate[h], fes_rate))
    for h in range(E):
        env_fes.addTransition('Hour%02d' % h, 'Hour%02d' % ((h + 1) % E),
                              Exp(1.0 / daily_duration_hr[h]))
    env_fes.init()

    Qsf, Usf, _, Tsf = ENV(env_fes, ctmc_factory, opt_statevec).getAvg()[:4]
    Qff, Uff, _, Tff = ENV(env_fes, ctmc_factory, opt_meanfield).getAvg()[:4]
    Qxf, Uxf, Txf = exact_joint_metrics(env_fes, ctmc_factory)

    fes = 1   # 0-based station index of the load-dependent FES
    print('%-10s %12s %12s %12s' % ('analyzer', 'QLen', 'Util', 'Tput'))
    for label, Q, U, T in (('exact', Qxf, Uxf, Txf), ('statevec', Qsf, Usf, Tsf),
                           ('meanfield', Qff, Uff, Tff)):
        print('%-10s %12.5f %12.5f %12.5f'
              % (label, np.sum(Q[fes, :]), np.sum(U[fes, :]), np.sum(T[fes, :])))
    print('\nQLen error vs exact:  statevec = %.3e , meanfield = %.3e'
          % (abs(np.sum(Qsf[fes, :] - Qxf[fes, :])),
             abs(np.sum(Qff[fes, :] - Qxf[fes, :]))))
    print('Util error vs exact:  statevec = %.3e , meanfield = %.3e'
          % (abs(np.sum(Usf[fes, :] - Uxf[fes, :])),
             abs(np.sum(Uff[fes, :] - Uxf[fes, :]))))

    # ===== Open system via MAM (MMPP(24)/M/c, infinite buffer) =============
    print('\n=== Open system via MAM (MMPP(24)/M/c, infinite buffer) ===')
    Qenv = np.zeros((E, E))
    for h in range(E):
        Qenv[h, (h + 1) % E] = 1.0 / daily_duration_hr[h]   # cyclic hour transition
        Qenv[h, h] = -1.0 / daily_duration_hr[h]
    mmpp = MAP(Qenv - np.diag(daily_rate_hr), np.diag(daily_rate_hr))
    c_open = 8                                              # cranes (open regime)

    open_model = Network('OpenTerminal')
    ships = Source(open_model, 'Ships')
    quay_open = Queue(open_model, 'QuayCranes', SchedStrategy.FCFS)
    done = Sink(open_model, 'Departures')
    oc = OpenClass(open_model, 'Containers')
    ships.setArrival(oc, mmpp)
    quay_open.setService(oc, Exp(crane_rate))
    quay_open.setNumberOfServers(c_open)
    open_model.link(Network.serialRouting(ships, quay_open, done))

    mean_lambda = float(np.sum(daily_duration_hr * daily_rate_hr) / np.sum(daily_duration_hr))
    print('MMPP mean arrival = %.2f/hr, %d cranes @ %.0f/hr, rho = %.3f'
          % (mean_lambda, c_open, crane_rate, mean_lambda / (c_open * crane_rate)))
    print(MAM(open_model).getAvgTable())
