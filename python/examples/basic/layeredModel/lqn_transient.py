#!/usr/bin/env python3
"""LQN_TRANSIENT  Transient (time-dependent) analysis of a layered network.

Demonstrates LN.getTranAvg, which returns the transient mean queue length,
utilization and throughput of every ensemble layer over time. The traces are
assembled block-diagonally: layer e occupies a disjoint block of rows (its
stations) and columns (its classes). Transient output is only produced by
transient-capable per-layer solvers (Fluid, CTMC, SSA); here the layers are
solved with the fluid ODE solver. Do not set a 'timespan' on the per-layer
factory: the steady-state fixed point rejects it, and getTranAvg auto-selects
the timespan per layer.

getTranAvg first converges the LN fixed point (getAvg), pinning the inter-layer
demands to equilibrium, then runs each layer's transient with those demands
frozen. The initial point is the layer's default state (all closed jobs at the
reference station, via State.initDefault), NOT the converged occupancy, so each
curve relaxes from all-at-reference to the layer steady state. Use setState to
seed a different start.

The model is a simple 3-layer LQN: a reference task T1 on processor P1 whose
activity synchronously calls entry E2 of task T2 on processor P2. The three
ensemble layers are the two processor (host) layers and the T2 task layer.
"""

import numpy as np

from line_solver import (Activity, Entry, Exp, FLD, GlobalConstants,
                         LayeredNetwork, LN, Processor, SchedStrategy, Task,
                         VerboseLevel)


def lqn_transient():
    model = LayeredNetwork('lqn_transient')

    P1 = Processor(model, 'P1', 1, SchedStrategy.PS)
    P2 = Processor(model, 'P2', 1, SchedStrategy.PS)

    T1 = Task(model, 'T1', 5, SchedStrategy.REF).on(P1).setThinkTime(Exp(1.0))
    T2 = Task(model, 'T2', 5, SchedStrategy.FCFS).on(P2).setThinkTime(Exp(1.0))

    E1 = Entry(model, 'E1').on(T1)
    E2 = Entry(model, 'E2').on(T2)

    Activity(model, 'A1', Exp(2.0)).on(T1).boundTo(E1).synchCall(E2, 1)
    Activity(model, 'A2', Exp(3.0)).on(T2).boundTo(E2).repliesTo(E2)

    return model


def _trace_arrays(trace):
    """(t, metric, station name, class name) of a transient trace entry, or None."""
    if trace is None:
        return None
    t = getattr(trace, 't', None)
    metric = getattr(trace, 'metric', None)
    if t is None and isinstance(trace, dict):
        t, metric = trace.get('t'), trace.get('metric')
    if t is None or metric is None:
        return None
    t = np.asarray(t, dtype=float).ravel()
    metric = np.asarray(metric, dtype=float).ravel()
    if metric.size == 0 or not np.any(metric > 1e-6):
        return None
    handle = getattr(trace, 'handle', None)
    if handle is None and isinstance(trace, dict):
        handle = trace.get('handle')
    try:
        stname = handle[0].getName()
        clname = handle[1].getName()
    except Exception:
        stname, clname = 'station', 'class'
    return t, metric, stname, clname


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    # Solve the ensemble with a fluid solver on each layer, then obtain the
    # per-layer transient averages in a single call.
    solver = LN(lqn_transient(), lambda m: FLD(m, verbose=False), verbose=False)
    QNt, UNt, TNt = solver.getTranAvg()[:3]

    E = solver.nlayers
    print('LN.getTranAvg returned transient traces for %d layers.' % E)

    # Plot the transient mean queue length E[N](t) of each station, one panel per
    # layer. The block-diagonal offsets (r0,c0) advance by each layer's station
    # and class counts, mirroring how getTranAvg stacks the per-layer blocks.
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        plt = None

    if plt is not None:
        fig, axes = plt.subplots(E, 1, figsize=(8, 3 * E), squeeze=False)
        fig.canvas.manager.set_window_title('LN.getTranAvg: per-layer transient E[N](t)')
    r0 = c0 = 0
    for e in range(E):
        sn = solver.ensemble[e].getStruct()
        M, K = int(sn.nstations), int(sn.nclasses)
        ax = axes[e][0] if plt is not None else None
        tsettle = 0.0
        for i in range(M):
            for r in range(K):
                got = _trace_arrays(QNt[r0 + i][c0 + r])
                if got is None:
                    continue
                t, metric, stname, clname = got
                if ax is not None:
                    ax.plot(t, metric, linewidth=1.5,
                            label='%s [%s]' % (stname, clname))
                # Track when this trace stops changing so the axis can be zoomed
                # onto the transient (getTranAvg auto-selects a wide timespan and
                # the curves are flat once settled).
                if metric.size > 1:
                    tol = 0.01 * (metric.max() - metric.min()) + 1e-9
                    moving = np.flatnonzero(np.abs(metric - metric[-1]) > tol)
                    if moving.size:
                        klast = int(moving[-1])
                        tsettle = max(tsettle, float(t[min(klast + 1, metric.size - 1)]))
        r0 += M
        c0 += K
        if ax is not None:
            if tsettle > 0:
                ax.set_xlim(0, 1.3 * tsettle)
            ax.set_title('Layer %d' % (e + 1))
            ax.set_xlabel('time t')
            ax.set_ylabel('E[N](t)')
            ax.grid(True)
            ax.legend(loc='best')
    if plt is not None:
        fig.tight_layout()
        fig.savefig('lqn_transient.png', dpi=110)
        print('Wrote lqn_transient.png')
