"""
EXAMPLE_FLUID_MOMENTCLOSURE Second-order moment closures in SolverFLD.

The default fluid methods close the moment hierarchy at first order: the drift
of the mean uses min(E[X],c) in place of E[min(X,c)], so no second moment is
ever computed and the mean is biased where min() bends. This example runs the
three second-order methods against the exact CTMC on a closed two-station model
swept through saturation, where that bias peaks.

  'minnormal'  min-normal closure (Guenther, Stefanek, Bradley), mean and
               covariance solved self-consistently
  'refined'    O(1/N) refined mean field on the mean-field fixed point
"""

import numpy as np

from line_solver import (ClosedClass, CTMC, Delay, Exp, FLD, GlobalConstants,
                         Network, Queue, SchedStrategy, VerboseLevel)


def local_model(N):
    model = Network('momentclosure')
    delay = Delay(model, 'Think')
    queue = Queue(model, 'Server', SchedStrategy.PS)
    queue.setNumberOfServers(2)
    jobclass = ClosedClass(model, 'Class1', N, delay)
    delay.setService(jobclass, Exp(1.0))
    queue.setService(jobclass, Exp(1.0))
    model.link(Network.serialRouting(delay, queue))
    return model


def local_ld_model(N, alpha):
    model = Network('momentclosure_ld')
    delay = Delay(model, 'Think')
    queue = Queue(model, 'Server', SchedStrategy.PS)
    queue.setLoadDependence(np.asarray(alpha, dtype=float))
    jobclass = ClosedClass(model, 'Class1', N, delay)
    delay.setService(jobclass, Exp(1.0))
    queue.setService(jobclass, Exp(1.0))
    model.link(Network.serialRouting(delay, queue))
    return model


def _two_class_model(name, sched, w2):
    model = Network(name)
    delay = Delay(model, 'Think')
    queue = Queue(model, 'Server', sched)
    queue.setNumberOfServers(1)
    c1 = ClosedClass(model, 'Class1', 2, delay)
    c2 = ClosedClass(model, 'Class2', 2, delay)
    delay.setService(c1, Exp(1.0))
    delay.setService(c2, Exp(1.0))
    queue.setService(c1, Exp(1.0))
    queue.setService(c2, Exp(1.0))
    queue.setStrategyParam(c1, 1)
    queue.setStrategyParam(c2, w2)
    P = model.initRoutingMatrix()
    P.set(c1, c1, Network.serialRouting(delay, queue))
    P.set(c2, c2, Network.serialRouting(delay, queue))
    model.link(P)
    return model


def local_dps_model(w2):
    return _two_class_model('momentclosure_dps', SchedStrategy.DPS, w2)


def local_gps_model(w2):
    return _two_class_model('momentclosure_gps', SchedStrategy.GPS, w2)


def _qlen(solver_result, station):
    """Station-total mean queue length off a (station x class) matrix."""
    return float(np.sum(np.asarray(solver_result)[station, :]))


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    Nvals = [2, 4, 6, 8]
    methods = ['closing', 'minnormal', 'refined']

    print('Closed model: Delay(Z=1) -> Queue(PS, mu=1, c=2), sweeping population\n')
    print('%4s %10s %10s %10s %10s' % ('N', 'CTMC Q2', 'closing', 'minnormal', 'refined'))
    for N in Nvals:
        model = local_model(N)
        qexact = _qlen(CTMC(model).getAvg()[0], 1)
        row = [_qlen(FLD(model, method=m).getAvg()[0], 1) for m in methods]
        print('%4d %10.4f %10.4f %10.4f %10.4f' % (N, qexact, row[0], row[1], row[2]))

    # the covariance is only produced by the second-order methods; the reference
    # is the stationary distribution of the birth-death chain this model reduces
    # to, with birth rate (N-n) and death rate min(n,2)
    print('\nQueue-length standard deviation at the server (exact vs closures)')
    print('%4s %10s %10s %10s' % ('N', 'exact', 'minnormal', 'refined'))
    for N in Nvals:
        model = local_model(N)
        p = np.ones(N + 1)
        for n in range(1, N + 1):
            p[n] = p[n - 1] * (N - n + 1) / min(n, 2)
        p = p / p.sum()
        n = np.arange(N + 1)
        std_exact = float(np.sqrt(np.sum(p * n ** 2) - np.sum(p * n) ** 2))

        row = []
        for m in ('minnormal', 'refined'):
            solver = FLD(model, method=m)
            solver.getAvg()
            row.append(float(solver.getMoments()['QStd'][1, 0]))
        print('%4d %10.4f %10.4f %10.4f' % (N, std_exact, row[0], row[1]))

    # the same closure carries a limited load dependence: alpha(n) multiplies the
    # scheduling share, so the closed term becomes E[min(X,c)*alpha(X)]. Only the
    # closing family evaluates alpha; the other FLD methods refuse the model.
    alpha = [1.0, 1.7, 2.2, 2.5, 2.6, 2.65]
    print('\nLoad-dependent server, alpha = %s' % np.array2string(np.array(alpha)))
    print('%4s %10s %10s %10s %10s' % ('N', 'exact', 'closing', 'minnormal', 'refined'))
    for N in range(2, 7):
        model = local_ld_model(N, alpha)
        qexact = _qlen(CTMC(model).getAvg()[0], 1)
        row = [_qlen(FLD(model, method=m).getAvg()[0], 1)
               for m in ('closing', 'minnormal', 'refined')]
        print('%4d %10.4f %10.4f %10.4f %10.4f' % (N, qexact, row[0], row[1], row[2]))

    # min() is not the only non-linear rate term. The capacity share of a DPS
    # station, w_k*X_k/sum_j w_j*X_j, is a RATIO of populations, so evaluating it
    # at the mean is a second closure: it biases the split towards the class with
    # the larger weight while leaving the station total correct. The same
    # covariance closes it, so the per-class utilization improves without the
    # aggregate moving.
    print('\nDPS server, per-class utilization split (weights [1 w2])')
    print('%4s %21s %21s %21s' % ('w2', 'exact', 'closing', 'minnormal'))
    for w2 in (1, 2, 4, 8):
        model = local_dps_model(w2)
        ue = np.asarray(CTMC(model).getAvg()[1])[1, :]
        uc = np.asarray(FLD(model, method='closing').getAvg()[1])[1, :]
        ug = np.asarray(FLD(model, method='minnormal').getAvg()[1])[1, :]
        print('%4d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f'
              % (w2, ue[0], ue[1], uc[0], uc[1], ug[0], ug[1]))

    # GPS is the discipline where the second moment is not a correction but the
    # ENTIRE mechanism. GPS divides the server by weight among the BACKLOGGED
    # classes, so its share depends on the backlog INDICATOR, not on populations.
    # A first-order closure cannot express it: with continuous x_k > 0 every class
    # is always backlogged and the share collapses to the constant w_k/sum(w),
    # which is the heavy-traffic limit and is wrong at any other load. The closure
    # enumerates the 2^K backlog patterns weighted by P(N_k >= 1).
    print('\nGPS server, per-class utilization split (weights [1 w2])')
    print('%4s %21s %21s %21s' % ('w2', 'exact', 'minnormal', 'first-order const'))
    for w2 in (1, 2, 4, 8):
        model = local_gps_model(w2)
        ue = np.asarray(CTMC(model).getAvg()[1])[1, :]
        ug = np.asarray(FLD(model, method='minnormal').getAvg()[1])[1, :]
        uc = np.array([1.0, float(w2)]) / (1.0 + w2)
        print('%4d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f'
              % (w2, ue[0], ue[1], ug[0], ug[1], uc[0], uc[1]))
