"""Regression tests for the SolverFLD non-homogeneous Poisson (NHPP) transient.

A cyclic NHPP source intensity lambda(t) is injected into the closing fluid ODE
as a per-event rate multiplier (solver_fluid_ratemult, via
options.config.nhpp_sched set by SolverFLD.getTranAvg). With a fast server the
queue is never the bottleneck, so its throughput must track lambda(t) segment by
segment. The steady-state getAvg is deliberately unaffected and reports the
time-average rate.

The SOURCE throughput is a fluid open-model artifact and is not asserted on;
observe the QUEUE. Mirrors line-test.git/test/testsFLD/test_solver_fld_nhpp.m.
"""

import unittest

import numpy as np

from line_solver import (Network, Source, Queue, Sink, OpenClass, Exp, NHPP,
                         SchedStrategy, SolverFLD)


BREAKPOINTS = [0, 3, 4, 6]
RATES = [2, 8, 4]
TEND = 12.0  # two full periods


def _build(arrival):
    model = Network('fld_nhpp')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    jobclass = OpenClass(model, 'Class1')
    source.setArrival(jobclass, arrival)
    queue.setService(jobclass, Exp(50))  # fast server: queue tracks the input
    model.link(Network.serialRouting(source, queue, sink))
    return model


class TestSolverFLDNHPP(unittest.TestCase):

    def test_transient_tracks_intensity(self):
        nhpp = NHPP(BREAKPOINTS, RATES, True)
        model = _build(nhpp)
        solver = SolverFLD(model, timespan=[0, TEND], verbose=0)
        _, _, TNt = solver.getTranAvg()
        t = np.asarray(TNt[1][0].t).ravel()
        y = np.asarray(TNt[1][0].metric).ravel()

        # Mid-point of every segment of both periods, avoiding the transition
        # ramps that the step-faithful grid leaves at the breakpoints.
        mids = []
        for k in range(2):
            for j in range(len(RATES)):
                mids.append(k * nhpp.getPeriod()
                            + (BREAKPOINTS[j] + BREAKPOINTS[j + 1]) / 2.0)
        mids = [tt for tt in mids if 0.5 < tt < TEND]

        for tt in mids:
            expected = nhpp.getRateAt(tt)
            actual = float(np.interp(tt, t, y))
            self.assertAlmostEqual(actual / expected, 1.0, delta=0.05,
                                   msg='queue Tput %g at t=%g does not track '
                                       'lambda(t)=%g' % (actual, tt, expected))

        # The transient must actually vary: a solver that silently used the
        # time-average rate would give a flat trajectory.
        mask = t > 0.5
        spread = (y[mask].max() - y[mask].min()) / nhpp.getTimeAverageRate()
        self.assertGreater(spread, 0.5,
                           'trajectory is flat; lambda(t) was ignored')

    def test_steady_state_is_time_average(self):
        nhpp = NHPP(BREAKPOINTS, RATES, True)
        model = _build(nhpp)
        Q, U, R, T, A, W = SolverFLD(model, verbose=0).getAvg()
        lam_avg = nhpp.getTimeAverageRate()
        self.assertAlmostEqual(float(T[1, 0]) / lam_avg, 1.0, delta=0.05)

    def test_homogeneous_model_unperturbed(self):
        lam_avg = NHPP(BREAKPOINTS, RATES, True).getTimeAverageRate()
        model = _build(Exp(lam_avg))
        Q, U, R, T, A, W = SolverFLD(model, verbose=0).getAvg()
        self.assertAlmostEqual(float(T[1, 0]), lam_avg, places=6)
        # rho = lambda/mu is the FIRST-ORDER queue length, which is what the
        # 'matrix' method returns exactly. The default resolves to 'minnormal',
        # whose second-order closure adds the variance term and lands slightly
        # above it, so this pins the closure and not a Python-only value.
        #
        # RE-RECORDED 2026-09-01. The old 0.0733548838158009 was MATLAB's answer
        # before the alternation and its inner mean solve were tightened from
        # CoarseTol to mom_tol = 1e-6, which moves the converged value by 4.8e-8
        # -- past the 5e-10 that places=9 allows. MATLAB R2026a now reports
        # 0.07335493204019336 and native python 0.07335493204025656, i.e. the
        # two agree to 6.3e-14, which is what licenses the re-record.
        self.assertAlmostEqual(float(Q[1, 0]), 0.07335493204019336, places=9)
        Qm = SolverFLD(model, method='matrix', verbose=0).getAvg()[0]
        self.assertAlmostEqual(float(Qm[1, 0]), lam_avg / 50.0, places=6)

    def test_featureset_declares_nhpp(self):
        self.assertIn('NHPP', SolverFLD.getFeatureSet())


if __name__ == '__main__':
    unittest.main()
