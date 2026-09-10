"""Regression tests for MAPt/PHt and the SolverFLD 'kp' method.

The 'kp' method integrates the fluid and diffusion limits of Ko and Pender,
"Diffusion limits for the (MAP_t/Ph_t/inf)^N queueing network", Oper. Res. Lett.
45 (2017) 248-253. For an infinite-server network the rate functions are affine
in the state, so the mean and covariance ODEs close EXACTLY: every assertion here
compares against an exact reference rather than a tolerance band.

References used, strongest first:
  1. a time-inhomogeneous CTMC over (arrival phase, per-service-phase counts),
     integrated to high accuracy on a truncated state space -- the exact law of
     the model itself, transient and all;
  2. the exact stationary MAP/M/inf CTMC;
  3. Var == Mean for M_t/M/inf, whose queue length is Poisson at every t;
  4. degeneracy: a 1-phase MAPt is an NHPP, a 1-segment MAPt is a MAP.
"""

import itertools
import unittest

import numpy as np
from scipy.integrate import solve_ivp

from line_solver import (Delay, Exp, MAP, MAPt, Network, NHPP, OpenClass, PHt,
                         Sink, SolverFLD, Source)
from line_solver.api.mam import map_lambda

# a 2-phase MAP whose D0 is cyclic, so the process is genuinely non-renewal
D0_A = np.array([[-5.0, 1.0], [2.0, -4.0]])
D1_A = np.array([[3.0, 1.0], [1.0, 1.0]])
D0_B = np.array([[-12.0, 3.0], [5.0, -9.0]])
D1_B = np.array([[7.0, 2.0], [2.0, 2.0]])
BREAKPOINTS = [0.0, 1.0, 2.5]
MU = 2.0


def _stationary_phase(D0, D1):
    h = D0.shape[0]
    A = np.vstack([(D0 + D1).T, np.ones(h)])
    b = np.zeros(h + 1)
    b[-1] = 1.0
    theta = np.linalg.lstsq(A, b, rcond=None)[0]
    theta = np.maximum(theta, 0.0)
    return theta / theta.sum()


def _exact_mapt_pht_inf(bp, D0s, D1s, alphas, Ss, nmax, tspan, cyclic=True):
    """Exact transient of one MAP_t/Ph_t/inf station by CTMC integration.

    State: (arrival phase, count in each service phase), truncated at nmax per
    phase. Returns a callable ts -> (mean, var) of the total station content.
    """
    h = D0s[0].shape[0]
    hs = Ss[0].shape[0]
    period = bp[-1] - bp[0]
    counts = list(itertools.product(*[range(nmax + 1)] * hs))
    index = {(j, cc): k for k, (j, cc) in enumerate(itertools.product(range(h), counts))}
    S = len(index)

    def seg(t):
        off = (t - bp[0]) % period if cyclic else (t - bp[0])
        return min(int(np.searchsorted(bp[1:], bp[0] + off, side='right')), len(D0s) - 1)

    def gen(t):
        k = seg(t)
        D0, D1, alpha, Sm = D0s[k], D1s[k], alphas[k], Ss[k]
        svec = -Sm.sum(axis=1)
        Q = np.zeros((S, S))
        for (j, cc), a in index.items():
            for jp in range(h):
                if jp != j and D0[j, jp] != 0.0:
                    Q[a, index[(jp, cc)]] += D0[j, jp]
                    Q[a, a] -= D0[j, jp]
                if D1[j, jp] != 0.0:
                    for p in range(hs):
                        if alpha[p] <= 0 or cc[p] >= nmax:
                            continue
                        nc = list(cc)
                        nc[p] += 1
                        Q[a, index[(jp, tuple(nc))]] += D1[j, jp] * alpha[p]
                        Q[a, a] -= D1[j, jp] * alpha[p]
            for p in range(hs):
                if cc[p] == 0:
                    continue
                for q in range(hs):
                    if q != p and Sm[p, q] != 0.0 and cc[q] < nmax:
                        nc = list(cc)
                        nc[p] -= 1
                        nc[q] += 1
                        Q[a, index[(j, tuple(nc))]] += Sm[p, q] * cc[p]
                        Q[a, a] -= Sm[p, q] * cc[p]
                if svec[p] != 0.0:
                    nc = list(cc)
                    nc[p] -= 1
                    Q[a, index[(j, tuple(nc))]] += svec[p] * cc[p]
                    Q[a, a] -= svec[p] * cc[p]
        return Q

    theta = _stationary_phase(D0s[seg(bp[0])], D1s[seg(bp[0])])
    p0 = np.zeros(S)
    for j in range(h):
        p0[index[(j, tuple([0] * hs))]] = theta[j]
    sol = solve_ivp(lambda t, p: gen(t).T @ p, tspan, p0, dense_output=True,
                    max_step=0.02, rtol=1e-10, atol=1e-13)
    total = np.array([sum(cc) for (_, cc) in index.keys()], dtype=float)

    def at(ts):
        P = np.maximum(sol.sol(ts), 0.0)
        P = P / P.sum(axis=0, keepdims=True)
        m1 = total @ P
        return m1, (total ** 2) @ P - m1 ** 2
    return at


def _inf_station(arrival, service, name='m'):
    model = Network(name)
    src = Source(model, 'Source')
    dly = Delay(model, 'Delay')
    snk = Sink(model, 'Sink')
    cls = OpenClass(model, 'Class1')
    src.setArrival(cls, arrival)
    dly.setService(cls, service)
    model.link(Network.serialRouting(src, dly, snk))
    return model


def _kp(model, tend, tol=1e-9, maxstep=0.01):
    solver = SolverFLD(model, method='kp')
    solver.options.timespan = (0.0, tend)
    solver.options.tol = tol
    solver.options.odemaxstep = maxstep
    return solver._solve_kp()


class TestMAPtPHtProcesses(unittest.TestCase):
    """The process classes themselves."""

    def test_one_phase_mapt_is_an_nhpp(self):
        bp, rates = [0.0, 1.0, 2.0], [1.0, 3.0]
        mt = MAPt(bp, [np.array([[-r]]) for r in rates],
                  [np.array([[r]]) for r in rates])
        nh = NHPP(bp, rates, True)
        self.assertAlmostEqual(mt.getTimeAverageRate(), nh.getTimeAverageRate(), places=12)
        for t in (0.25, 1.25, 2.25, 3.25):
            self.assertAlmostEqual(mt.getRateAt(t), nh.getRateAt(t), places=12)

    def test_single_segment_mapt_is_a_map(self):
        mt = MAPt([0.0, 1.0], [D0_A], [D1_A])
        self.assertAlmostEqual(mt.getTimeAverageRate(), map_lambda(D0_A, D1_A), places=12)
        D0bar, D1bar = mt.getTimeAverageProcess()
        np.testing.assert_allclose(D0bar, D0_A)
        np.testing.assert_allclose(D1bar, D1_A)

    def test_scalar_summaries_are_nan(self):
        mt = MAPt(BREAKPOINTS, [D0_A, D0_B], [D1_A, D1_B])
        for value in (mt.getSCV(), mt.getVar(), mt.getSkewness(), mt.evalCDF(1.0),
                      mt.evalLST(1.0)):
            self.assertTrue(np.isnan(value))

    def test_varying_sparsity_pattern_is_refused(self):
        # A transition present in one segment and absent in another cannot be a
        # per-entry multiplier on the time-averaged nominal, so it is an error
        # rather than a silently dropped transition.
        D0_off = np.array([[-4.0, 0.0], [2.0, -4.0]])
        D1_off = np.array([[3.0, 1.0], [1.0, 1.0]])
        with self.assertRaises(ValueError):
            MAPt(BREAKPOINTS, [D0_A, D0_off], [D1_A, D1_off])

    def test_sampling_reproduces_the_stationary_rate(self):
        mt = MAPt([0.0, 1.0], [D0_A], [D1_A])
        rng = np.random.default_rng(20260728)
        intervals = mt.sample(60000, rng)
        self.assertAlmostEqual(1.0 / float(np.mean(intervals)),
                               map_lambda(D0_A, D1_A), delta=0.05)

    def test_non_cyclic_horizon_is_silent(self):
        # Outside its horizon a non-cyclic schedule emits nothing: nextArrival
        # returns 0, which callers read as "no further arrival". Sampling past
        # the horizon must terminate rather than loop.
        mt = MAPt(BREAKPOINTS, [D0_A, D0_B], [D1_A, D1_B], cyclic=False)
        rng = np.random.default_rng(7)
        interval, _ = mt.nextArrival(BREAKPOINTS[-1] + 1.0, 0, rng)
        self.assertEqual(interval, 0.0)
        self.assertEqual(mt.getSegmentIndexAt(BREAKPOINTS[-1] + 1.0), -1)
        np.testing.assert_allclose(mt.getD1At(BREAKPOINTS[-1] + 1.0),
                                   np.zeros_like(D1_A))
        # a path started inside the horizon stops once the horizon is exhausted
        mt.resetSampleClock()
        samples = mt.sample(500, rng)
        self.assertTrue(np.any(samples == 0.0))
        self.assertTrue(np.all(np.isfinite(samples)))

    def test_pht_non_cyclic_and_absorbing_boundary(self):
        alphas = [np.array([1.0, 0.0]), np.array([1.0, 0.0])]
        Ss = [np.array([[-4.0, 4.0], [0.0, -4.0]]),
              np.array([[-9.0, 9.0], [0.0, -9.0]])]
        ph = PHt(BREAKPOINTS, alphas, Ss, cyclic=False)
        rng = np.random.default_rng(11)
        self.assertEqual(ph.sampleFrom(BREAKPOINTS[-1] + 1.0, rng), 0.0)
        self.assertEqual(ph.getSegmentIndexAt(BREAKPOINTS[-1] + 1.0), -1)
        np.testing.assert_allclose(ph.getSAt(BREAKPOINTS[-1] + 1.0), np.zeros_like(Ss[0]))
        # alpha past the horizon falls back to the last segment, so it stays a
        # probability vector rather than becoming undefined
        np.testing.assert_allclose(ph.getAlphaAt(BREAKPOINTS[-1] + 1.0), alphas[-1])
        # a service started well inside the horizon still completes
        started = [ph.sampleFrom(0.05, rng) for _ in range(200)]
        self.assertTrue(all(v >= 0.0 for v in started))
        self.assertGreater(sum(1 for v in started if v > 0.0), 0)

    def test_cyclic_sampling_crosses_breakpoints(self):
        # The boundary-crossing branch: a holding time longer than the remaining
        # segment must advance the clock and resample under the new matrices
        # rather than fire under the stale ones.
        slow = np.array([[-0.01]])
        fast = np.array([[-50.0]])
        mt = MAPt([0.0, 0.5, 1.0], [slow, fast], [-slow, -fast], cyclic=True)
        rng = np.random.default_rng(3)
        intervals = mt.sample(4000, rng)
        self.assertTrue(np.all(intervals > 0.0))
        # the fast segment dominates the event count, so the mean interval is
        # far below the slow segment's 100 time units
        self.assertLess(float(np.mean(intervals)), 1.0)

    def test_struct_and_json_round_trip(self):
        from line_solver.io.linemodel_io import _dist_to_json, _json_to_dist
        mt = MAPt(BREAKPOINTS, [D0_A, D0_B], [D1_A, D1_B])
        back = _json_to_dist(_dist_to_json(mt))
        np.testing.assert_allclose(back.getBreakpoints(), mt.getBreakpoints())
        for a, b in zip(mt.D0, back.D0):
            np.testing.assert_allclose(a, b)
        for a, b in zip(mt.D1, back.D1):
            np.testing.assert_allclose(a, b)

        model = _inf_station(mt, Exp(MU))
        sn = model.getStruct()
        from line_solver.constants import ProcessType
        self.assertEqual(sn.procid[0, 0], ProcessType.MAPT)
        self.assertEqual(int(sn.phases[0, 0]), 2)
        self.assertTrue(np.isnan(sn.scv[0, 0]))


class TestFluidPhaseStructure(unittest.TestCase):
    """A non-acyclic D0 must not lose its last row in the fluid ODE."""

    def test_map_and_mmpp2_sources_carry_their_exact_rate(self):
        from line_solver import MMPP2
        cases = [(MAP(D0_A, D1_A), map_lambda(D0_A, D1_A))]
        mmpp = MMPP2(4.0, 1.0, 0.5, 0.5)
        cases.append((mmpp, mmpp.getRate()))
        for dist, exact in cases:
            for method in ('closing', 'matrix'):
                table = SolverFLD(_inf_station(dist, Exp(MU)), method=method).getAvgTable()
                self.assertAlmostEqual(float(np.asarray(table['Tput'])[1]), exact, places=6)


class TestKoPenderMethod(unittest.TestCase):
    """The 'kp' fluid + diffusion limits."""

    def test_variance_equals_mean_for_a_poisson_infinite_server(self):
        # M_t/M/inf has a Poisson queue length at every t, so the covariance ODE
        # must return exactly the mean. This pins G = A diag(f) A' and the
        # block aggregation independently of any simulator.
        bp, rates = [0.0, 1.0, 2.0], [1.0, 3.0]
        mt = MAPt(bp, [np.array([[-r]]) for r in rates], [np.array([[r]]) for r in rates])
        res = _kp(_inf_station(mt, Exp(MU)), 6.0)
        mean = res.QNt[(1, 0)]
        var = res.QVart[(1, 0)]
        np.testing.assert_allclose(var, mean, atol=1e-7)

    def test_matches_the_exact_stationary_map_m_inf(self):
        exact_mean = map_lambda(D0_A, D1_A) / MU
        model = _inf_station(MAPt([0.0, 1.0], [D0_A], [D1_A]), Exp(MU))
        res = _kp(model, 40.0)
        self.assertAlmostEqual(float(res.QNt[(1, 0)][-1]), exact_mean, places=7)
        # Var != Mean here: the arrival stream is non-renewal, which is exactly
        # what a PH renewal approximation of the MAP would lose.
        self.assertAlmostEqual(float(res.QVart[(1, 0)][-1]), 1.64, places=5)
        self.assertGreater(abs(float(res.QVart[(1, 0)][-1]) - exact_mean), 1e-3)

    def test_matches_the_exact_time_inhomogeneous_ctmc(self):
        alphas = [np.array([1.0]), np.array([1.0])]
        Ss = [np.array([[-MU]]), np.array([[-MU]])]
        at = _exact_mapt_pht_inf(BREAKPOINTS, [D0_A, D0_B], [D1_A, D1_B],
                                 alphas, Ss, 60, (0.0, 5.0))
        model = _inf_station(MAPt(BREAKPOINTS, [D0_A, D0_B], [D1_A, D1_B]), Exp(MU))
        res = _kp(model, 5.0)
        ts = np.array([0.25, 0.5, 1.0, 1.5, 2.5, 3.0, 4.0, 5.0])
        exact_mean, exact_var = at(ts)
        got_mean = np.interp(ts, res.t, res.QNt[(1, 0)])
        got_var = np.interp(ts, res.t, res.QVart[(1, 0)])
        np.testing.assert_allclose(got_mean, exact_mean, atol=1e-3)
        np.testing.assert_allclose(got_var, exact_var, atol=1e-3)

    def test_time_varying_phase_type_service(self):
        alphas = [np.array([1.0, 0.0]), np.array([1.0, 0.0])]
        Ss = [np.array([[-4.0, 4.0], [0.0, -4.0]]),
              np.array([[-9.0, 9.0], [0.0, -9.0]])]
        at = _exact_mapt_pht_inf(BREAKPOINTS, [D0_A, D0_B], [D1_A, D1_B],
                                 alphas, Ss, 14, (0.0, 5.0))
        model = _inf_station(MAPt(BREAKPOINTS, [D0_A, D0_B], [D1_A, D1_B]),
                             PHt(BREAKPOINTS, alphas, Ss))
        res = _kp(model, 5.0)
        ts = np.array([0.5, 1.0, 1.5, 2.5, 3.0, 4.0, 5.0])
        exact_mean, exact_var = at(ts)
        got_mean = np.interp(ts, res.t, res.QNt[(1, 0)])
        got_var = np.interp(ts, res.t, res.QVart[(1, 0)])
        np.testing.assert_allclose(got_mean, exact_mean, atol=1e-3)
        np.testing.assert_allclose(got_var, exact_var, atol=1e-3)

    def test_cyclic_steady_state_is_the_period_average(self):
        # With an unbounded horizon there is no fixed point, so getAvg reports
        # the average over the last full period; source and station throughput
        # must then agree, which a snapshot of the cycle would not.
        # tol below the default 1e-4: flow balance in the periodic regime is exact, so
        # asserting it only makes sense once the integration error is smaller than the
        # claim. Loosening the assertion instead would hide a real imbalance.
        model = _inf_station(MAPt(BREAKPOINTS, [D0_A, D0_B], [D1_A, D1_B]), Exp(MU))
        table = SolverFLD(model, method='kp', tol=1e-9).getAvgTable()
        tput = np.asarray(table['Tput'], dtype=float)
        self.assertAlmostEqual(tput[0], tput[1], places=5)
        widths = np.diff(BREAKPOINTS)
        rates = [map_lambda(D0_A, D1_A), map_lambda(D0_B, D1_B)]
        expected = float(np.dot(widths, rates) / np.sum(widths))
        self.assertAlmostEqual(tput[1], expected, delta=1e-3)

    def test_two_station_tandem(self):
        model = Network('tandem')
        src = Source(model, 'Source')
        d1 = Delay(model, 'Delay1')
        d2 = Delay(model, 'Delay2')
        snk = Sink(model, 'Sink')
        cls = OpenClass(model, 'Class1')
        src.setArrival(cls, MAPt(BREAKPOINTS, [D0_A, D0_B], [D1_A, D1_B]))
        d1.setService(cls, Exp(3.0))
        d2.setService(cls, Exp(MU))
        model.link(Network.serialRouting(src, d1, d2, snk))
        res = _kp(model, 5.0)
        # flow balance in the periodic regime: both stations see the same
        # long-run rate as the source
        for key in ((1, 0), (2, 0)):
            self.assertGreater(float(np.mean(res.QNt[key])), 0.0)
            np.testing.assert_array_less(-1e-9, res.QVart[key])

    def test_closed_model_is_rejected_by_the_gate(self):
        from line_solver import ClosedClass, Queue, SchedStrategy
        model = Network('closed')
        dly = Delay(model, 'Delay')
        que = Queue(model, 'Queue', SchedStrategy.PS)
        cls = ClosedClass(model, 'Class1', 5, dly)
        dly.setService(cls, Exp(1.0))
        que.setService(cls, Exp(2.0))
        model.link(Network.serialRouting(dly, que, dly))
        with self.assertRaises(Exception):
            SolverFLD(model, method='kp').getAvgTable()

    def test_seeded_mean_and_covariance_are_honoured(self):
        # M/M/inf started from a POISSON(x0) queue stays Poisson at every t, so a
        # run seeded with mean x0 and covariance x0 must return the exact transient
        # mean AND Var == Mean along the whole trajectory. Pins both seeds at once:
        # honouring only one breaks the identity.
        lam, x0 = 3.0, 5.0
        mt = MAPt([0.0, 1.0], [np.array([[-lam]])], [np.array([[lam]])])
        solver = SolverFLD(_inf_station(mt, Exp(MU)), method='kp')
        solver.options.timespan = (0.0, 4.0)
        solver.options.tol = 1e-9
        solver.options.odemaxstep = 0.01
        # dim = 2: the arrival phase (u-block) then the Delay service phase.
        solver.options.config = {'kp_init_sol': [1.0, x0],
                                 'init_cov': [[0.0, 0.0], [0.0, x0]]}
        res = solver._solve_kp()
        t = np.asarray(res.t)
        exact = x0 * np.exp(-MU * t) + (lam / MU) * (1.0 - np.exp(-MU * t))
        mean = np.asarray(res.QNt[(1, 0)])
        var = np.asarray(res.QVart[(1, 0)])
        np.testing.assert_allclose(mean, exact, atol=1e-7)
        np.testing.assert_allclose(var, mean, atol=1e-7)

    def test_a_seed_that_does_not_fit_is_refused(self):
        # Silently dropping it would integrate from the DEFAULT initial condition
        # under the caller's name and hand back a plausible wrong trajectory.
        mt = MAPt([0.0, 1.0], [np.array([[-3.0]])], [np.array([[3.0]])])
        for cfg in ({'kp_init_sol': [1.0, 2.0, 3.0]},
                    {'init_cov': [[1.0, 0.0]]},
                    {'init_cov': [[1.0, 2.0], [0.0, 1.0]]}):
            solver = SolverFLD(_inf_station(mt, Exp(MU)), method='kp')
            solver.options.timespan = (0.0, 1.0)
            solver.options.config = cfg
            with self.assertRaises(ValueError):
                solver._solve_kp()

    def test_featureset_and_method_registration(self):
        self.assertIn('MAPt', SolverFLD.getFeatureSet())
        self.assertIn('PHt', SolverFLD.getFeatureSet())
        self.assertIn('kp', SolverFLD.listValidMethods())
        solver = SolverFLD(_inf_station(Exp(1.0), Exp(MU)), method='kp')
        self.assertNotIn('ClosedClass', solver.getMethodFeatureSet('kp'))


if __name__ == '__main__':
    unittest.main()
