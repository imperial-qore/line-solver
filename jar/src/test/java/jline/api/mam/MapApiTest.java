/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validation of the MAP calculus APIs (jline.api.mam) against closed forms:
 * - the Poisson process as a 1-phase MAP has mean 1/lambda, SCV 1, zero
 *   autocorrelation and unit index of dispersion;
 * - map_erlang / map_hyperexp constructors must reproduce the requested
 *   moments exactly;
 * - structural transforms (scale, super, sum, normalize) obey exact
 *   composition rules for rates and moments.
 */
public class MapApiTest {

    private static final double TOL = 1e-9;

    private static MatrixCell poissonMap(double lambda) {
        Matrix d0 = new Matrix(1, 1);
        d0.set(0, 0, -lambda);
        Matrix d1 = new Matrix(1, 1);
        d1.set(0, 0, lambda);
        return new MatrixCell(d0, d1);
    }

    /** Two-phase MMPP-like MAP used for structural checks. */
    private static MatrixCell testMap() {
        Matrix d0 = new Matrix(2, 2);
        d0.set(0, 0, -4.0); d0.set(0, 1, 1.0);
        d0.set(1, 0, 0.5);  d0.set(1, 1, -2.0);
        Matrix d1 = new Matrix(2, 2);
        d1.set(0, 0, 2.0); d1.set(0, 1, 1.0);
        d1.set(1, 0, 1.0); d1.set(1, 1, 0.5);
        return new MatrixCell(d0, d1);
    }

    @Test
    public void poissonMapHasExactClosedForms() {
        double lambda = 2.0;
        MatrixCell map = poissonMap(lambda);
        assertEquals(lambda, Map_lambda.map_lambda(map), TOL, "Poisson MAP rate");
        assertEquals(1.0 / lambda, Map_mean.map_mean(map), TOL, "Poisson MAP mean");
        assertEquals(1.0, Map_scv.map_scv(map), TOL, "Poisson MAP SCV");
        assertEquals(1.0, Map_idc.map_idc(map), TOL, "Poisson MAP index of dispersion");
        Matrix acf = Map_acf.map_acf(map);
        for (int i = 0; i < acf.getNumElements(); i++) {
            assertEquals(0.0, acf.get(i), 1e-8, "Poisson MAP autocorrelation lag " + (i + 1));
        }
        assertTrue(Map_isfeasible.map_isfeasible(map), "Poisson MAP must be feasible");
    }

    @Test
    public void erlangConstructorMatchesMoments() {
        double mean = 0.5;
        int k = 3;
        MatrixCell erl = Map_erlang.map_erlang(mean, k);
        assertEquals(mean, Map_mean.map_mean(erl), TOL, "Erlang MAP mean");
        assertEquals(1.0 / k, Map_scv.map_scv(erl), TOL, "Erlang MAP SCV = 1/k");
        assertTrue(Map_isfeasible.map_isfeasible(erl), "Erlang MAP must be feasible");
    }

    @Test
    public void hyperexpConstructorMatchesMoments() {
        // p=0.99 is the canonical branch-probability used by the fitters;
        // small p values are infeasible for this two-moment construction
        double mean = 1.0, scv = 4.0, p = 0.99;
        MatrixCell he = Map_hyperexp.map_hyperexp(mean, scv, p);
        assertEquals(mean, Map_mean.map_mean(he), 1e-6, "HyperExp MAP mean");
        assertEquals(scv, Map_scv.map_scv(he), 1e-6, "HyperExp MAP SCV");
        assertTrue(Map_isfeasible.map_isfeasible(he), "HyperExp MAP must be feasible");
    }

    @Test
    public void scaleSetsMeanExactly() {
        MatrixCell scaled = Map_scale.map_scale(testMap(), 2.5);
        assertEquals(2.5, Map_mean.map_mean(scaled), TOL, "scaled MAP mean");
        // Scaling preserves the coefficient of variation
        assertEquals(Map_scv.map_scv(testMap()), Map_scv.map_scv(scaled), TOL,
                "scaling must preserve SCV");
    }

    @Test
    public void superpositionAddsRates() {
        MatrixCell a = poissonMap(2.0);
        MatrixCell b = poissonMap(3.0);
        MatrixCell sup = Map_super.map_super(a, b);
        assertEquals(5.0, Map_lambda.map_lambda(sup), TOL,
                "superposition of Poisson streams adds rates");
        // Superposition of independent Poisson processes is Poisson
        assertEquals(1.0, Map_scv.map_scv(sup), TOL,
                "superposed Poisson must remain SCV 1");
    }

    @Test
    public void sumOfMapsAddsMeans() {
        MatrixCell erl = Map_sum.map_sum(poissonMap(2.0), 3);
        // Sum of 3 iid Exp(2) interarrivals = Erlang-3: mean 1.5, SCV 1/3
        assertEquals(1.5, Map_mean.map_mean(erl), TOL, "3-fold sum mean");
        assertEquals(1.0 / 3.0, Map_scv.map_scv(erl), TOL, "3-fold sum SCV");
    }

    @Test
    public void embeddedChainAndPieAreConsistent() {
        MatrixCell map = testMap();
        Matrix pie = Map_pie.map_pie(map);
        Matrix embedded = Map_embedded.map_embedded(map);
        // pie is the stationary vector of the embedded chain: pie P = pie
        double sum = 0.0;
        for (int i = 0; i < pie.getNumElements(); i++) {
            assertTrue(pie.get(i) >= -TOL, "pie entries must be nonnegative");
            sum += pie.get(i);
        }
        assertEquals(1.0, sum, 1e-8, "pie must be a probability vector");
        for (int j = 0; j < 2; j++) {
            double v = 0.0;
            for (int i = 0; i < 2; i++) {
                v += pie.get(i) * embedded.get(i, j);
            }
            assertEquals(pie.get(j), v, 1e-8, "pie must be stationary for embedded P");
        }
    }

    @Test
    public void infgenIsAProperGenerator() {
        Matrix q = Map_infgen.map_infgen(testMap().get(0), testMap().get(1));
        for (int i = 0; i < 2; i++) {
            double rowSum = 0.0;
            for (int j = 0; j < 2; j++) {
                rowSum += q.get(i, j);
                if (i != j) {
                    assertTrue(q.get(i, j) >= 0, "infgen off-diagonal must be nonnegative");
                }
            }
            assertEquals(0.0, rowSum, TOL, "infgen rows must sum to zero");
        }
    }

    @Test
    public void cdfIsMonotoneAndReachesOne() {
        MatrixCell map = poissonMap(2.0);
        Matrix points = new Matrix(1, 5);
        double[] ts = {0.1, 0.5, 1.0, 2.0, 5.0};
        for (int i = 0; i < ts.length; i++) {
            points.set(0, i, ts[i]);
        }
        Matrix cdf = Map_cdf.map_cdf(map.get(0), map.get(1), points);
        double prev = 0.0;
        for (int i = 0; i < ts.length; i++) {
            double v = cdf.get(i);
            // Exponential CDF closed form
            assertEquals(1.0 - Math.exp(-2.0 * ts[i]), v, 1e-8,
                    "Poisson MAP interarrival CDF at t=" + ts[i]);
            assertTrue(v >= prev, "CDF must be nondecreasing");
            prev = v;
        }
    }

    @Test
    public void normalizationRepairsPerturbedMap() {
        Matrix d0 = testMap().get(0);
        Matrix d1 = testMap().get(1);
        d0.set(0, 0, d0.get(0, 0) + 1e-9); // slight infeasibility
        MatrixCell fixed = Map_normalize.map_normalize(d0, d1);
        assertTrue(Map_isfeasible.map_isfeasible(fixed),
                "normalized MAP must be feasible");
    }

    @Test
    public void mixtureCombinesRatesConvexly() {
        MatrixCell a = poissonMap(2.0);
        MatrixCell b = poissonMap(4.0);
        MatrixCell mix = Map_mixture.map_mixture(new double[]{0.5, 0.5},
                new MatrixCell[]{a, b});
        // Mixture mean: 0.5*(1/2) + 0.5*(1/4) = 0.375
        assertEquals(0.375, Map_mean.map_mean(mix), 1e-8,
                "mixture mean is the convex combination of means");
    }

    @Test
    public void timeReversePreservesMarginalStatistics() {
        // Time reversal of a MAP keeps the marginal interarrival distribution
        // (same mean, SCV, rate); only the correlation structure is reversed.
        MatrixCell map = testMap();
        MatrixCell rev = Map_timereverse.map_timereverse(map);
        assertTrue(Map_isfeasible.map_isfeasible(rev), "reversed MAP must be feasible");
        assertEquals(Map_mean.map_mean(map), Map_mean.map_mean(rev), 1e-9,
                "time reversal must preserve the mean");
        assertEquals(Map_scv.map_scv(map), Map_scv.map_scv(rev), 1e-9,
                "time reversal must preserve the SCV");
        assertEquals(Map_lambda.map_lambda(map), Map_lambda.map_lambda(rev), 1e-9,
                "time reversal must preserve the arrival rate");
    }

    @Test
    public void normalizeRepairsWithoutChangingRate() {
        MatrixCell map = testMap();
        double rate = Map_lambda.map_lambda(map);
        MatrixCell norm = Map_normalize.map_normalize(map.get(0), map.get(1));
        assertTrue(Map_isfeasible.map_isfeasible(norm), "normalized MAP must be feasible");
        assertEquals(rate, Map_lambda.map_lambda(norm), 1e-9,
                "normalization must preserve the arrival rate");
    }

    @Test
    public void indexOfDispersionIsOneForPoissonAboveOneForCorrelated() {
        assertEquals(1.0, Map_idc.map_idc(poissonMap(2.0)), 1e-9,
                "Poisson process has index of dispersion 1");
        assertTrue(Map_idc.map_idc(testMap()) > 1.0 + 1e-6,
                "a positively correlated MAP has index of dispersion > 1");
    }

    @Test
    public void gammaIsTheAcfDecayRateNotTheSecondEigenvalue() {
        // order 1: a Poisson process is uncorrelated
        assertEquals(0.0, Map_gamma.map_gamma(poissonMap(2.0)), TOL,
                "Poisson process has no autocorrelation decay");

        // order 2 with real correlation: the ACF is geometric, so gamma is
        // acf(2)/acf(1) exactly, and there it coincides with the second
        // eigenvalue of the embedded DTMC
        MatrixCell m2 = Map_mmpp2.map_mmpp2(1.0, 2.0, -1.0, 0.2);
        double acf1 = Map_acf.map_acf(m2.get(0), m2.get(1), 1).value();
        double acf2 = Map_acf.map_acf(m2.get(0), m2.get(1), 2).value();
        assertEquals(acf2 / acf1, Map_gamma.map_gamma(m2), TOL,
                "for order 2 gamma must be acf(2)/acf(1)");
        assertEquals(Map_gamma2.map_gamma2(m2)[0], Map_gamma.map_gamma(m2), 1e-6,
                "for order 2 gamma must equal the second eigenvalue");

        // testMap has a rank-one D1, hence is a renewal process: its ACF is zero
        // to roundoff, so acf(2)/acf(1) is meaningless and the guard must win
        assertEquals(0.0, Map_gamma.map_gamma(testMap()), TOL,
                "a renewal MAP has no autocorrelation decay");

        // order > 2: the ACF is NOT geometric, so gamma comes from a curve fit
        // and acf(2)/acf(1) is a ratio of near-zero numbers that blows up.
        // Erlang-3 is phase-type, hence uncorrelated: acf(2)/acf(1) is 0/0 here.
        MatrixCell erlang3 = Map_erlang.map_erlang(1.0, 3);
        assertEquals(0.0, Map_gamma.map_gamma(erlang3), TOL,
                "Erlang-3 is phase-type, so its ACF decay rate is 0");
    }

    @Test
    public void gammaMatchesMatlabOnHigherOrderMaps() {
        // Reference values from MATLAB map_gamma (kpctoolbox). The order-4 MAP is
        // the superposition of two MMPP(2)s. Before 2026-07 the JAR returned a
        // clamped acf(2)/acf(1) here, which is a different quantity entirely.
        MatrixCell superposed = Map_super.map_super(
                Map_mmpp2.map_mmpp2(1.0, 2.0, -1.0, 0.2),
                Map_mmpp2.map_mmpp2(1.0, 4.0, -1.0, 0.3));
        assertEquals(0.8144297260, Map_gamma.map_gamma(superposed), 1e-6,
                "map_gamma of an order-4 MAP must match the MATLAB reference");
    }

    @Test
    public void mmpp2RejectsInfeasibleRequestsInsteadOfReturningANonMap() {
        // The closed form solves the moment-matching equations without
        // constraining its solution to be a MAP. Outside the feasible set it
        // used to return negative rates silently; D0+D1 still had zero rowsums,
        // so the usual generator check did not catch it.

        // ACF1 above the maximum RHO0 = (1-1/SCV)/2 = 1/3 attainable at SCV=3
        assertThrows(IllegalArgumentException.class,
                () -> Map_mmpp2.map_mmpp2(1.0, 3.0, 0.05, 0.4),
                "ACF1 beyond the feasible maximum must be rejected");
        // an MMPP(2) is over-dispersed
        assertThrows(IllegalArgumentException.class,
                () -> Map_mmpp2.map_mmpp2(1.0, 0.5, -1.0, 0.0),
                "SCV<1 must be rejected");
        // SCV=1 divides by 1-1/SCV=0 and returns NaN rates
        assertThrows(IllegalArgumentException.class,
                () -> Map_mmpp2.map_mmpp2(1.0, 1.0, -1.0, 0.0),
                "the SCV=1 Poisson boundary must be rejected");
        // an MMPP(2) cannot be negatively autocorrelated; this used to fall into
        // the G2<1e-6 branch and silently return an uncorrelated MAP
        assertThrows(IllegalArgumentException.class,
                () -> Map_mmpp2.map_mmpp2(1.0, 3.0, -1.0, -0.2),
                "negative ACF1 must be rejected, not silently ignored");
        // third moment below the infimum 1.5*E2^2/E1 of the class
        assertThrows(IllegalArgumentException.class,
                () -> Map_mmpp2.map_mmpp2(1.0, 3.0, 2.0, 0.1),
                "a skewness below the feasible minimum must be rejected");

        // the sentinel ACF1=-1 asks for the maximum feasible autocorrelation and
        // must still be honoured: it yields RHO0 exactly, via a near-uncoupled
        // chain whose smallest rate is ~1e-9
        MatrixCell maxAcf = Map_mmpp2.map_mmpp2(1.0, 25.0, -1.0, -1.0);
        assertTrue(Map_isfeasible.map_isfeasible(maxAcf),
                "the ACF1=-1 sentinel must produce a feasible MAP");
        assertEquals(0.5 * (1.0 - 1.0 / 25.0),
                Map_acf.map_acf(maxAcf.get(0), maxAcf.get(1), 1).value(), 1e-6,
                "ACF1=-1 must attain the maximum feasible lag-1 autocorrelation");
    }

    @Test
    public void countMeanGrowsAsRateTimesTime() {
        // Aggregate count mean E[N(t)] = lambda * t for a marked MAP
        MatrixCell mmap = new MatrixCell(new Matrix[]{
                testMap().get(0), testMap().get(1), testMap().get(1)});
        double lambda = Map_lambda.map_lambda(testMap());
        double t = 10.0;
        assertEquals(lambda * t,
                Mmap_count_mean.mmap_count_mean(mmap, t).get(0), 1e-6,
                "aggregate count mean must be lambda*t");
    }
}
