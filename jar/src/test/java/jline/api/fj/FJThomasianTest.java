/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.fj;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.function.DoubleUnaryOperator;

import org.junit.jupiter.api.Test;

/**
 * Known-answer tests for the formulas ported from A. Thomasian, "Analysis of
 * Fork/Join and Related Queueing Systems", ACM Computing Surveys 47(2),
 * Article 17, 2014.
 *
 * Every expected value is either the survey's own worked number, an exact
 * identity the formula must reproduce, or the value MATLAB's fj_* twin returns
 * (matlab/src/api/fj/), which is the reference implementation.
 */
public class FJThomasianTest {

    private static final double TOL = 1e-9;

    @Test
    public void qgbReducesToTheOrdinaryGeometricBound() {
        double[] D = {1.0, 2.0, 3.0};
        int[] P = {1, 1, 1};
        FJ_closed.FJQgbResult g = FJ_closed.fj_qgb(D, P, 7, 1.5);
        double denom = 1.5 + 6.0 + 3.0 * 7.0;
        for (int i = 0; i < 3; i++) {
            double y = D[i] * 7.0 / denom;
            assertEquals(y / (1 - y) - Math.pow(y, 8) / (1 - y), g.Q[i], 1e-12);
        }
    }

    @Test
    public void amvaWithUnitForkDegreesIsExactMva() {
        double[] D = {1.0, 2.0, 3.0};
        int[] P = {1, 1, 1};
        FJ_closed.FJAmvaResult a = FJ_closed.fj_amva(D, P, 7, 1.5);
        double[] Q = new double[3];
        double X = 0;
        for (int m = 1; m <= 7; m++) {
            double[] R = new double[3];
            double tot = 0;
            for (int i = 0; i < 3; i++) {
                R[i] = D[i] * (1 + Q[i]);
                tot += R[i];
            }
            X = m / (1.5 + tot);
            for (int i = 0; i < 3; i++) {
                Q[i] = X * R[i];
            }
        }
        assertEquals(X, a.X, 1e-12);
        for (int i = 0; i < 3; i++) {
            assertEquals(Q[i], a.Q[i], 1e-12);
        }
    }

    @Test
    public void resptClosedIsTightAtTwoBranches() {
        FJ_closed.FJResptClosedResult r = FJ_closed.fj_respt_closed(2, 0.4, 5);
        assertEquals(0.4 * (1.5 + 4), r.R, TOL);
        assertTrue(r.exact);
    }

    @Test
    public void xmaxHetMatchesTheTextbookTwoVariableAnswer() {
        assertEquals(1 + 0.5 - 1.0 / 3.0, FJ_maxima.fj_xmax_het(new double[] {1.0, 2.0}), TOL);
        assertEquals(FJ_harmonic.fj_harmonic(4) / 2.0,
                FJ_maxima.fj_xmax_het(new double[] {2.0, 2.0, 2.0, 2.0}), TOL);
    }

    @Test
    public void momentRecurrenceAgreesWithInclusionExclusion() {
        double[] lam = {1.0, 2.0, 3.0, 5.0};
        double[] m = FJ_maxima.fj_xmax_moments_het(lam, 3);
        for (int n = 1; n <= 3; n++) {
            assertEquals(FJ_maxima.fj_xmax_het(lam, n), m[n - 1], 1e-10);
        }
    }

    @Test
    public void transformIsOneAtTheOriginAndItsSlopeIsTheMean() {
        double[] lam = {1.0, 2.0, 3.0, 5.0};
        assertEquals(1.0, FJ_maxima.fj_lst_max_het(lam, 0.0), 1e-12);
        double h = 1e-5;
        assertEquals(FJ_maxima.fj_xmax_het(lam),
                -(FJ_maxima.fj_lst_max_het(lam, h) - 1.0) / h, 1e-4);
    }

    @Test
    public void harrisonZertalIsExactForTheExponential() {
        double m1 = 0.7;
        assertEquals(FJ_harmonic.fj_harmonic(6) * m1, FJ_maxima.fj_xmax_hz(m1, 2 * m1 * m1, 6), TOL);
    }

    @Test
    public void harrisonZertalRecurrenceIsExactForIidExponentials() {
        int K = 4;
        final double lam = 1.3;
        double[] m1 = new double[K];
        double[] m2 = new double[K];
        DoubleUnaryOperator[] cdf = new DoubleUnaryOperator[K];
        for (int i = 0; i < K; i++) {
            m1[i] = 1 / lam;
            m2[i] = 2 / (lam * lam);
            cdf[i] = new DoubleUnaryOperator() {
                @Override
                public double applyAsDouble(double t) {
                    return t > 0 ? 1 - Math.exp(-lam * t) : 0.0;
                }
            };
        }
        assertEquals(FJ_harmonic.fj_harmonic(K) / lam,
                FJ_maxima.fj_xmax_hz_het(m1, m2, cdf), 1e-6);
    }

    @Test
    public void characteristicMaximumBoundsTheLatticeMaximum() {
        FJ_charmax_ext.FJCharMaxDiscreteResult g =
                FJ_charmax_ext.fj_char_max_discrete(8, FJ_charmax_ext.DiscreteDist.GEOMETRIC, 0.6);
        assertTrue(g.MK >= g.exact - 1e-9);
        FJ_charmax_ext.FJCharMaxDiscreteResult p =
                FJ_charmax_ext.fj_char_max_discrete(8, FJ_charmax_ext.DiscreteDist.POISSON, 4.0);
        assertTrue(p.MK >= p.exact - 1e-9);
    }

    @Test
    public void blomPositionSitsInsideTheKruskalWeissBracket() {
        FJ_charmax_ext.FJCharMaxBlomResult b = FJ_charmax_ext.fj_char_max_blom(20);
        assertTrue(b.bracketAvailable);
        assertTrue(b.mK > b.lo && b.mK < b.hi);
    }

    @Test
    public void coxianFitReproducesItsTargetMomentsAndDegeneratesToTheExponential() {
        FJ_maxima.FJCoxFitResult f = FJ_maxima.fj_cox_fit(2.0, 1.5);
        FJ_maxima.FJXmaxCoxianResult c1 = FJ_maxima.fj_xmax_coxian(1, f.mu1, f.mu2, f.q);
        assertEquals(2.0, c1.m1, 1e-12);
        assertEquals(1.5, c1.c2, 1e-12);
        assertEquals(2.0, c1.Xmax, 1e-12);
        assertEquals(FJ_harmonic.fj_harmonic(5),
                FJ_maxima.fj_xmax_coxian(5, 1.0, 1.0, 0.0).Xmax, 1e-10);
    }

    @Test
    public void dispersionOfTwoExponentialBranches() {
        double mu = 1.7;
        FJ_dispersion.FJDispersionResult d =
                FJ_dispersion.fj_dispersion(new int[] {1, 1}, new double[] {mu, mu});
        assertEquals(1.5 / mu, d.Emax, 1e-7);
        assertEquals(0.5 / mu, d.Emin, 1e-7);
        assertEquals(1.0 / mu, d.Edisp, 1e-7);
    }

    @Test
    public void delayingNeverIncreasesTheDispersion() {
        int[] shape = {1, 3, 2};
        double[] rate = {1.0, 2.0, 0.8};
        double d0 = FJ_dispersion.fj_dispersion(shape, rate).Edisp;
        FJ_dispersion.FJDelayOptResult o = FJ_dispersion.fj_delay_opt(shape, rate);
        assertTrue(o.Edisp <= d0 + 1e-9);
    }

    @Test
    public void noSplittingCollapsesToMm1AtOneTask() {
        assertEquals(1.0 / 0.9, FJ_parallel.fj_respt_nosplit(1, 0.5, 1.4).R, TOL);
    }

    @Test
    public void bulkArrivalsCollapseToMm1AtUnitBatchAndOneServer() {
        FJ_parallel.FJResptBulkResult b = FJ_parallel.fj_respt_bulk(1, 0.5, 1.4, 1);
        double rho = 0.5 / 1.4;
        assertEquals(rho / (1 - rho), b.Q, 1e-7);
        assertEquals(1.0 / 0.9, b.Rreq, 1e-7);
    }

    @Test
    public void independentServerModelReproducesTheSurveysUtilization() {
        FJ_parallel.FJIsmGreenResult g =
                FJ_parallel.fj_ism_green(1.0, 1.4, 4, new double[] {0.1, 0.2, 0.3, 0.4});
        // The survey quotes rho = 0.9286 for this instance
        assertEquals(0.9285714285714286, g.rho, 1e-9);
        assertTrue(g.W > 0);
        assertTrue(g.pq > 0 && g.pq < 1);
        assertTrue(g.pd > 0 && g.pd <= 1);
    }

    @Test
    public void teamServiceCapacityMatchesTheSurveysWorkedExample() {
        double[] f = {0.25, 0.25, 0.25, 0.25};
        int[] r = {1, 2, 3, 4};
        double[] x = {1, 1, 1, 1};
        FJ_parallel.FJTsmCapacityResult t = FJ_parallel.fj_tsm_capacity(4, f, r, x);
        // Lambda_max = s / sum f r x = 4 / 2.5 = 1.6, attained because every state
        // carrying probability is full capacity
        assertEquals(1.6, t.Lmax, 1e-12);
        assertEquals(1.6, t.Llp, 1e-7);
        // The skewed frequency vector of the survey reaches only 4/3
        double eps = 1e-6;
        double[] fs = {eps, 0.5 - eps, 0.5 - eps, eps};
        FJ_parallel.FJTsmCapacityResult t2 = FJ_parallel.fj_tsm_capacity(4, fs, r, x);
        assertEquals(4.0 / 3.0, t2.Llp, 1e-4);
    }

    @Test
    public void serializationBlockingProbability() {
        FJ_parallel.FJSerializationResult s =
                FJ_parallel.fj_serialization(new double[] {0.2, 0.3}, 1.0, 5);
        assertEquals(1 - Math.pow(1 - 0.2 / 1.5, 4), s.P[0], 1e-12);
    }

    @Test
    public void dagMakespanReproducesBothWorkedExamples() {
        boolean[][] pred = new boolean[2][2];
        double[][] rate = {{1.0 / 10, 1.0 / 15}, {1.0 / 20, 1.0 / 30}};
        FJ_parallel.FJDagMakespanResult d = FJ_parallel.fj_dag_makespan(pred, rate);
        // The survey's coupled two-task example gives 26.77 by hand, 80/3 exactly
        assertEquals(80.0 / 3.0, d.C, 1e-9);
        assertEquals(40.0 / 3.0, d.Cend[0], 1e-9);
        assertEquals(70.0 / 3.0, d.Cend[1], 1e-9);
        double[][] rate2 = {{1.0, 1.0}, {0.1, 0.1}};
        FJ_parallel.FJDagMakespanResult d2 = FJ_parallel.fj_dag_makespan(pred, rate2);
        // The state-truncation example gives 10.1
        assertEquals(1 / 1.1 + (1 / 1.1) * 10 + (0.1 / 1.1) * 1, d2.C, 1e-9);
    }

    @Test
    public void dagMakespanOfAChainIsTheSumOfTheMeans() {
        boolean[][] pred = new boolean[3][3];
        pred[0][1] = true;
        pred[1][2] = true;
        double[][] rate = new double[3][3];
        for (int k = 0; k < 3; k++) {
            rate[0][k] = 1.0;
            rate[1][k] = 2.0;
            rate[2][k] = 4.0;
        }
        FJ_parallel.FJDagMakespanResult d = FJ_parallel.fj_dag_makespan(pred, rate);
        assertEquals(1 + 0.5 + 0.25, d.C, 1e-12);
        assertEquals(1.5, d.I[2], 1e-12);
    }
}
