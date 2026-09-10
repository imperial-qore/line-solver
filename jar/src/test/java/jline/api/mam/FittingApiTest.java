/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.mam;

import jline.io.Ret;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Fit-consistency validation of the moment/autocorrelation fitting APIs
 * (jline.api.mam): a fitted process must reproduce the moments (and, where
 * fitted, the lag-1 autocorrelation decay rate gamma) it was fitted to.
 * Target moments are generated from known analytic processes.
 */
public class FittingApiTest {

    private static final double FIT_RTOL = 1e-6;

    // Moments of a HyperExp(p=0.4, l1=4, l2=1): mk = k! (p/l1^k + (1-p)/l2^k)
    private static double heMoment(int k) {
        double fact = 1;
        for (int f = 2; f <= k; f++) {
            fact *= f;
        }
        return fact * (0.4 / Math.pow(4.0, k) + 0.6 / Math.pow(1.0, k));
    }


    @Test
    public void aphFitReproducesHyperExpMoments() {
        double m1 = heMoment(1), m2 = heMoment(2), m3 = heMoment(3);
        MatrixCell aph = Aph_fit.aph_fit(m1, m2, m3);
        assertNotNull(aph, "aph_fit returned null");
        assertEquals(m1, Map_moment.map_moment(aph, 1), 1e-4 * m1, "fitted first moment");
        assertEquals(m2, Map_moment.map_moment(aph, 2), 1e-3 * m2, "fitted second moment");
        assertEquals(m3, Map_moment.map_moment(aph, 3), 1e-2 * m3, "fitted third moment");
    }

    @Test
    public void aph2FitReproducesMoments() {
        double m1 = heMoment(1), m2 = heMoment(2), m3 = heMoment(3);
        Ret.mamAPH2Fit fit = Aph2_fit.aph2_fit(m1, m2, m3);
        assertNotNull(fit, "aph2_fit returned null");
        assertNotNull(fit.APH, "aph2_fit returned null APH");
        assertFalse(fit.APH.isEmpty(), "aph2_fit returned empty APH");
        assertEquals(m1, Map_moment.map_moment(fit.APH, 1), 1e-4 * m1,
                "APH(2) fitted first moment");
        assertEquals(m2, Map_moment.map_moment(fit.APH, 2), 1e-3 * m2,
                "APH(2) fitted second moment");
    }

    @Test
    public void amap2FitGammaReproducesMomentsAndGamma() {
        // Target: HyperExp renewal moments with a mild positive autocorrelation
        double m1 = heMoment(1), m2 = heMoment(2), m3 = heMoment(3);
        double gamma = 0.2;
        Pair<MatrixCell, List<MatrixCell>> fit =
                Amap2_fit_gamma.amap2_fit_gamma(m1, m2, m3, gamma);
        assertNotNull(fit, "amap2_fit_gamma returned null");
        MatrixCell map = fit.getLeft();
        assertNotNull(map, "amap2_fit_gamma returned null MAP");
        assertTrue(Map_isfeasible.map_isfeasible(map), "fitted AMAP(2) must be feasible");
        assertEquals(m1, Map_mean.map_mean(map), 1e-3 * m1, "AMAP(2) fitted mean");
        assertEquals(gamma, Map_gamma.map_gamma(map), 0.05,
                "AMAP(2) fitted autocorrelation decay rate");
    }

    @Test
    public void aphRandProducesFeasibleDistributions() {
        for (int k = 1; k <= 3; k++) {
            MatrixCell aph = Aph_rand.aph_rand(k);
            assertNotNull(aph, "aph_rand(" + k + ") returned null");
            assertTrue(Map_isfeasible.map_isfeasible(aph),
                    "aph_rand(" + k + ") must produce a feasible process");
            double mean = Map_mean.map_mean(aph);
            assertTrue(mean > 0 && Double.isFinite(mean),
                    "aph_rand(" + k + ") mean must be positive and finite");
        }
    }

    @Test
    public void mapRenewalProcessHasZeroGamma() {
        // A PH-renewal process (no correlation) must have gamma ~ 0
        MatrixCell exp = Map_exponential.map_exponential(0.5);
        assertEquals(0.0, Map_gamma.map_gamma(exp), 1e-8,
                "renewal process autocorrelation decay must be zero");
    }

    // ------------------------------------------------------------------
    // Measure-refit-remeasure self-consistency loops
    // ------------------------------------------------------------------

    /** A feasible 2-phase MMPP used as the measurement source. */
    private static MatrixCell sourceMmpp2() {
        // MMPP2 with rates (4, 1) and switching rates (0.5, 0.2)
        jline.util.matrix.Matrix d0 = new jline.util.matrix.Matrix(2, 2);
        d0.set(0, 0, -4.5); d0.set(0, 1, 0.5);
        d0.set(1, 0, 0.2);  d0.set(1, 1, -1.2);
        jline.util.matrix.Matrix d1 = new jline.util.matrix.Matrix(2, 2);
        d1.set(0, 0, 4.0);
        d1.set(1, 1, 1.0);
        return new MatrixCell(d0, d1);
    }

    @Test
    public void mapMmpp2RefitReproducesMeasuredStatistics() {
        MatrixCell src = sourceMmpp2();
        double mean = Map_mean.map_mean(src);
        double scv = Map_scv.map_scv(src);
        double skew = Map_skew.map_skew(src);
        double acf1 = Map_acf.map_acf(src.get(0), src.get(1), 1).get(0);

        MatrixCell fit = Map_mmpp2.map_mmpp2(mean, scv, skew, acf1);
        assertNotNull(fit, "map_mmpp2 returned null");
        assertTrue(Map_isfeasible.map_isfeasible(fit), "fitted MMPP2 must be feasible");
        assertEquals(mean, Map_mean.map_mean(fit), 1e-6 * mean, "MMPP2 refit mean");
        assertEquals(scv, Map_scv.map_scv(fit), 1e-5 * scv, "MMPP2 refit SCV");
        assertEquals(acf1, Map_acf.map_acf(fit.get(0), fit.get(1), 1).get(0), 1e-5,
                "MMPP2 refit lag-1 autocorrelation");
    }

    /**
     * map_mmpp2 must construct a feasible MMPP(2) hitting the requested mean,
     * SCV and lag-1 autocorrelation exactly. Reference values verified against
     * MATLAB map_mmpp2 for each row (both the G2>=1e-6 general branch and the
     * G2->0 boundary branch at acf1=0).
     */
    @ParameterizedTest(name = "mean={0} scv={1} acf1={3}")
    @CsvSource({
            "1.0, 2.0, -1.0, 0.1",
            "1.0, 3.0, -1.0, 0.2",
            "2.0, 2.0, -1.0, 0.05",
            "1.0, 4.0, -1.0, 0.15",
            "1.0, 2.0, -1.0, 0.0"
    })
    public void mapMmpp2HitsTargetStatistics(double mean, double scv, double skew, double acf1) {
        MatrixCell m = Map_mmpp2.map_mmpp2(mean, scv, skew, acf1);
        assertNotNull(m, "map_mmpp2 returned null");
        assertTrue(Map_isfeasible.map_isfeasible(m), "map_mmpp2 must produce a feasible MAP");
        assertEquals(mean, Map_mean.map_mean(m), 1e-6 * mean, "map_mmpp2 target mean");
        assertEquals(scv, Map_scv.map_scv(m), 1e-5 * scv, "map_mmpp2 target SCV");
        assertEquals(acf1, Map_acf.map_acf(m.get(0), m.get(1), 1).get(0), 1e-5,
                "map_mmpp2 target lag-1 autocorrelation");
    }

    @Test
    public void mapDistanceIsNonnegativeAndZeroToSelf() {
        // Squared L2 distance: zero to itself, strictly positive between
        // distinct MAPs. Cross value verified 3-way (MATLAB, native-Python).
        jline.util.matrix.Matrix a0 = new jline.util.matrix.Matrix(2, 2);
        a0.set(0, 0, -4.5); a0.set(0, 1, 0.5); a0.set(1, 0, 0.2); a0.set(1, 1, -1.2);
        jline.util.matrix.Matrix a1 = new jline.util.matrix.Matrix(2, 2);
        a1.set(0, 0, 4.0); a1.set(1, 1, 1.0);
        jline.util.matrix.Matrix b0 = new jline.util.matrix.Matrix(2, 2);
        b0.set(0, 0, -3.0); b0.set(0, 1, 0.3); b0.set(1, 0, 0.4); b0.set(1, 1, -2.0);
        jline.util.matrix.Matrix b1 = new jline.util.matrix.Matrix(2, 2);
        b1.set(0, 0, 2.5); b1.set(0, 1, 0.2); b1.set(1, 0, 0.3); b1.set(1, 1, 1.3);

        assertEquals(0.0, Map_dist.map_dist(a0, a1, a0, a1, 3), 1e-9, "map_dist self must be 0");
        double cross = Map_dist.map_dist(a0, a1, b0, b1, 3);
        assertTrue(cross > 0, "map_dist between distinct MAPs must be positive");
        assertEquals(1.1890107514295, cross, 1e-9, "map_dist cross vs MATLAB reference");

        assertEquals(0.0, Map_dist_acf.map_dist_acf(a0, a1, a0, a1), 1e-9,
                "map_dist_acf self must be 0");
        double crossAcf = Map_dist_acf.map_dist_acf(a0, a1, b0, b1);
        assertTrue(crossAcf > 0, "map_dist_acf between distinct MAPs must be positive");
        assertEquals(0.0430873486890663, crossAcf, 1e-9,
                "map_dist_acf cross vs MATLAB reference");
    }

    @Test
    public void map2FitReproducesMomentsAndGamma() {
        MatrixCell src = sourceMmpp2();
        double e1 = Map_moment.map_moment(src, 1);
        double e2 = Map_moment.map_moment(src, 2);
        double e3 = Map_moment.map_moment(src, 3);
        double g2 = Map_gamma.map_gamma(src);

        Ret.mamMAPFitReturn fit = Map2_fit.map2_fit(e1, e2, e3, g2);
        assertNotNull(fit.MAP, "map2_fit returned null MAP");
        assertFalse(fit.MAP.isEmpty(), "map2_fit returned empty MAP");
        assertEquals(e1, Map_moment.map_moment(fit.MAP, 1), 1e-6 * e1, "MAP(2) refit m1");
        assertEquals(e2, Map_moment.map_moment(fit.MAP, 2), 1e-5 * e2, "MAP(2) refit m2");
        assertEquals(g2, Map_gamma.map_gamma(fit.MAP), 1e-5, "MAP(2) refit gamma");
    }

    @Test
    public void mapMaxMatchesMatlabReference() {
        // For renewal inputs map_max is the max of independent random variables,
        // so E[max(X,Y)] = 1/l1 + 1/l2 - 1/(l1+l2) for exp(l1), exp(l2)
        MatrixCell a = Map_exponential.map_exponential(1.0 / 2.0);
        MatrixCell b = Map_exponential.map_exponential(1.0 / 3.0);
        MatrixCell mx = Map_max.map_max(a, b);
        assertEquals(1.0 / 2.0 + 1.0 / 3.0 - 1.0 / 5.0, Map_mean.map_mean(mx), 1e-8,
                "map_max mean must match the MATLAB reference");
    }

    @Test
    public void mapMaxIsAValidGeneratorMatchingTheMaxDistribution() {
        // A mean alone cannot detect a malformed generator, so pin the whole
        // distribution: for exp(l1), exp(l2) the max has closed-form moments
        // E[T^k] = k! * (1/l1^k + 1/l2^k - 1/(l1+l2)^k)
        double l1 = 2.0;
        double l2 = 3.0;
        MatrixCell mx = Map_max.map_max(Map_exponential.map_exponential(1.0 / l1),
                Map_exponential.map_exponential(1.0 / l2));

        double factorial = 1.0;
        for (int k = 1; k <= 4; k++) {
            factorial *= k;
            double expected = factorial * (Math.pow(1.0 / l1, k) + Math.pow(1.0 / l2, k)
                    - Math.pow(1.0 / (l1 + l2), k));
            assertEquals(expected, Map_moment.map_moment(mx, k), 1e-8 * expected,
                    "map_max moment " + k);
        }

        // D0+D1 must be an infinitesimal generator; the pre-2026-07 construction
        // used a rateless D1 = pie'*ones and failed this
        Matrix gen = mx.get(0).add(1.0, mx.get(1));
        for (int i = 0; i < gen.getNumRows(); i++) {
            double rowsum = 0.0;
            for (int j = 0; j < gen.getNumCols(); j++) rowsum += gen.get(i, j);
            assertEquals(0.0, rowsum, 1e-12, "row " + i + " of D0+D1 must sum to zero");
        }

        // D1 = d*pie is rank one, so the output is a renewal process
        assertEquals(0.0, Map_gamma.map_gamma(mx), 1e-8, "map_max output must be renewal");

        // orders must be allowed to differ: this threw before the fix
        MatrixCell mixed = Map_max.map_max(Map_erlang.map_erlang(1.5, 2),
                Map_exponential.map_exponential(0.7));
        assertEquals(1.6631391201, Map_mean.map_mean(mixed), 1e-6,
                "map_max must accept MAPs of different order");
    }

    @Test
    public void mmapSuperpositionAddsAggregateRates() {
        // Superposition of independent streams is additive in the aggregate
        // rate (mmap_super, distinct from the mmap_sum convolution below)
        MatrixCell m1 = sourceMarkedMmap();
        MatrixCell m2 = sourceMarkedMmap();
        double r1 = Map_lambda.map_lambda(new MatrixCell(m1.get(0), m1.get(1)));
        MatrixCell sup = Mmap_super.mmap_super(m1, m2, "match");
        double rs = Map_lambda.map_lambda(new MatrixCell(sup.get(0), sup.get(1)));
        assertEquals(2 * r1, rs, 1e-6 * 2 * r1,
                "superposed marked MAP must have additive aggregate rate");
    }

    @Test
    public void mmapSumConvolvesInterarrivalTimes() {
        // mmap_sum(mmap, n): interarrival = sum of n base interarrivals, so
        // the mean scales by n and the aggregate rate divides by n. The class
        // marking is inherited from the last of the n base arrivals, so the
        // per-class split is unchanged.
        MatrixCell base = sourceMarkedMmap();
        double baseMean = Map_mean.map_mean(new MatrixCell(base.get(0), base.get(1)));
        jline.util.matrix.Matrix baseRates = Mmap_lambda.mmap_lambda(base);

        int n = 3;
        MatrixCell summed = Mmap_sum.mmap_sum(base, n);
        assertEquals(base.size(), summed.size(), "class count must be preserved");
        assertEquals(n * base.get(0).getNumRows(), summed.get(0).getNumRows(),
                "order must scale by n");
        assertTrue(Map_isfeasible.map_isfeasible(new MatrixCell(summed.get(0), summed.get(1))),
                "summed MMAP must be a feasible MAP");
        assertEquals(n * baseMean,
                Map_mean.map_mean(new MatrixCell(summed.get(0), summed.get(1))), 1e-6,
                "summed interarrival mean must be n times the base mean");

        // Per-class rate = (per-class split) * (aggregate rate / n)
        jline.util.matrix.Matrix summedRates = Mmap_lambda.mmap_lambda(summed);
        for (int k = 0; k < baseRates.getNumElements(); k++) {
            assertEquals(baseRates.get(k) / n, summedRates.get(k), 1e-6 * baseRates.get(k),
                    "summed class " + k + " rate must be the base rate divided by n");
        }
    }

    @Test
    public void dmapMomentMatchesBruteForceRawMoments() {
        // Discrete MAP; raw interarrival moments verified by brute-force sum
        // over P(X=k)=al*D0^(k-1)*D1*e (matches MATLAB dmap_moment after fix)
        jline.util.matrix.Matrix d0 = new jline.util.matrix.Matrix(2, 2);
        d0.set(0, 0, 0.2); d0.set(0, 1, 0.1); d0.set(1, 0, 0.05); d0.set(1, 1, 0.3);
        jline.util.matrix.Matrix d1 = new jline.util.matrix.Matrix(2, 2);
        d1.set(0, 0, 0.4); d1.set(0, 1, 0.3); d1.set(1, 0, 0.25); d1.set(1, 1, 0.4);
        assertEquals(1.489362, Dmap_moment.dmap_moment(d0, d1, 1), 1e-5, "DMAP raw m1");
        assertEquals(2.957638, Dmap_moment.dmap_moment(d0, d1, 2), 1e-5, "DMAP raw m2");
        assertEquals(8.103912, Dmap_moment.dmap_moment(d0, d1, 3), 1e-5, "DMAP raw m3");

        // dmap_pie must be the stationary vector of the embedded chain
        // P=(I-D0)^-1 D1: a probability vector satisfying pi P = pi
        jline.util.matrix.Matrix pie = Dmap_pie.dmap_pie(d0, d1);
        double sum = 0.0;
        for (int i = 0; i < pie.getNumElements(); i++) {
            assertTrue(pie.get(i) >= -1e-12, "dmap_pie entries must be nonnegative");
            sum += pie.get(i);
        }
        assertEquals(1.0, sum, 1e-9, "dmap_pie must be a probability vector");
    }

    /** Deterministic 2-class marked MMAP built on sourceMmpp2. */
    private static MatrixCell sourceMarkedMmap() {
        MatrixCell base = sourceMmpp2();
        return new MatrixCell(new jline.util.matrix.Matrix[]{
                base.get(0), base.get(1),
                base.get(1).scale(0.6), base.get(1).scale(0.4)});
    }

    @Test
    public void mamap2mFitMmapPreservesClassRates() {
        // MATLAB reference on the same input: class rates (1.114286, 0.742857)
        MatrixCell src = sourceMarkedMmap();
        jline.util.matrix.Matrix srcRates = Mmap_lambda.mmap_lambda(src);
        assertEquals(1.114286, srcRates.get(0), 1e-5, "source class 1 rate sanity");
        MatrixCell fit = Mamap2m_fit_mmap.mamap2m_fit_mmap(src);
        assertNotNull(fit, "mamap2m_fit_mmap returned null");
        jline.util.matrix.Matrix fitRates = Mmap_lambda.mmap_lambda(fit);
        for (int k = 0; k < srcRates.getNumElements(); k++) {
            assertEquals(srcRates.get(k), fitRates.get(k), 1e-4 * srcRates.get(k),
                    "MAMAP(2,m) fit must preserve class " + k + " rate");
        }
    }

    @Test
    public void aph2AdjustRestoresFeasibility() {
        // Moments infeasible for APH(2): scv < 0.5 needs more phases; the
        // adjuster must move (M2, M3) to the nearest feasible set with M1 kept
        double m1 = 1.0;
        double m2 = 1.2;  // scv = 0.2, infeasible for APH(2)
        double m3 = 1.6;
        java.util.Map<Integer, Double> adjusted = Aph2_adjust.aph2_adjust(m1, m2, m3);
        assertNotNull(adjusted, "aph2_adjust returned null");
        double m2a = adjusted.get(0);
        double m3a = adjusted.get(1);
        Ret.mamAPH2Fit fit = Aph2_fit.aph2_fit(m1, m2a, m3a);
        assertFalse(fit.APH.isEmpty(), "adjusted moments must be APH(2)-feasible");
        assertEquals(m1, Map_moment.map_moment(fit.APH, 1), 1e-6,
                "aph2_adjust must preserve the first moment");
    }
}
