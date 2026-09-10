/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api;

import jline.api.mam.Amap2_fit_gamma;
import jline.api.mam.Aph2_fit;
import jline.api.mam.Aph_fit;
import jline.api.mam.Map2_fit;
import jline.api.mam.Map_erlang;
import jline.api.mam.Map_gamma2;
import jline.api.mam.Map_hyperexp;
import jline.api.mam.Map_mmpp2;
import jline.api.mam.Map_moment;
import jline.api.mam.Map_scv;
import jline.api.mam.Maph2m_fit;
import jline.api.mam.Mmap_pc;
import jline.api.mam.Mmpp2_fit;
import jline.io.Ret;
import jline.util.Pair;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Drift guard for the closed-form MAP/MMPP/APH fitting formulas.
 *
 * Two kinds of assertion, both required:
 * - ROUND TRIP: the fitted process must reproduce the statistics it was asked
 *   to match. This catches a mis-transcribed expression regardless of which
 *   codebase drifts, and needs no golden value.
 * - PINNED VALUES: a handful of rates taken from the MATLAB reference (which
 *   was verified symbolically in SageMath for the MMPP(2) closed form). This
 *   catches a change that is self-consistent but no longer the same formula.
 *
 * See _kb/03-api-layer.md, "MMPP(2) moment fitting".
 */
public class MamFittingFormulaTest {

    private static final double TOL = 1e-9;

    private static void assertMoments(MatrixCell map, double e1, double e2, double e3, String what) {
        assertEquals(e1, Map_moment.map_moment(map, 1), Math.abs(e1) * TOL, what + ": E1");
        assertEquals(e2, Map_moment.map_moment(map, 2), Math.abs(e2) * TOL, what + ": E2");
        assertEquals(e3, Map_moment.map_moment(map, 3), Math.abs(e3) * TOL, what + ": E3");
    }

    @Test
    public void mmpp2FitMatchesMomentsAndDecayRate() {
        double[][] cases = {{1, 5, 45, 0.3}, {0.5, 1, 3.6, 0.8}, {2, 12, 216, 0.05}, {1, 21, 1000, 1e-4}};
        for (double[] c : cases) {
            MatrixCell map = Mmpp2_fit.mmpp2_fit(c[0], c[1], c[2], c[3]);
            String what = "mmpp2_fit" + java.util.Arrays.toString(c);
            assertMoments(map, c[0], c[1], c[2], what);
            assertEquals(c[3], Map_gamma2.map_gamma2(map)[0], 1e-8, what + ": gamma2");
            // an MMPP(2) has non-negative rates by construction
            assertTrue(map.get(1).get(0, 0) >= 0 && map.get(1).get(1, 1) >= 0
                    && map.get(0).get(0, 1) >= 0 && map.get(0).get(1, 0) >= 0, what + ": not a MAP");
        }
    }

    @Test
    public void mmpp2FitMatchesTheMatlabClosedForm() {
        // Pinned against MATLAB mmpp2_fit3.m at (E1,E2,E3,G2) = (1,5,45,0.3);
        // the same values are produced by native Python and the C++ header.
        MatrixCell map = Mmpp2_fit.mmpp2_fit(1, 5, 45, 0.3);
        assertEquals(0.11835708875252297, map.get(1).get(0, 0), 1e-12, "mu00");
        assertEquals(3.0416429112474770, map.get(1).get(1, 1), 1e-12, "mu11");
        assertEquals(0.25333822637151965, map.get(0).get(0, 1), 1e-12, "q01");
        assertEquals(0.58666177362848038, map.get(0).get(1, 0), 1e-12, "q10");
    }

    @Test
    public void mapMmpp2RealizesTheRequestedLag1Autocorrelation() {
        // rho1 = gamma2 * (1 - 1/SCV)/2 holds for every MMPP(2) (proved in Sage),
        // so map_mmpp2 must return exactly the ACF1 it is given.
        double mean = 1, scv = 4, skew = 3.5, acf1 = 0.2;
        MatrixCell map = Map_mmpp2.map_mmpp2(mean, scv, skew, acf1);
        assertEquals(mean, Map_moment.map_moment(map, 1), TOL, "mean");
        assertEquals(scv, Map_scv.map_scv(map), 1e-8, "scv");
        assertEquals(acf1 / ((1 - 1 / scv) / 2), Map_gamma2.map_gamma2(map)[0], 1e-8, "gamma2");
    }

    @Test
    public void map2FitMatchesMomentsAndDecayRate() {
        double[][] cases = {{1, 5, 45, 0.3}, {1, 3, 20, 0.1}, {0.5, 1, 3.6, 0.5}, {2, 12, 216, 0.05}};
        for (double[] c : cases) {
            Ret.mamMAPFitReturn r = Map2_fit.map2_fit(c[0], c[1], c[2], c[3]);
            String what = "map2_fit" + java.util.Arrays.toString(c);
            assertMoments(r.MAP, c[0], c[1], c[2], what);
            assertEquals(c[3], Map_gamma2.map_gamma2(r.MAP)[0], 1e-8, what + ": gamma2");
        }
    }

    @Test
    public void aphFitMatchesThreeMoments() {
        double[][] cases = {{1, 3, 15}, {1, 2, 6.5}, {1, 5, 60}, {2, 9, 60}, {1, 1.2, 1.7}};
        for (double[] c : cases) {
            MatrixCell map = Aph_fit.aph_fit(c[0], c[1], c[2]);
            assertMoments(map, c[0], c[1], c[2], "aph_fit" + java.util.Arrays.toString(c));
        }
    }

    @Test
    public void aph2FitMatchesFeasibleMoments() {
        double[][] cases = {{1, 3, 15}, {1, 2.5, 10}, {0.5, 0.75, 1.8}};
        for (double[] c : cases) {
            Ret.mamAPH2Fit r = Aph2_fit.aph2_fit(c[0], c[1], c[2]);
            assertMoments(r.APH, c[0], c[1], c[2], "aph2_fit" + java.util.Arrays.toString(c));
        }
    }

    @Test
    public void amap2FitGammaMatchesMomentsAndDecayRate() {
        double[][] cases = {{1, 3, 15, 0.3}, {1, 5, 45, 0.5}, {1, 2.5, 10, -0.2}};
        for (double[] c : cases) {
            Pair<MatrixCell, List<MatrixCell>> r =
                    Amap2_fit_gamma.amap2_fit_gamma(c[0], c[1], c[2], c[3]);
            String what = "amap2_fit_gamma" + java.util.Arrays.toString(c);
            assertMoments(r.getLeft(), c[0], c[1], c[2], what);
            assertEquals(c[3], Map_gamma2.map_gamma2(r.getLeft())[0], 1e-8, what + ": gamma2");
        }
    }

    @Test
    public void mapErlangAndHyperexpMatchTheirMoments() {
        for (int k : new int[]{2, 3, 5}) {
            MatrixCell map = Map_erlang.map_erlang(1.5, k);
            assertEquals(1.5, Map_moment.map_moment(map, 1), TOL, "map_erlang mean");
            assertEquals(1.0 / k, Map_scv.map_scv(map), 1e-8, "map_erlang scv");
        }
        // p = 0.5 with SCV = 4 is only reachable after the fallback to a smaller
        // p, and only through the second root of the quadratic: a regression here
        // means the feasible root is being discarded again
        MatrixCell hyp = Map_hyperexp.map_hyperexp(1.0, 4.0, 0.5);
        assertNotNull(hyp, "map_hyperexp(1,4,0.5) must return a feasible fit");
        assertEquals(1.0, Map_moment.map_moment(hyp, 1), 1e-8, "map_hyperexp mean");
        assertEquals(4.0, Map_scv.map_scv(hyp), 1e-8, "map_hyperexp scv");
    }

    @Test
    public void maph2mFitMatchesClassProbabilitiesExactly() {
        // maph2m_fit fits the class probabilities exactly and the backward
        // moments as closely as the APH(2) form allows. Backward moments must be
        // consistent with the mean: sum_c p_c b_c = M1.
        Matrix P = new Matrix(2, 1);
        P.set(0, 0, 0.6);
        P.set(1, 0, 0.4);
        Matrix B = new Matrix(2, 1);
        B.set(0, 0, 1.0769230769230769);
        B.set(1, 0, 0.88461538461538458);
        MatrixCell maph = Maph2m_fit.maph2m_fit(1, 3, 15, P, B);
        Matrix pc = Mmap_pc.mmap_pc(maph);
        assertEquals(0.6, pc.get(0), 1e-6, "class 1 probability");
        assertEquals(0.4, pc.get(1), 1e-6, "class 2 probability");
        // the class matrices must sum to D1
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                assertEquals(maph.get(1).get(i, j), maph.get(2).get(i, j) + maph.get(3).get(i, j),
                        1e-6, "class matrices must sum to D1 at (" + i + "," + j + ")");
            }
        }
        // pinned against MATLAB maph2m_fit.m for the same input
        assertEquals(1.107692307692308, maph.get(2).get(0, 0), 1e-6, "D11(0,0)");
        assertEquals(0.89230769230769258, maph.get(3).get(0, 0), 1e-6, "D12(0,0)");
    }
}
