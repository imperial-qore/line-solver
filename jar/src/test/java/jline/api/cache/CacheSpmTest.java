/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.cache;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.Timeout;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Boundaries of the SPM cache normalizing constant.
 *
 * <p>cache_spm returns Ehat(m) = E(m) * prod_l m_l!, the same constant
 * cache_erec evaluates exactly, so cache_erec is the oracle throughout.
 * Two boundaries used to be mishandled in every codebase:</p>
 *
 * <ul>
 * <li>m_l == 0: list l has xi(l)=0, a boundary of the Laplace integral rather
 * than a direction of it. Left in, the -(1/2) sum_l log xi(l) prefactor gained
 * ~+17 per empty list -- n=12 at gamma=(.8,.6,.4) and m=(0,3,3) returned
 * lZ=24.814 against an exact 9.127.</li>
 * <li>n == sum(m): every item is cached, the multipliers diverge and the
 * fixed-point iteration cannot converge. It used to be run anyway, and did not
 * terminate; the @Timeout on those cases is the regression guard.</li>
 * </ul>
 */
public class CacheSpmTest {

    /** Uniform popularity: 12 items, three lists. */
    private static Matrix uniform() {
        Matrix g = new Matrix(12, 3);
        for (int i = 0; i < 12; i++) {
            g.set(i, 0, 0.8);
            g.set(i, 1, 0.6);
            g.set(i, 2, 0.4);
        }
        return g;
    }

    /** Non-uniform popularity, monotone across the list index as SPM assumes. */
    private static Matrix skewed() {
        Matrix g = new Matrix(12, 3);
        for (int i = 0; i < 12; i++) {
            double gi = 0.3 + 2.7 * ((i + 1.0) / 12.0);
            g.set(i, 0, gi);
            g.set(i, 1, 0.6 * gi);
            g.set(i, 2, 0.3 * gi);
        }
        return g;
    }

    private static Matrix row(double... v) {
        Matrix out = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) { out.set(0, i, v[i]); }
        return out;
    }

    /** gamma restricted to the given columns, in order. */
    private static Matrix cols(Matrix g, int... keep) {
        Matrix out = new Matrix(g.getNumRows(), keep.length);
        for (int i = 0; i < g.getNumRows(); i++) {
            for (int a = 0; a < keep.length; a++) { out.set(i, a, g.get(i, keep[a])); }
        }
        return out;
    }

    private static double exact(Matrix g, Matrix m) {
        return Math.log(Cache_erec.cache_erec(g, m).get(0));
    }

    // ----------------------------------------------------------------------
    // m_l == 0: a zero-capacity list must leave the expansion
    // ----------------------------------------------------------------------

    @Test
    public void zeroCapacityListEqualsDroppingIt() {
        Matrix[] profiles = {uniform(), skewed()};
        double[][] caps = {{0, 3, 3}, {3, 0, 3}, {0, 0, 6}, {4, 4, 0}};
        int[][] keeps = {{1, 2}, {0, 2}, {2}, {0, 1}};
        for (Matrix g : profiles) {
            for (int c = 0; c < caps.length; c++) {
                Matrix m = row(caps[c]);
                Matrix mk = new Matrix(1, keeps[c].length);
                for (int a = 0; a < keeps[c].length; a++) { mk.set(0, a, caps[c][keeps[c][a]]); }

                Ret.cacheSpm withEmpty = Cache_spm.cache_spm(g, m);
                Ret.cacheSpm without = Cache_spm.cache_spm(cols(g, keeps[c]), mk);
                // dropping an empty list is exact, so the two agree bit for bit
                assertEquals(without.lZ, withEmpty.lZ, 0.0);
                // and SPM's own error, not a 1e8 blow-up, separates it from exact
                assertEquals(exact(g, m), withEmpty.lZ, 0.30);
                // the dropped lists carry xi = 0, the root of their capacity equation
                for (int l = 0; l < 3; l++) {
                    assertEquals(caps[c][l] == 0.0, withEmpty.xi.get(l) == 0.0);
                }

                // the deprecated rayint alias carries the same fix
                assertEquals(Cache_rayint.cache_rayint(cols(g, keeps[c]), mk).lZ,
                        Cache_rayint.cache_rayint(g, m).lZ, 0.0);
            }
        }
    }

    @Test
    public void emptyCacheIsTheUnitConstant() {
        Ret.cacheSpm r = Cache_spm.cache_spm(uniform(), row(0, 0, 0));
        assertEquals(1.0, r.Z, 0.0);
        assertEquals(0.0, r.lZ, 0.0);
        for (int l = 0; l < 3; l++) { assertEquals(0.0, r.xi.get(l), 0.0); }
    }

    // ----------------------------------------------------------------------
    // n == sum(m): degenerate saddle, exact fallback, no iteration
    // ----------------------------------------------------------------------

    @Test
    @Timeout(60)
    public void fullCacheReturnsTheExactConstantAndTerminates() {
        Matrix[] profiles = {uniform(), skewed()};
        double[][] caps = {{4, 4, 4}, {12, 0, 0}, {6, 5, 1}};
        for (Matrix g : profiles) {
            for (double[] cap : caps) {
                Matrix m = row(cap);
                Ret.cacheSpm r = Cache_spm.cache_spm(g, m);
                assertEquals(exact(g, m), r.lZ, 1e-9);
                assertEquals(Cache_erec.cache_erec(g, m).get(0), r.Z, 1e-6 * Math.abs(r.Z));
                for (int l = 0; l < 3; l++) {
                    // the multipliers' limit, not a number
                    assertTrue(Double.isInfinite(r.xi.get(l)));
                }
                assertEquals(exact(g, m), Cache_rayint.cache_rayint(g, m).lZ, 1e-9);
            }
        }
    }

    // ----------------------------------------------------------------------
    // interior points must be untouched by the boundary handling
    // ----------------------------------------------------------------------

    @Test
    public void interiorSaddleIsUnchanged() {
        Matrix g = uniform();
        double[][] cases = {{3, 3, 3, 13.34844416820355}, {4, 3, 2, 14.04836686545269},
                            {2, 2, 2, 10.23839884632498}, {6, 4, 1, 15.878606855246987}};
        for (double[] c : cases) {
            Matrix m = row(c[0], c[1], c[2]);
            assertEquals(c[3], Cache_spm.cache_spm(g, m).lZ, 1e-9);
            assertEquals(c[3], Cache_rayint.cache_rayint(g, m).lZ, 1e-9);
        }
    }

    // ----------------------------------------------------------------------
    // cache_prob_spm reaches m_l == 0 on its own: it evaluates E at oner(m,l)
    // ----------------------------------------------------------------------

    @Test
    @Timeout(120)
    public void probSpmIsADistributionWhenAListHoldsOneItem() {
        Matrix[] profiles = {uniform(), skewed()};
        Matrix m = row(1, 3, 2); // list 0 sends cache_spm to m_0 = 0
        for (Matrix g : profiles) {
            Matrix prob = Cache_prob_spm.cache_prob_spm(g, m);
            for (int i = 0; i < prob.getNumRows(); i++) {
                double rowsum = 0.0;
                for (int j = 0; j < prob.getNumCols(); j++) {
                    assertTrue(prob.get(i, j) >= 0.0 && prob.get(i, j) <= 1.0,
                            "probability out of range at (" + i + "," + j + "): " + prob.get(i, j));
                    rowsum += prob.get(i, j);
                }
                assertEquals(1.0, rowsum, 1e-9);
            }
            // sum_i P(item i in list l) = m_l holds exactly for E, closely for SPM
            for (int l = 0; l < 3; l++) {
                double occ = 0.0;
                for (int i = 0; i < prob.getNumRows(); i++) { occ += prob.get(i, l + 1); }
                assertEquals(m.get(l), occ, 0.20 * m.get(l));
            }
        }
    }
}
