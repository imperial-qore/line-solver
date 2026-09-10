/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.cache;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Ray (WKB) expansion of the cost-capped cache normalizing constant. The profile
 * is the one the derivation note validates against its own dynamic program:
 * gamma_i = 0.3 + 2.7 i/n, sizes {1,2,3} in equal thirds of the item index (so
 * gcd = 1 and the sizes are genuinely diverse), m = n/4, k = 1.8 m.
 */
public class CacheSpmSizeTest {

    private static double[][] profile(int n) {
        double[][] gamma = new double[n][1];
        for (int i = 0; i < n; i++) { gamma[i][0] = 0.3 + 2.7 * ((i + 1.0) / n); }
        return gamma;
    }

    private static double[] sizes(int n) {
        double[] sigma = new double[n];
        for (int i = 0; i < n; i++) { sigma[i] = 1 + (3 * i) / n; }
        return sigma;
    }

    private static Matrix row(double[] v) {
        Matrix out = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) { out.set(0, i, v[i]); }
        return out;
    }

    private static Matrix mat(double[][] v) {
        Matrix out = new Matrix(v.length, v[0].length);
        for (int i = 0; i < v.length; i++) {
            for (int j = 0; j < v[0].length; j++) { out.set(i, j, v[i][j]); }
        }
        return out;
    }

    private static double logFactorial(int m) {
        double s = 0;
        for (int i = 2; i <= m; i++) { s += Math.log(i); }
        return s;
    }

    @Test
    public void exactCostModeReproducesTheNotesTable() {
        int[] ns = {50, 100, 200, 400};
        int[] ms = {12, 25, 50, 100};
        double[] want = {27.234323, 57.620291, 118.438233, 240.843731};
        for (int t = 0; t < 4; t++) {
            int n = ns[t];
            int m = ms[t];
            int k = (int) Math.ceil(1.8 * m);
            Cache_spm_size.Result r = Cache_spm_size.cache_spm_size(
                    profile(n), new double[]{m}, sizes(n), new double[]{k}, Cache_spm_size.EXACT);
            assertEquals(want[t], r.logE - logFactorial(m), 1e-5);
        }
    }

    @Test
    public void cumulativeCapsTrackCacheErecAndDegenerateToSizeFree() {
        int n = 120;
        double[][] gamma = profile(n);
        double[] sigma = sizes(n);
        double[] m = {30};

        // binding cap: the expansion sits close to the exact recursion
        Cache_spm_size.Result bindres = Cache_spm_size.cache_spm_size(
                gamma, m, sigma, new double[]{54});
        double exactBind = Math.log(Cache_erec.cache_erec(
                mat(gamma), row(m), row(sigma), row(new double[]{54})).get(0));
        assertTrue(bindres.binding[0]);
        assertTrue(bindres.zeta[0] < 1.0);
        assertTrue(Math.abs(bindres.logE - exactBind) < 0.15, "binding cap error too large");
        assertEquals("spm-size", bindres.method);
        assertTrue(bindres.iter < 30, "the saddle should converge in a handful of Newton steps");

        // slack cap: zeta returns to 1 and the cost coordinate leaves the saddle
        Cache_spm_size.Result slack = Cache_spm_size.cache_spm_size(
                gamma, m, sigma, new double[]{200});
        assertFalse(slack.binding[0]);
        assertEquals(1.0, slack.zeta[0], 1e-12);
        assertTrue(slack.iter < 30, "the slack branch should converge too");

        // a cap below the cheapest m items admits no state at all
        Cache_spm_size.Result none = Cache_spm_size.cache_spm_size(
                gamma, m, sigma, new double[]{10});
        assertEquals(0.0, none.E, 0.0);
        assertEquals("boundary", none.method);
    }

    @Test
    public void sizeLatticeIsDividedOutAndUniformSizesAreExact() {
        int n = 120;
        double[][] gamma = profile(n);
        double[] sigma = sizes(n);
        double[] doubled = new double[n];
        for (int i = 0; i < n; i++) { doubled[i] = 2 * sigma[i]; }
        double[] m = {30};

        Cache_spm_size.Result base = Cache_spm_size.cache_spm_size(
                gamma, m, sigma, new double[]{54});
        // sizes 2/4/6 at cap 108 is the same problem; 109 snaps down to the same lattice point
        Cache_spm_size.Result twice = Cache_spm_size.cache_spm_size(
                gamma, m, doubled, new double[]{108});
        Cache_spm_size.Result odd = Cache_spm_size.cache_spm_size(
                gamma, m, doubled, new double[]{109});
        assertEquals(2, twice.span);
        assertEquals(base.logE, twice.logE, 1e-12);
        assertEquals(base.logE, odd.logE, 1e-12);

        // a single item size: the cost of the list is sigma*m identically, so the cap
        // carries no information and the 2h saddle would be singular (note Sec. 5)
        double[] flat = new double[n];
        for (int i = 0; i < n; i++) { flat[i] = 2; }
        Cache_spm_size.Result uniform = Cache_spm_size.cache_spm_size(
                gamma, m, flat, new double[]{60});
        assertEquals("uniform-size", uniform.method);
        assertEquals(1.0, uniform.zeta[0], 1e-12);
        Cache_spm_size.Result tight = Cache_spm_size.cache_spm_size(
                gamma, m, flat, new double[]{59});
        assertEquals(0.0, tight.E, 0.0);   // 2*30 = 60 > 59, no state fits
    }

    @Test
    public void occupancyAndMeanCostSatisfyTheSaddleConditions() {
        int n = 80;
        Cache_spm_size.Result r = Cache_spm_size.cache_spm_size(
                profile(n), new double[]{20}, sizes(n), new double[]{36});
        double[] sigma = sizes(n);
        double occ = 0;
        double cost = 0;
        for (int i = 0; i < n; i++) {
            assertTrue(r.pij[i][1] >= 0.0);
            assertTrue(r.pij[i][0] >= 0.0);
            assertEquals(1.0, r.pij[i][0] + r.pij[i][1], 1e-12);
            occ += r.pij[i][1];
            cost += sigma[i] * r.pij[i][1];
        }
        // the saddle conditions are exactly sum_i pi_i = m and sum_i sigma_i pi_i = k
        assertEquals(20.0, occ, 1e-9);
        assertEquals(36.0, cost, 1e-9);
        assertEquals(cost, r.K[0], 1e-12);
    }
}
