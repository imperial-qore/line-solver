/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.nc;

import jline.api.pfqn.ld.Pfqn_ncld;
import jline.io.Ret;
import jline.solvers.SolverOptions;
import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Saddle-tilted Edgeworth normalizing constant.
 *
 * The fixtures are MATLAB pfqn_nre values. The point of the method is the
 * tilt: on the same model the untilted contour of Pfqn_nrl is 9.3e-2 off the
 * exact constant and the tilted one 1.2e-4, so a port that drops the tilt or
 * the Edgeworth term still looks plausible in isolation and is caught only by
 * comparing the two errors.
 */
public class Pfqn_nreTest {

    private static final double TOL = 1e-9;

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) m.set(0, i, v[i]);
        return m;
    }

    private static Matrix mat(double[][] a) {
        Matrix m = new Matrix(a.length, a[0].length);
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) m.set(i, j, a[i][j]);
        }
        return m;
    }

    private static Matrix ones(int rows, int cols) {
        Matrix m = new Matrix(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) m.set(i, j, 1.0);
        }
        return m;
    }

    /** Two classes, three single-server stations and a delay. */
    @Test
    public void twoClassModelMatchesMatlabAndBeatsTheUntiltedContour() {
        SolverOptions options = new SolverOptions();
        Matrix L = mat(new double[][]{{1.0, 0.5}, {0.7, 1.2}, {0.3, 0.9}});
        Matrix N = row(2, 3);
        Matrix Z = row(1, 0.5);
        double nre = Pfqn_nre.pfqn_nre(L, N, Z, ones(3, 5), options);
        assertEquals(4.783188572136291, nre, TOL);
        double nrl = Pfqn_nrl.pfqn_nrl(L, N, Z, ones(3, 5), options);
        // MATLAB 4.338130313502534; this port's numerical Hessian puts it 5e-6
        // off (Python-native lands on the same 4.338125491314063), which is
        // immaterial next to the 4.5e-1 that separates it from nre here
        assertEquals(4.338130313502534, nrl, 1e-5);
        double exact = 4.783752287503922; // MATLAB pfqn_ncld, method exact
        assertTrue(Math.abs(nre - exact) < Math.abs(nrl - exact),
                "the saddle-tilted contour must be closer to the exact constant than the untilted one");
    }

    /**
     * One class quotients the torus down to dimension zero, so the routine
     * returns pfqn_gldsingle itself rather than any expansion of it.
     */
    @Test
    public void singleClassIsExact() {
        SolverOptions options = new SolverOptions();
        Matrix L = mat(new double[][]{{0.5}, {1.0 / 3.0}, {0.2}});
        double nre = Pfqn_nre.pfqn_nre(L, row(5), row(0), ones(3, 5), options);
        assertEquals(-1.995145480198095, nre, 1e-12);
    }

    /** Station 1 is a 2-server queue, mu(1,n) = min(n,2); station 2 single. */
    @Test
    public void loadDependentRateRowIsCarried() {
        SolverOptions options = new SolverOptions();
        Matrix L = mat(new double[][]{{1.0, 0.6}, {0.5, 1.1}});
        Matrix mu = ones(2, 8);
        for (int j = 1; j < 8; j++) mu.set(0, j, 2.0);
        Matrix N = row(4, 4);
        Matrix Z = row(0.5, 0.5);
        assertEquals(3.911235290463940, Pfqn_nre.pfqn_nre(L, N, Z, mu, options), TOL);
        // and the same value through the load-dependent dispatcher
        options.method = "nre";
        Ret.pfqnNc ret = Pfqn_ncld.pfqn_ncld(L, N, Z, mu, options);
        assertEquals(3.911235290463940, ret.lG, TOL);
    }
}
