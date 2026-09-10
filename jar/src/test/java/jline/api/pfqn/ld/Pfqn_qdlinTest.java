/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.ld;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * pfqn_qdlin: the Linearizer arm of AMVA-LD on a plain demand matrix.
 *
 * <p>Every golden here was measured with the native-Python pfqn_qdlin, which was itself validated
 * against {@code SolverMVA(model,'qdlin')} on 640 random closed models (single server, multiserver,
 * with and without think time, load dependent) with every metric agreeing to 3e-16 relative. The
 * MATLAB and C++ twins reproduce the same figures to all ten printed digits. So these numbers are
 * the solver's, not an independent approximation of it.
 *
 * <p>The ITERATION COUNTS are pinned alongside the values, and are the sharper check: a count
 * agrees only if the convergence test, the initial guess and every intermediate iterate agree too.
 *
 * <p>Two of the checks pin PROPERTIES rather than numbers, and those are the ones a bad
 * transcription trips: the reported utilization switches between the iterated Uchain and the
 * analytic T*S/c depending on whether the model carries lld scaling, and mu and nservers are
 * different mechanisms rather than two spellings of one station.
 */
public class Pfqn_qdlinTest {

    /** The convergence tolerance the fixed point stops on; goldens are only meaningful with it. */
    private static final double TOL = 1e-9;

    private static Matrix mat(double[][] a) {
        Matrix m = new Matrix(a.length, a[0].length);
        for (int i = 0; i < a.length; i++) {
            for (int j = 0; j < a[0].length; j++) {
                m.set(i, j, a[i][j]);
            }
        }
        return m;
    }

    private static Matrix row(double... a) {
        Matrix m = new Matrix(1, a.length);
        for (int j = 0; j < a.length; j++) {
            m.set(0, j, a[j]);
        }
        return m;
    }

    private static Matrix col(double... a) {
        Matrix m = new Matrix(a.length, 1);
        for (int j = 0; j < a.length; j++) {
            m.set(j, 0, a[j]);
        }
        return m;
    }

    /** Two stations, one class, demands 0.5 and 0.3, three jobs, think time one. */
    private static Matrix oneClassDemands() {
        return mat(new double[][] {{0.5}, {0.3}});
    }

    @Test
    public void oneClassMatchesTheSolver() {
        Pfqn_qdlin.Result r =
                Pfqn_qdlin.pfqn_qdlin(oneClassDemands(), row(3.0), row(1.0), null, null);
        assertEquals(1.1055055966, r.Q.get(0, 0), TOL);
        assertEquals(0.5464847591, r.Q.get(1, 0), TOL);
        assertEquals(0.6740047306, r.U.get(0, 0), TOL);
        assertEquals(0.4044028384, r.U.get(1, 0), TOL);
        assertEquals(0.8201022533, r.R.get(0, 0), TOL);
        assertEquals(1.3480094612, r.X.get(0, 0), TOL);
        assertEquals(2.2255034473, r.C.get(0, 0), TOL);
        assertEquals(397, r.iter);
    }

    @Test
    public void twoClassesAndThreeStationsMatchTheSolver() {
        Matrix L = mat(new double[][] {{0.4, 0.9}, {0.6, 0.1}, {0.2, 0.5}});
        Matrix N = row(4.0, 3.0);
        Matrix Z = row(2.0, 1.0);
        Pfqn_qdlin.Result r = Pfqn_qdlin.pfqn_qdlin(L, N, Z, null, null);
        assertEquals(1.2023471976, r.Q.get(0, 0), TOL);
        assertEquals(0.8455406444, r.Q.get(1, 0), TOL);
        assertEquals(0.2939223883, r.Q.get(2, 0), TOL);
        assertEquals(1.7138257835, r.Q.get(0, 1), TOL);
        assertEquals(0.1321216124, r.Q.get(1, 1), TOL);
        assertEquals(0.5172388252, r.Q.get(2, 1), TOL);
        assertEquals(0.8290948878, r.X.get(0, 0), TOL);
        assertEquals(0.6368140023, r.X.get(0, 1), TOL);
        assertEquals(853, r.iter);

        // Little's law over the whole network, to the tolerance the fixed point was
        // actually converged to rather than to machine precision.
        for (int c = 0; c < 2; c++) {
            double q = 0.0;
            for (int i = 0; i < 3; i++) {
                q += r.Q.get(i, c);
            }
            assertEquals(N.get(c), q + r.X.get(0, c) * Z.get(c), 1e-5);
            // The cycle time is the residence times summed, plus the think time.
            double res = 0.0;
            for (int i = 0; i < 3; i++) {
                res += r.R.get(i, c);
            }
            assertEquals(r.C.get(0, c), res + Z.get(c), 1e-5);
        }
    }

    @Test
    public void withoutLldScalingTheReportedUtilizationIsAnalytic() {
        Matrix L = oneClassDemands();
        Pfqn_qdlin.Result r = Pfqn_qdlin.pfqn_qdlin(L, row(3.0), row(1.0), null, col(2.0, 1.0));
        assertEquals(0.7844085389, r.Q.get(0, 0), TOL);
        assertEquals(0.6467257168, r.Q.get(1, 0), TOL);
        assertEquals(1.5688658520, r.X.get(0, 0), TOL);
        assertEquals(300, r.iter);
        // sn_deaggregate_chain_results recomputes T*S/c from the NOMINAL demand when the
        // model carries no lld, cd or jd scaling.
        assertEquals(r.X.get(0, 0) * L.get(0, 0) / 2.0, r.U.get(0, 0), 1e-12);
        assertEquals(r.X.get(0, 0) * L.get(1, 0), r.U.get(1, 0), 1e-12);
    }

    @Test
    public void muAndNserversAreDifferentMechanisms() {
        Matrix L = oneClassDemands();
        Pfqn_qdlin.Result viaServers =
                Pfqn_qdlin.pfqn_qdlin(L, row(3.0), row(1.0), null, col(2.0, 1.0));
        // The same two-server station spelled as a load-dependent rate lattice.
        Matrix mu = mat(new double[][] {{1.0, 2.0, 2.0, 2.0}, {1.0, 1.0, 1.0, 1.0}});
        Pfqn_qdlin.Result viaMu = Pfqn_qdlin.pfqn_qdlin(L, row(3.0), row(1.0), mu, null);
        assertEquals(0.7844030387, viaMu.Q.get(0, 0), TOL);
        assertEquals(300, viaMu.iter);

        // A load-dependent model reports the ITERATED utilization instead, because the
        // analyzer forwards Uchain to the deaggregation only under lld/cd/jd scaling.
        assertEquals(0.5088714356, viaMu.U.get(0, 0), TOL);
        assertEquals(0.4705485024, viaMu.U.get(1, 0), TOL);
        assertNotEquals(viaMu.X.get(0, 0) * L.get(0, 0), viaMu.U.get(0, 0), 1e-6);

        // Close, because both describe the same station, but NOT equal: the server count
        // goes through the softmin and the lattice through the interpolation.
        assertNotEquals(viaServers.Q.get(0, 0), viaMu.Q.get(0, 0));
        assertTrue(Math.abs(viaServers.Q.get(0, 0) - viaMu.Q.get(0, 0)) < 1e-4);
    }

    @Test
    public void thinkTimeMayBeAbsentAndAClassMayBeEmpty() {
        Matrix L = mat(new double[][] {{0.5, 0.2}, {0.3, 0.7}});
        Matrix N = row(3.0, 2.0);
        // No delay station is appended when there is no think time.
        Pfqn_qdlin.Result nz = Pfqn_qdlin.pfqn_qdlin(L, N, row(0.0, 0.0), null, null);
        assertEquals(1.6035096548, nz.Q.get(0, 0), TOL);
        assertEquals(1.4989846125, nz.Q.get(1, 1), TOL);
        assertEquals(1.3103309302, nz.X.get(0, 0), TOL);
        for (int c = 0; c < 2; c++) {
            double q = 0.0;
            for (int i = 0; i < 2; i++) {
                q += nz.Q.get(i, c);
            }
            assertEquals(N.get(c), q, 1e-5);
        }

        // An empty population is answered with zeros, not refused, and no sweep runs.
        Pfqn_qdlin.Result e = Pfqn_qdlin.pfqn_qdlin(L, row(0.0, 0.0), row(0.0, 0.0), null, null);
        assertEquals(0, e.iter);
        for (int i = 0; i < 2; i++) {
            for (int c = 0; c < 2; c++) {
                assertEquals(0.0, e.Q.get(i, c), 0.0);
            }
        }

        // One class present, one empty: the empty one stays at zero throughout.
        Pfqn_qdlin.Result p = Pfqn_qdlin.pfqn_qdlin(L, row(2.0, 0.0), row(0.0, 0.0), null, null);
        assertEquals(0.0, p.X.get(0, 1), 0.0);
        for (int i = 0; i < 2; i++) {
            assertEquals(0.0, p.Q.get(i, 1), 0.0);
        }
        assertTrue(p.X.get(0, 0) > 0.0);
    }

    @Test
    public void theRefusalsAreByName() {
        Matrix L = mat(new double[][] {{0.5, 0.2}, {0.3, 0.7}});
        assertThrows(IllegalArgumentException.class,
                () -> Pfqn_qdlin.pfqn_qdlin(L, row(1.0, 1.0, 1.0), row(1.0, 1.0), null, null));
        assertThrows(IllegalArgumentException.class,
                () -> Pfqn_qdlin.pfqn_qdlin(L, row(1.0, 1.0), row(1.0, 1.0, 1.0), null, null));
        assertThrows(IllegalArgumentException.class,
                () -> Pfqn_qdlin.pfqn_qdlin(L, row(1.0, 1.0), row(1.0, 1.0), null,
                        col(1.0, 1.0, 1.0)));
        assertThrows(IllegalArgumentException.class,
                () -> Pfqn_qdlin.pfqn_qdlin(L, row(Double.POSITIVE_INFINITY, 1.0), row(1.0, 1.0),
                        null, null));
    }
}
