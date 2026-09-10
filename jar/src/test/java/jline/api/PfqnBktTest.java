package jline.api;

import jline.api.pfqn.nc.Pfqn_ca;
import jline.api.pfqn.nc.Pfqn_kt;
import jline.api.pfqn.nc.Pfqn_bkt;
import jline.api.pfqn.nc.Pfqn_le;
import jline.api.pfqn.nc.Pfqn_ble;
import jline.util.matrix.Matrix;
import org.apache.commons.math3.special.Gamma;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * BKT subtracts the exact Stirling remainder s(N) = log(N!) - (N log N - N + log(2 pi N)/2)
 * of every class direction pfqn_kt Laplaces. With a think time the corrected expansion
 * is the SAME estimator as pfqn_ble (LE corrected by M units): the two saddle points
 * are one point in dual coordinates and Sylvester's identity exchanges the R x R Hessian
 * determinant for the M x M one. Without a think time they differ by the single
 * constant kappa - r(N+M), the remainder of the radial Gamma(N+M) direction that LE
 * integrates exactly and KT does not.
 */
public class PfqnBktTest {

    private static final double KAPPA = 1.0 - Math.log(2 * Math.PI) / 2;

    /** Stirling remainder of a Gamma(a) direction, r(1) = kappa. */
    private static double r(double a) {
        return Gamma.logGamma(a) - (a - 0.5) * Math.log(a) + a - 0.5 * Math.log(2 * Math.PI);
    }

    private static Matrix demands(int M, int R, long seed) {
        Matrix L = new Matrix(M, R);
        long x = seed;
        for (int i = 0; i < M; i++) {
            for (int c = 0; c < R; c++) {
                x = (x * 6364136223846793005L + 1442695040888963407L);
                double u = ((x >>> 11) / (double) (1L << 53));
                L.set(i, c, 0.5 + u);
            }
        }
        return L;
    }

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) m.set(0, i, v[i]);
        return m;
    }

    @Test
    public void theRemainderAtOneIsTheBleConstant() {
        assertEquals(KAPPA, Pfqn_bkt.stirlingRemainder(1.0), 1e-15);
        assertEquals(1.0 / (12 * 200.0), Pfqn_bkt.stirlingRemainder(200.0), 1e-7);
    }

    @Test
    public void theCorrectionIsTheSumOfTheClassRemainders() {
        int[][] shapes = {{2, 1}, {3, 2}, {5, 3}, {4, 4}};
        for (int[] s : shapes) {
            int M = s[0], R = s[1];
            Matrix L = demands(M, R, 10L * M + R);
            Matrix N = new Matrix(1, R);
            double sum = 0.0;
            for (int c = 0; c < R; c++) {
                N.set(0, c, 2.0 + c);
                sum += Pfqn_bkt.stirlingRemainder(2.0 + c);
            }
            Matrix[] Zs = {new Matrix(1, R), new Matrix(1, R)};
            for (int c = 0; c < R; c++) Zs[1].set(0, c, 3.0);
            for (Matrix Z : Zs) {
                double gap = Pfqn_bkt.pfqn_bkt(L, N, Z).lG - Pfqn_kt.pfqn_kt(L, N, Z).lG;
                assertEquals(-sum, gap, 1e-12);
            }
        }
    }

    @Test
    public void withAThinkTimeBktIsBle() {
        // Proposition: for Z > 0 the two corrected expansions coincide; 1e-5 covers the
        // accuracy of the two saddle-point solvers, not of the identity.
        int[][] shapes = {{3, 2}, {4, 3}, {6, 2}};
        for (int[] s : shapes) {
            int M = s[0], R = s[1];
            Matrix L = demands(M, R, 300L + M);
            Matrix N = new Matrix(1, R), Z = new Matrix(1, R);
            for (int c = 0; c < R; c++) {
                N.set(0, c, 6.0);
                Z.set(0, c, 2.5);
            }
            double d = Pfqn_bkt.pfqn_bkt(L, N, Z).lG - Pfqn_ble.pfqn_ble(L, N, Z).lG;
            assertTrue(Math.abs(d) < 1e-5, "M=" + M + " R=" + R + " |bkt - ble| = " + d);
        }
    }

    @Test
    public void withoutAThinkTimeTheGapToLeIsAKnownConstant() {
        int[][] shapes = {{3, 2}, {4, 3}, {6, 2}};
        for (int[] s : shapes) {
            int M = s[0], R = s[1];
            Matrix L = demands(M, R, 400L + M);
            Matrix N = new Matrix(1, R), Z = new Matrix(1, R);
            for (int c = 0; c < R; c++) N.set(0, c, 6.0);
            double eta = 6.0 * R + M;
            double gap = Pfqn_bkt.pfqn_bkt(L, N, Z).lG - Pfqn_le.pfqn_le(L, N, Z).lG;
            assertEquals(M * KAPPA - r(eta), gap, 1e-5);
        }
    }

    @Test
    public void lightLoadAgainstExactConvolution() {
        // At Z = 100 the class directions are Poisson and the remainder is the whole error.
        Matrix L = demands(4, 2, 7L);
        Matrix N = row(8.0, 6.0), Z = row(100.0, 100.0);
        double exact = Pfqn_ca.pfqn_ca(L, N, Z).lG;
        double errKt = Math.abs(Pfqn_kt.pfqn_kt(L, N, Z).lG - exact);
        double errPp = Math.abs(Pfqn_bkt.pfqn_bkt(L, N, Z).lG - exact);
        assertTrue(errKt > 0.01, "kt error " + errKt);
        assertTrue(errPp < 1e-3 && errPp < errKt / 20, "bkt error " + errPp + " vs kt " + errKt);
    }

    @Test
    public void matchesMatlab() {
        Matrix L = new Matrix(3, 2);
        L.set(0, 0, 1.0); L.set(0, 1, 0.5);
        L.set(1, 0, 0.7); L.set(1, 1, 1.2);
        L.set(2, 0, 0.3); L.set(2, 1, 0.9);
        Matrix N = row(2.0, 3.0), Z = row(1.0, 0.5);
        // MATLAB pfqn_bkt(L,[2 3],[1 .5]) and pfqn_bkt(L,[2 3])
        assertEquals(4.77242756160185, Pfqn_bkt.pfqn_bkt(L, N, Z).lG, 1e-8);
        assertEquals(4.13902507926115, Pfqn_bkt.pfqn_bkt(L, N).lG, 1e-8);
        double[][] refs = {{10, 14.1321407982845}, {20, 26.141862103174}, {40, 50.201513735287}};
        for (double[] nr : refs) {
            assertEquals(nr[1], Pfqn_bkt.pfqn_bkt(L, row(nr[0], nr[0]), Z).lG, 1e-6);
        }
        Matrix Ls = new Matrix(3, 2);
        Ls.set(0, 0, 1.0); Ls.set(1, 0, 0.7); Ls.set(2, 0, 0.3); Ls.set(1, 1, 1.2);
        assertEquals(2.81538680478809, Pfqn_bkt.pfqn_bkt(Ls, row(3.0, 2.0)).lG, 1e-8);
        assertEquals(1.97269392277537, Pfqn_bkt.pfqn_bkt(L, row(3.0, 0.0), Z).lG, 1e-8);
    }

    @Test
    public void emptyAndSelfLoopingClassesCarryNoRemainder() {
        Matrix L = demands(3, 2, 11L);
        Matrix N = row(5.0, 0.0), Z0 = new Matrix(1, 2);
        assertEquals(-Pfqn_bkt.stirlingRemainder(5.0),
                Pfqn_bkt.pfqn_bkt(L, N, Z0).lG - Pfqn_kt.pfqn_kt(L, N, Z0).lG, 1e-12);
        // class 2 visits one station with no think time: extracted exactly by pfqn_kt
        Matrix Ls = new Matrix(3, 2);
        Ls.set(0, 0, 1.0); Ls.set(1, 0, 0.7); Ls.set(2, 0, 0.3);
        Ls.set(1, 1, 1.2);
        Matrix N2 = row(3.0, 2.0);
        assertEquals(-Pfqn_bkt.stirlingRemainder(3.0),
                Pfqn_bkt.pfqn_bkt(Ls, N2, Z0).lG - Pfqn_kt.pfqn_kt(Ls, N2, Z0).lG, 1e-12);
        // with a think time it is Laplaced like any other
        Matrix Z = row(0.0, 1.0);
        double both = Pfqn_bkt.stirlingRemainder(3.0) + Pfqn_bkt.stirlingRemainder(2.0);
        assertEquals(-both, Pfqn_bkt.pfqn_bkt(Ls, N2, Z).lG - Pfqn_kt.pfqn_kt(Ls, N2, Z).lG, 1e-12);
    }
}
