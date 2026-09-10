package jline.api;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import jline.api.pfqn.nc.Pfqn_bk;
import jline.api.pfqn.nc.Pfqn_ca;
import jline.api.pfqn.nc.Pfqn_kt;
import jline.io.Ret;
import jline.util.matrix.Matrix;

/**
 * Birman-Kogan asymptotics (Stochastic Models 8(3):543-563, 1992).
 *
 * The multiprogramming model of the paper's Table 3 is fully specified in print
 * -- J = 2 device groups of 10 and 50 stations, service times .015 and .030,
 * five jobs per chain, branching .50/.50 and .6/.4 alternating, cpu rate
 * mu_k = M*mu0 with M = 50 -- and the paper tabulates four algorithms on it, so
 * the expected values below pin the algorithms rather than this implementation.
 */
public class PfqnBkTest {

    private static final double[] TH = {1 / 0.015, 1 / 0.030};
    private static final int[] MG = {10, 50};
    private static final int MPAR = 50;

    private static double branch(int k, int j) {
        return (k % 2 == 0) ? 0.5 : (j == 0 ? 0.6 : 0.4);
    }

    /** The Table 3 network: K dedicated cpus followed by the two device groups. */
    private static Matrix multiprogramming(int K, double mu0) {
        int M = K + MG[0] + MG[1];
        Matrix L = new Matrix(M, K);
        for (int k = 0; k < K; k++) L.set(k, k, 1.0 / (MPAR * mu0));
        int row = K;
        for (int j = 0; j < 2; j++) {
            for (int c = 0; c < MG[j]; c++) {
                for (int k = 0; k < K; k++) L.set(row, k, branch(k, j) / (TH[j] * MG[j]));
                row++;
            }
        }
        return L;
    }

    private static Matrix pop(int K, double n) {
        Matrix N = new Matrix(1, K);
        for (int k = 0; k < K; k++) N.set(0, k, n);
        return N;
    }

    @Test
    public void saddlePointReproducesTable3Column() {
        // U_k = x_k^0 / mu_0k by Corollary 1, capped at one for a saturated chain.
        // x^0 does not depend on mu_0k, so the mu_0 = 4 row fixes every other row
        // of the column by a pure rescaling.
        int[] Ks = {2, 3, 4, 5};
        double[][] want = {{90.9, 95.0}, {83.2, 86.0}, {75.6, 76.9}, {69.7, 70.2}};
        for (int i = 0; i < Ks.length; i++) {
            int K = Ks[i];
            Matrix L = multiprogramming(K, 4.0);
            Ret.pfqnNc sp = Pfqn_bk.pfqn_bk(L, pop(K, 5.0), new Matrix(1, K));
            for (int k = 0; k < 2; k++) {
                double U = Math.min(sp.X.get(0, k) * L.get(k, k), 1.0) * 100;
                assertEquals(want[i][k], U, 0.05);
            }
        }
    }

    @Test
    public void algorithmOnePinsEverySaturatedChain() {
        Matrix L = multiprogramming(2, 2.0);
        Ret.pfqnNc sat = Pfqn_bk.pfqn_bk(L, pop(2, 5.0), new Matrix(1, 2));
        for (int k = 0; k < 2; k++) {
            assertTrue(sat.X.get(0, k) * L.get(k, k) >= 1.0,
                    "chain " + k + " must sit on the pole of its dedicated station");
        }
    }

    @Test
    public void loadConcealmentReproducesBothTable3Columns() {
        int[] Ks = {2, 3, 2, 5};
        double[] mus = {4.0, 4.0, 2.0, 1.5};
        double[][] mva = {{70.2, 72.3}, {67.2, 69.1}, {93.1, 94.0}, {95.9, 96.4}};
        double[][] ue = {{69.8, 72.0}, {66.8, 68.7}, {93.4, 94.4}, {96.3, 96.8}};
        for (int i = 0; i < Ks.length; i++) {
            Matrix L = multiprogramming(Ks[i], mus[i]);
            Matrix N = pop(Ks[i], 5.0);
            Matrix Z = new Matrix(1, Ks[i]);
            Ret.pfqnBkLc tm = Pfqn_bk.pfqn_bklc(L, N, Z, "mva", 1e-12, 2000);
            Ret.pfqnBkLc tu = Pfqn_bk.pfqn_bklc(L, N, Z, "ue", 1e-12, 2000);
            for (int k = 0; k < 2; k++) {
                assertEquals(mva[i][k], tm.X.get(0, k) * L.get(k, k) * 100, 0.2);
                assertEquals(ue[i][k], tu.X.get(0, k) * L.get(k, k) * 100, 0.3);
            }
        }
    }

    @Test
    public void saddlePointEqualsKtWithoutDedicatedStations() {
        // With a think time in every chain there are no dedicated stations to keep
        // out of the exponent, so Proposition 1 IS the multidimensional saddle
        // point that pfqn_kt already computes.
        Matrix L = new Matrix(3, 2);
        double[][] v = {{1.0, 0.6}, {0.8, 1.2}, {0.5, 0.9}};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 2; j++) L.set(i, j, v[i][j]);
        }
        Matrix N = new Matrix(1, 2);
        N.set(0, 0, 10.0);
        N.set(0, 1, 8.0);
        Matrix Z = new Matrix(1, 2);
        Z.set(0, 0, 2.0);
        Z.set(0, 1, 1.0);
        assertEquals(Pfqn_kt.pfqn_kt(L, N, Z).lG, Pfqn_bk.pfqn_bk(L, N, Z).lG, 1e-9);
    }

    @Test
    public void uniformExpansionIsExactOnOneSlowStationAgainstAGroup() {
        // The regime the uniform expansion is built for: a single dominant pole
        // against M >> 1 identical stations.
        Matrix L = new Matrix(13, 1);
        L.set(0, 0, 0.9);
        for (int i = 1; i < 13; i++) L.set(i, 0, 0.1);
        Matrix N = new Matrix(1, 1);
        N.set(0, 0, 20.0);
        Matrix Z = new Matrix(1, 1);
        double exact = Pfqn_ca.pfqn_ca(L, N, Z).lG;
        assertEquals(exact, Pfqn_bk.pfqn_bkue(L, 20.0, 0.0).lG, 1e-4);
        assertEquals(exact, Pfqn_bk.pfqn_bk(L, N, Z).lG, 1e-4);
    }

    @Test
    public void uniformExpansionDegeneratesWithoutAGroup() {
        Matrix L = new Matrix(3, 1);
        L.set(0, 0, 0.5);
        L.set(1, 0, 0.4);
        L.set(2, 0, 0.3);
        Matrix N = new Matrix(1, 1);
        N.set(0, 0, 60.0);
        Matrix Z = new Matrix(1, 1);
        Z.set(0, 0, 5.0);
        assertEquals(Pfqn_kt.pfqn_kt(L, N, Z).lG, Pfqn_bk.pfqn_bkue(L, 60.0, 5.0).lG, 1e-9);
    }
}
