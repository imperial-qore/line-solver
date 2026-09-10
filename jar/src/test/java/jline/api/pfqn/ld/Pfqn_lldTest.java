/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.ld;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Pfqn_lld must reproduce Pfqn_gld BIT FOR BIT.
 *
 * Saturating the rate shift at the threshold returns the same value, because
 * past it a further shift leaves the row unchanged over the columns the
 * recursion can still read; memoising the resulting repeated state changes what
 * is recomputed, never what is computed. Every terminal case is delegated back
 * to Pfqn_gld on the materialised block, so exact equality is the right
 * assertion and a tolerance would hide the regression these tests exist for.
 */
public class Pfqn_lldTest {

    private static Matrix mat(double[][] v) {
        Matrix m = new Matrix(v.length, v[0].length);
        for (int i = 0; i < v.length; i++)
            for (int j = 0; j < v[i].length; j++) m.set(i, j, v[i][j]);
        return m;
    }

    private static Matrix row(double... v) {
        Matrix m = new Matrix(1, v.length);
        for (int i = 0; i < v.length; i++) m.set(0, i, v[i]);
        return m;
    }

    /** Rate matrices spanning every threshold regime the saturation handles. */
    private static Matrix rates(int M, int Nt, int regime) {
        Matrix mu = new Matrix(M, Nt);
        for (int k = 0; k < M; k++) {
            for (int j = 1; j <= Nt; j++) {
                double v;
                switch (regime) {
                    case 0: v = Math.min(j, 2 + (k % 3)); break;          // multiserver
                    case 1: v = 1.0; break;                               // load independent
                    case 2: v = j; break;                                 // never settles
                    case 3: v = (k % 3 == 0) ? j : (k % 3 == 1 ? 1.0 : Math.min(j, 2)); break;
                    case 4: v = (j < 4) ? j : 4.0 + k; break;             // settles late
                    default: v = 0.5 + ((k * 7 + j * 13) % 7) * 0.25; break; // arbitrary
                }
                mu.set(k, j - 1, v);
            }
        }
        return mu;
    }

    @Test
    public void matchesGldBitwiseAcrossRegimes() {
        double[][] demands = {
                {2.0, 1.0, 3.0}, {1.0, 3.0, 0.5}, {4.0, 2.0, 1.5}, {0.7, 1.1, 2.2}};
        for (int regime = 0; regime <= 5; regime++) {
            for (int M = 1; M <= 4; M++) {
                for (int R = 1; R <= 3; R++) {
                    for (int Ntot : new int[]{1, 3, 6}) {
                        double[][] Lv = new double[M][R];
                        for (int i = 0; i < M; i++)
                            for (int r = 0; r < R; r++) Lv[i][r] = demands[i][r];
                        Matrix L = mat(Lv);
                        double[] nv = new double[R];
                        for (int r = 0; r < R; r++) nv[r] = Ntot / R;
                        nv[0] += Ntot - (Ntot / R) * R;
                        Matrix N = row(nv);
                        Matrix mu = rates(M, Ntot, regime);
                        double gld = Pfqn_gld.pfqn_gld(L, N, mu, null).G;
                        double lld = Pfqn_lld.pfqn_lld(L, N, mu, null).G;
                        assertEquals(gld, lld, 0.0,
                                "regime=" + regime + " M=" + M + " R=" + R + " Ntot=" + Ntot);
                    }
                }
            }
        }
    }

    /** A null rate matrix must default the same way in both routines. */
    @Test
    public void nullRatesDefaultIdentically() {
        Matrix L = mat(new double[][]{{2.0, 1.0}, {1.0, 3.0}});
        Matrix N = row(2.0, 2.0);
        assertEquals(Pfqn_gld.pfqn_gld(L, N, null, null).G,
                Pfqn_lld.pfqn_lld(L, N, null, null).G, 0.0);
    }

    /** A delay row never settles and must still give the reference value. */
    @Test
    public void delayRowIsHandled() {
        for (int Ntot : new int[]{3, 6}) {
            Matrix mu = new Matrix(3, Ntot);
            for (int j = 1; j <= Ntot; j++) {
                mu.set(0, j - 1, Math.min(j, 2));
                mu.set(1, j - 1, Math.min(j, 3));
                mu.set(2, j - 1, j);               // delay, s = Ntot
            }
            Matrix L = mat(new double[][]{{2.0, 1.0}, {1.0, 3.0}, {0.5, 0.5}});
            Matrix N = row(Ntot / 2, Ntot - Ntot / 2);   // integer populations
            assertEquals(Pfqn_gld.pfqn_gld(L, N, mu, null).G,
                    Pfqn_lld.pfqn_lld(L, N, mu, null).G, 0.0, "Ntot=" + Ntot);
        }
    }

    /** Zero population is the empty product in both. */
    @Test
    public void zeroPopulation() {
        Matrix L = mat(new double[][]{{2.0, 1.0}, {1.0, 3.0}});
        Matrix N = row(0.0, 0.0);
        Matrix mu = rates(2, 1, 0);
        assertEquals(Pfqn_gld.pfqn_gld(L, N, mu, null).G,
                Pfqn_lld.pfqn_lld(L, N, mu, null).G, 0.0);
    }
}
