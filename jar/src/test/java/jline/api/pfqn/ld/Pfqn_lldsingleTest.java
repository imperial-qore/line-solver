/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */

package jline.api.pfqn.ld;

import jline.util.matrix.Matrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Pfqn_lldsingle must reproduce Pfqn_gldsingle BIT FOR BIT.
 *
 * The capped sweep performs a subset of the full sweep's operations, never a
 * rearrangement of them, so the two constants agree in the last bit and not
 * merely to a tolerance. Every assertion here is exact equality for that
 * reason: a tolerance would hide precisely the regression these tests exist to
 * catch. The rate shapes span every threshold regime the cap has to handle,
 * from s_k = 1 (load independent) through s_k = N (a row that never settles).
 */
public class Pfqn_lldsingleTest {

    private static Matrix col(double... v) {
        Matrix m = new Matrix(v.length, 1);
        for (int i = 0; i < v.length; i++) m.set(i, 0, v[i]);
        return m;
    }

    /** Row k of the rate matrix, by regime index. */
    private static Matrix rates(int M, int Nt, int regime) {
        Matrix mu = new Matrix(M, Nt);
        for (int k = 0; k < M; k++) {
            for (int j = 1; j <= Nt; j++) {
                double v;
                switch (regime) {
                    case 0: v = 1.0; break;                              // s = 1
                    case 1: v = Math.min(j, 2 + (k % 3)); break;         // multiserver
                    case 2: v = j; break;                                // never settles
                    case 3: v = (k % 2 == 0) ? 1.0 : Math.min(j, 3); break;
                    case 4: v = (j < 4) ? j : 4.0 + k; break;            // settles late
                    default: v = 1.0 + ((k * 7 + j * 13) % 5) * 0.25; break; // arbitrary
                }
                mu.set(k, j - 1, v);
            }
        }
        return mu;
    }

    @Test
    public void matchesGldsingleBitwiseAcrossRegimes() {
        double[] demands = {0.7, 0.4, 0.9, 1.3, 0.2, 1.1, 0.55, 0.85, 0.15};
        for (int regime = 0; regime <= 5; regime++) {
            for (int M : new int[]{1, 2, 5, 9}) {
                for (int Nt : new int[]{1, 2, 5, 16, 33}) {
                    double[] d = new double[M];
                    System.arraycopy(demands, 0, d, 0, M);
                    Matrix L = col(d);
                    Matrix mu = rates(M, Nt, regime);
                    Matrix N = Matrix.singleton(Nt);
                    double lld = Pfqn_lldsingle.pfqn_lldsingle(L, N, mu, null).lG;
                    double gld = Pfqn_gldsingle.pfqn_gldsingle(L, N, mu, null).lG;
                    assertEquals(gld, lld, 0.0,
                            "regime=" + regime + " M=" + M + " N=" + Nt);
                }
            }
        }
    }

    /** A zero demand sends log(L) to -Inf; the cap must not turn that into NaN. */
    @Test
    public void zeroDemandStation() {
        for (int Nt : new int[]{5, 16}) {
            Matrix L = col(0.0, 0.7, 1.2, 0.4);
            Matrix mu = rates(4, Nt, 1);
            Matrix N = Matrix.singleton(Nt);
            assertEquals(Pfqn_gldsingle.pfqn_gldsingle(L, N, mu, null).lG,
                    Pfqn_lldsingle.pfqn_lldsingle(L, N, mu, null).lG, 0.0);
        }
    }

    /** Negative demands leave the log domain, which is a separate code path. */
    @Test
    public void negativeDemandUsesLinearBranch() {
        for (int Nt : new int[]{5, 16}) {
            Matrix L = col(-0.7, -1.3, -0.4);
            Matrix mu = rates(3, Nt, 1);
            Matrix N = Matrix.singleton(Nt);
            assertEquals(Pfqn_gldsingle.pfqn_gldsingle(L, N, mu, null).G,
                    Pfqn_lldsingle.pfqn_lldsingle(L, N, mu, null).G, 0.0);
        }
    }

    /** N = 0 is the empty product, as Pfqn_gldsingle returns from unrun loops. */
    @Test
    public void emptyPopulation() {
        Matrix L = col(1.0, 2.0);
        Matrix mu = new Matrix(2, 1);
        mu.set(0, 0, 1.0);
        mu.set(1, 0, 1.0);
        assertEquals(1.0, Pfqn_lldsingle.pfqn_lldsingle(L, Matrix.singleton(0), mu, null).G, 0.0);
        assertEquals(0.0, Pfqn_lldsingle.pfqn_lldsingle(L, Matrix.singleton(0), mu, null).lG, 0.0);
    }

    @Test
    public void multiclassIsRefused() {
        Matrix L = new Matrix(2, 3);
        Matrix mu = new Matrix(2, 4);
        assertThrows(RuntimeException.class,
                () -> Pfqn_lldsingle.pfqn_lldsingle(L, Matrix.singleton(4), mu, null));
    }

    /**
     * A row that ends in +Inf must not tie its FINITE entries to that tail.
     * With tail = Inf the bound ulp*max(abs(tail),1) is itself Inf and
     * abs(prev-Inf) &lt;= Inf holds for every finite prev, which would collapse
     * the whole row onto alpha = Inf and zero the constant. Pfqn_rd reaches
     * exactly this input, mapping NaN rates to Inf before calling in.
     */
    @Test
    public void infiniteTailDoesNotSwallowFiniteRates() {
        for (int Nt : new int[]{4, 9, 17}) {
            Matrix mu = new Matrix(2, Nt);
            for (int j = 1; j <= Nt; j++) {
                mu.set(0, j - 1, Math.min(j, 2));
                mu.set(1, j - 1, j < Nt - 1 ? (double) j : Double.POSITIVE_INFINITY);
            }
            Matrix L = col(0.7, 0.4);
            Matrix N = Matrix.singleton(Nt);
            assertEquals(Pfqn_gldsingle.pfqn_gldsingle(L, N, mu, null).lG,
                    Pfqn_lldsingle.pfqn_lldsingle(L, N, mu, null).lG, 0.0, "N=" + Nt);
        }
    }

    /** An entirely infinite row ties with itself and must reach s_k = 1. */
    @Test
    public void allInfiniteRowStillCollapses() {
        int Nt = 8;
        Matrix mu = new Matrix(2, Nt);
        for (int j = 1; j <= Nt; j++) {
            mu.set(0, j - 1, Math.min(j, 3));
            mu.set(1, j - 1, Double.POSITIVE_INFINITY);
        }
        Matrix L = col(0.7, 0.4);
        Matrix N = Matrix.singleton(Nt);
        assertEquals(Pfqn_gldsingle.pfqn_gldsingle(L, N, mu, null).lG,
                Pfqn_lldsingle.pfqn_lldsingle(L, N, mu, null).lG, 0.0);
    }

    /**
     * A FALSE tie would be a wrong answer, so the scan must not collapse a row
     * whose tail only nearly repeats.
     */
    @Test
    public void thresholdScanDoesNotMergeDistinctRates() {
        int Nt = 8;
        Matrix mu = new Matrix(1, Nt);
        for (int j = 0; j < Nt; j++) mu.set(0, j, 3.0);
        mu.set(0, Nt - 3, 3.0 * (1 + 1e-9));
        Matrix L = col(0.7);
        Matrix N = Matrix.singleton(Nt);
        assertEquals(Pfqn_gldsingle.pfqn_gldsingle(L, N, mu, null).lG,
                Pfqn_lldsingle.pfqn_lldsingle(L, N, mu, null).lG, 0.0);
    }
}
