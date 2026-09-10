/**
 * @file Convolutional Bound Hierarchy (Dowdy, Eager, Gordon, Saxton 1984)
 *
 * Level-parameterized throughput bounds: fill column c=M-level of Buzen's g-array with
 * BJB-derived estimates, then convolve the remaining `level` servers exactly (and the IS
 * think-time station exactly). Ported at parity from MATLAB pfqn_cbh.m.
 *
 * @since LINE 3.0
 */
package jline.api.pfqn.mva;

import jline.util.matrix.Matrix;

public final class Pfqn_cbh {
    private Pfqn_cbh() {}

    /**
     * Level-`level` convolutional bound hierarchy on throughput.
     *
     * @param L     service demand vector (M x 1)
     * @param N     population
     * @param Z     think time (folded as an exactly-convolved IS station)
     * @param level number of exactly-convolved servers, 1..M
     * @return {Xlo, Xhi} lower/upper throughput bounds
     */
    public static double[] pfqn_cbh(Matrix L, double N, double Z, int level) {
        int M = L.length();
        if (level < 1) level = 1;
        if (level > M) level = M;
        int c = Math.max(1, M - level);   // c=1 (single-server column) is already exact
        double Xlo = cbhHier(L, N, Z, c, false);
        double Xhi = cbhHier(L, N, Z, c, true);
        return new double[]{Xlo, Xhi};
    }

    private static double cbhHier(Matrix L, double N, double Z, int c, boolean upper) {
        int M = L.length();
        int Ni = (int) N;
        double Rc = 0.0, Lbc = 0.0, sum = 0.0;
        for (int i = 0; i < c; i++) {
            Rc += L.get(i);
            if (L.get(i) > Lbc) Lbc = L.get(i);
            sum += L.get(i);
        }
        double Lac = sum / c;
        double[] e = new double[Ni + 1];
        e[0] = 1.0;
        for (int i = 1; i <= Ni; i++) {
            double B = upper ? i / (Rc + (i - 1) * Lac) : i / (Rc + (i - 1) * Lbc);
            e[i] = e[i - 1] / B;
        }
        if (c == 1) {                     // single-server column is exact
            e[0] = 1.0;
            for (int i = 1; i <= Ni; i++) e[i] = e[i - 1] * L.get(0);
        }
        double[] g = new double[Ni + 1];
        System.arraycopy(e, 0, g, 0, Ni + 1);
        for (int m = c; m < M; m++) {     // exact convolution of servers c+1..M (0-based m>=c)
            for (int n = 1; n <= Ni; n++) {
                g[n] = g[n] + L.get(m) * g[n - 1];
            }
        }
        if (Z > 0) {                      // convolve the IS delay exactly: g_Z(j)=Z^j/j!
            double[] gd = new double[Ni + 1];
            double fact = 1.0;
            for (int j = 0; j <= Ni; j++) {
                if (j > 0) fact *= j;
                gd[j] = Math.pow(Z, j) / fact;
            }
            double[] gfull = new double[Ni + 1];
            for (int n = 0; n <= Ni; n++) {
                double acc = 0.0;
                for (int j = 0; j <= n; j++) acc += g[j] * gd[n - j];
                gfull[n] = acc;
            }
            g = gfull;
        }
        return g[Ni - 1] / g[Ni];
    }
}
