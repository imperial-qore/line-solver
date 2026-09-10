/**
 * @file Cache Miss Analysis via Mean Value Analysis
 *
 * Implements Mean Value Analysis (MVA) algorithms for computing cache miss
 * probabilities in multi-level cache systems.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import org.apache.commons.math3.util.FastMath;

import jline.util.Pair;
import jline.util.matrix.Matrix;

public final class Cache_mva_miss {
    private Cache_mva_miss() {}

    /**
     * Compute cache miss probabilities using Mean Value Analysis approach.
     *
     * @param p Popularity vector
     * @param m Cache sizes vector
     * @param R Routing matrix
     * @return Pair of overall miss rate M and per-item miss probability Mk
     */
    public static Pair<Double, Matrix> cache_mva_miss(Matrix p, Matrix m, Matrix R) {
        int n = p.length();
        int h = m.length();

        if (m.elementSum() == 0.0 || m.elementMin() < 0.0) {
            Matrix Mk = Matrix.ones(1, n);
            double M = p.mult(Mk.transpose()).get(0, 0);
            return new Pair<Double, Matrix>(M, Mk);
        }

        Matrix w = Matrix.zeros(n, h);

        for (int j = 0; j < h; j++) {
            Matrix onerM = Matrix.oner(m, j);
            Pair<Double, Matrix> rec = cache_mva_miss(p, onerM, R);
            Matrix Mj = rec.getRight();

            for (int k = 0; k < n; k++) {
                double prod = 1.0;
                for (int i = 0; i <= j; i++) {
                    prod *= R.get(i, k);
                }
                double pPowJ = FastMath.pow(p.get(k), (double) (j + 1));
                w.set(k, j, prod * pPowJ * Math.abs(Mj.get(0, k)));
            }
        }

        double[] x = new double[h];
        for (int j = 0; j < h; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += Math.abs(w.get(k, j));
            }
            x[j] = 1.0 / sum;
        }

        Matrix Mk = Matrix.zeros(1, n);
        for (int k = 0; k < n; k++) {
            double value = 1.0;
            for (int j = 0; j < h; j++) {
                value -= x[j] * m.get(j) * w.get(k, j);
            }
            Mk.set(0, k, Math.abs(value));
        }

        double M = p.mult(Mk.transpose()).get(0, 0);
        return new Pair<Double, Matrix>(M, Mk);
    }
}
