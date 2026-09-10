/**
 * @file Exact Recursive Cache Analysis
 *
 * Implements exact recursive (EREC) algorithms for cache system analysis.
 * Provides numerically exact solutions for cache performance metrics in
 * systems where computational complexity allows for precise evaluation.
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import jline.io.InputOutput;
import jline.util.matrix.Matrix;

public final class Cache_erec {
    private Cache_erec() {}

    /**
     * Computes the cache miss rate using an exact recursive method.
     *
     * @param gamma Matrix representing the cache access factors.
     * @param m     Matrix representing the cache capacity vector.
     * @return Matrix containing the computed cache miss rates.
     */
    public static Matrix cache_erec(Matrix gamma, Matrix m) {
        return cache_erec_aux(gamma, m, gamma.getNumRows());
    }

    /**
     * Computes the normalizing constant under per-list storage cost caps,
     * E(m,k) = E_i(m,k) + sum_j m_j gamma_ij E_i(m-1_j, k-sigma_i 1_j), with the
     * extra boundary E(m,k)=0 whenever some residual cap is negative. See
     * Casale-Gast, IEEE/ACM Trans. Networking 29(2), 2021, Sec. IX.
     *
     * @param gamma Matrix representing the cache access factors.
     * @param m     Matrix representing the cache capacity vector.
     * @param sigma Item storage costs (sizes), one per item; null or empty for none.
     * @param k     Per-list storage cost caps, one per list; null or empty for none.
     * @return Matrix containing the computed normalizing constant.
     */
    public static Matrix cache_erec(Matrix gamma, Matrix m, Matrix sigma, Matrix k) {
        if (sigma == null || k == null || sigma.isEmpty() || k.isEmpty()) {
            return cache_erec(gamma, m);
        }
        return Matrix.singleton(cache_erec_cost(gamma, m, sigma, k));
    }

    /**
     * Dynamic program over the (residual capacity, residual cost cap) lattice.
     */
    private static double cache_erec_cost(Matrix gamma, Matrix m, Matrix sigma, Matrix k) {
        int n = gamma.getNumRows();
        int h = m.length();
        if (sigma.length() != n) {
            InputOutput.line_error(InputOutput.mfilename(new Object()),
                    "The item size vector must have one entry per item.");
        }
        if (k.length() != h) {
            InputOutput.line_error(InputOutput.mfilename(new Object()),
                    "The cost cap vector must have one entry per cache list.");
        }
        int[] mv = new int[h];
        int[] kv = new int[h];
        for (int j = 0; j < h; j++) {
            mv[j] = (int) Math.round(m.get(j));
            kv[j] = (int) Math.round(k.get(j));
            if (mv[j] < 0) {
                return 0.0;
            }
            if (kv[j] < 0) {
                return 0.0;
            }
        }
        int[] sv = new int[n];
        int mtot = 0;
        for (int j = 0; j < h; j++) {
            mtot += mv[j];
        }
        for (int i = 0; i < n; i++) {
            sv[i] = (int) Math.round(sigma.get(i));
            if (sv[i] <= 0) {
                InputOutput.line_error(InputOutput.mfilename(new Object()),
                        "Item sizes must be positive integers.");
            }
        }
        if (mtot > n) {
            return 0.0;
        }
        if (mtot == 0) {
            return 1.0;
        }
        int[] dims = new int[2 * h];
        for (int j = 0; j < h; j++) {
            dims[j] = mv[j] + 1;
            dims[h + j] = kv[j] + 1;
        }
        long lattice = 1;
        for (int d = 0; d < 2 * h; d++) {
            lattice *= dims[d];
            if (lattice > 10000000L) {
                InputOutput.line_error(InputOutput.mfilename(new Object()),
                        "The cost-constrained normalizing constant lattice exceeds the exact method limit; use the sampling method.");
            }
        }
        int size = (int) lattice;
        int[] stride = new int[2 * h];
        stride[0] = 1;
        for (int d = 1; d < 2 * h; d++) {
            stride[d] = stride[d - 1] * dims[d - 1];
        }
        int[] sub = new int[2 * h];
        double[] F = new double[size];
        double[] Fprev = new double[size];
        for (int idx = 0; idx < size; idx++) {
            decode(idx, stride, sub);
            int c = 0;
            for (int j = 0; j < h; j++) {
                c += sub[j];
            }
            F[idx] = (c == 0) ? 1.0 : 0.0;
        }
        for (int t = 1; t <= n; t++) {
            System.arraycopy(F, 0, Fprev, 0, size);
            for (int idx = 0; idx < size; idx++) {
                decode(idx, stride, sub);
                int c = 0;
                for (int j = 0; j < h; j++) {
                    c += sub[j];
                }
                if (c > t) {
                    F[idx] = 0.0;
                    continue;
                }
                double val = Fprev[idx];
                for (int j = 0; j < h; j++) {
                    int mj = sub[j];
                    int kj = sub[h + j];
                    double g = gamma.get(t - 1, j);
                    if (mj > 0 && kj >= sv[t - 1] && g != 0.0) {
                        val += g * mj * Fprev[idx - stride[j] - sv[t - 1] * stride[h + j]];
                    }
                }
                F[idx] = val;
            }
        }
        return F[size - 1];
    }

    private static void decode(int idx, int[] stride, int[] sub) {
        int rem = idx;
        for (int d = stride.length - 1; d >= 0; d--) {
            sub[d] = rem / stride[d];
            rem -= sub[d] * stride[d];
        }
    }

    /**
     * Auxiliary method for computing the cache miss rate using an exact recursive method.
     *
     * @param gamma Matrix representing the cache access factors.
     * @param m     Matrix representing the cache capacity vector.
     * @param k     Integer representing the current number of rows in the recursive step.
     * @return Matrix containing the computed cache miss rates.
     */
    public static Matrix cache_erec_aux(Matrix gamma, Matrix m, int k) {
        int h = m.getNumElements();

        if (m.elementSum() == 0.0) {
            return Matrix.singleton(1.0);
        }

        if (m.elementSum() > k || m.elementMin() < 0) {
            return Matrix.singleton(0.0);
        }

        if (k == 1 && m.elementSum() == 1.0) {
            int j = (int) m.find().value(); // Find the index of the non-zero element in m
            return Matrix.singleton(gamma.get(0, j));
        }

        Matrix E = cache_erec_aux(gamma, m, k - 1);
        for (int j = 0; j < h; j++) {
            if (m.get(j) > 0) {
                Matrix onerM = Matrix.oner(m, j);
                Matrix term = cache_erec_aux(gamma, onerM, k - 1).scale(gamma.get(k - 1, j) * m.get(j));
                E = E.add(term);
            }
        }
        return E;
    }
}
