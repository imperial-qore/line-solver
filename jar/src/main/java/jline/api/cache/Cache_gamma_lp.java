/**
 * @file Cache Access Factor Computation via Linear Programming
 *
 * @since LINE 3.0
 */
package jline.api.cache;

import java.util.ArrayList;

import jline.io.Ret;
import jline.util.matrix.Matrix;

public final class Cache_gamma_lp {
    private Cache_gamma_lp() {}

    /**
     * Computes access factors for the cache.
     */
    public static Ret.cacheGamma cache_gamma_lp(Matrix[] lambda, Matrix[][] R) {
        int u = lambda.length;
        int n = lambda[0].getNumRows();
        int h = lambda[0].getNumCols() - 1;
        Matrix gamma = new Matrix(n, h);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < h; j++) {
                Matrix Rvi = new Matrix(R[0][i].getNumRows(), R[0][i].getNumCols());
                for (int v = 0; v < u; v++) {
                    Rvi = Rvi.add(1.0, R[v][i]);
                }
                ArrayList<Integer> Pij = new ArrayList<Integer>();
                Pij.add(1 + j);
                ArrayList<Integer> pr_j = cache_par(Rvi, 1 + j);
                while (pr_j.size() > 0) {
                    Pij.add(0, pr_j.get(0));
                    pr_j = cache_par(Rvi, pr_j.get(0));
                }
                if (Pij.size() == 0) {
                    gamma.set(i, j, 0);
                } else {
                    gamma.set(i, j, 1);
                    for (int li = 1; li < Pij.size(); li++) {
                        double y = 0.0;
                        int l_1 = Pij.get(li - 1);
                        int l = Pij.get(li);
                        for (int v = 0; v < u; v++) {
                            for (int t = 0; t <= l_1; t++) {
                                y += lambda[v].get(i, t) * R[v][i].get(t, l);
                            }
                        }
                        gamma.set(i, j, gamma.get(i, j) * y);
                    }
                }
            }
        }
        int[] parent = new int[h];
        // tree structure read off item 0's routing matrix aggregated over users,
        // the same matrix the gamma loop walks
        Matrix Rtot = new Matrix(R[0][0].getNumRows(), R[0][0].getNumCols());
        for (int v = 0; v < u; v++) {
            Rtot = Rtot.add(1.0, R[v][0]);
        }
        for (int j = 0; j < h; j++) {
            ArrayList<Integer> pj = cache_par(Rtot, 1 + j);
            parent[j] = pj.isEmpty() ? -1 : (pj.get(0) - 1); // list indices, -1 = miss list
        }
        return new Ret.cacheGamma(gamma, u, n, h, parent);
    }

    /**
     * Finds the parent of a given list index.
     */
    public static ArrayList<Integer> cache_par(Matrix R, int j) {
        ArrayList<Integer> parent = new ArrayList<Integer>();
        for (int i = 0; i <= j - 1; i++) {
            if (R.get(i, j) != 0.0) {
                parent.add(i);
            }
        }
        if (parent.size() > 1) {
            throw new RuntimeException("A cache has a list with more than one parent, but the structure must be a tree.");
        }
        return parent;
    }
}
