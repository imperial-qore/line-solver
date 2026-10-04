/*
 * Ported from BUTools-family fluid tools (G. Horvath).
 */
package jline.lib.butools.reptrans;

import java.util.Arrays;
import java.util.Comparator;

import jline.util.matrix.Matrix;

/**
 * Similarity transformation mapping a closing vector to a vector of ones.
 */
public final class TransformToOnes {
    private TransformToOnes() {}

    /**
     * Returns the similarity transformation matrix B such that B*clovec = ones.
     * Works even if clovec has zero entries.
     *
     * @param clovec the original closing (column) vector, shape (M,1)
     * @return B, shape (M,M), with B*clovec = ones(M,1)
     */
    public static Matrix transformToOnes(Matrix clovec) {
        int m = clovec.getNumRows() * clovec.getNumCols();
        double[] cv = new double[m];
        int k = 0;
        for (int i = 0; i < clovec.getNumRows(); i++)
            for (int j = 0; j < clovec.getNumCols(); j++)
                cv[k++] = clovec.get(i, j);

        // stable descending sort of clovec (ascending of -clovec)
        Integer[] ix = new Integer[m];
        for (int i = 0; i < m; i++) ix[i] = i;
        final double[] cvf = cv;
        Arrays.sort(ix, new Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                int c = Double.compare(cvf[b], cvf[a]);
                return c != 0 ? c : Integer.compare(a, b);
            }
        });

        Matrix P = Matrix.zeros(m, m);
        for (int i = 0; i < m; i++) P.set(i, ix[i], 1.0);

        double[] cp = new double[m];
        for (int i = 0; i < m; i++) cp[i] = cv[ix[i]];

        Matrix B = Matrix.zeros(m, m);
        for (int i = 0; i < m; i++) {
            double s = 0;
            for (int j = 0; j <= i; j++) s += cp[j];
            for (int j = 0; j <= i; j++) B.set(i, j, 1.0 / s);
        }
        return B.mult(P);
    }
}
