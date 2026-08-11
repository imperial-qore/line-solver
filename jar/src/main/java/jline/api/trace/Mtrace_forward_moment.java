package jline.api.trace;

import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

import java.util.TreeSet;

public final class Mtrace_forward_moment {
    private Mtrace_forward_moment() {}

    /**
     * Computes the forward moments of a marked trace.
     *
     * @param T the inter-arrival times
     * @param A the class labels
     * @param orders vector with the orders of the moments to compute
     * @param norm 0 to return F_{i,c} unnormalized; 1 (default) to normalize
     * @return the forward moments as a matrix F(c,i) = F_{i,c}
     */
    public static Matrix mtrace_forward_moment(double[] T, int[] A, int[] orders, int norm) {
        TreeSet<Integer> uniq = new TreeSet<Integer>();
        for (int v : A) uniq.add(v);
        Integer[] marks = uniq.toArray(new Integer[0]);
        int C = marks.length;

        double[][] M = new double[C][orders.length];

        for (int j = 0; j < orders.length; j++) {
            int k = orders[j];
            for (int c = 0; c < C; c++) {
                double sum = 0.0;
                int count = 0;

                for (int i = 0; i < T.length - 1; i++) {
                    if (A[i] == marks[c].intValue()) {
                        sum += FastMath.pow(T[i + 1], k);
                        count++;
                    }
                }

                M[c][j] = count > 0 ? sum / count : 0.0;

                if (norm != 0) {
                    int countMarksInPrefix = 0;
                    for (int i = 0; i < A.length - 1; i++) {
                        if (A[i] == marks[c].intValue()) countMarksInPrefix++;
                    }
                    if (countMarksInPrefix > 0) {
                        M[c][j] = M[c][j] * (double) (T.length - 1) / countMarksInPrefix;
                    }
                }
            }
        }

        return new Matrix(M);
    }

    public static Matrix mtrace_forward_moment(double[] T, int[] A, int[] orders) {
        return mtrace_forward_moment(T, A, orders, 1);
    }
}
