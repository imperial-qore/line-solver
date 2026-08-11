/**
 * @file Multi-class trace moment computation
 *
 * @since LINE 3.0
 */
package jline.api.trace;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

import jline.util.matrix.Matrix;

public final class Mtrace_moment {
    private Mtrace_moment() {}

    /**
     * Computes the empirical class-dependent moments of a multi-class trace.
     *
     * @param T      vector of inter-arrival times
     * @param A      vector of class labels
     * @param orders vector with the orders of the moments to compute
     * @param after  0 for Horvath, 1 for Bucholz variables
     * @param norm   0 to not normalize, 1 to normalize
     * @return matrix of moments
     */
    public static Matrix mtrace_moment(double[] T, int[] A, int[] orders, int after, int norm) {
        TreeSet<Integer> distinct = new TreeSet<Integer>();
        for (int a : A) distinct.add(a);
        List<Integer> marks = new ArrayList<Integer>(distinct);
        Collections.sort(marks);
        int C = marks.size();

        double[][] M = new double[C][orders.length];

        for (int j = 0; j < orders.length; j++) {
            int k = orders[j];
            for (int c = 0; c < C; c++) {
                int mark = marks.get(c);

                if (after == 1) {
                    double sum = 0.0;
                    int count = 0;
                    for (int i = 0; i < T.length - 1; i++) {
                        if (A[i] == mark) {
                            sum += Math.pow(T[i + 1], k);
                            count++;
                        }
                    }
                    M[c][j] = (count > 0) ? sum / count : 0.0;

                    if (norm == 1) {
                        int totalCount = A.length - 1;
                        int classCount = 0;
                        for (int i = 0; i < A.length - 1; i++) {
                            if (A[i] == mark) classCount++;
                        }
                        if (classCount > 0) {
                            M[c][j] = M[c][j] * (double) totalCount / classCount;
                        }
                    }
                } else {
                    double sum = 0.0;
                    int count = 0;
                    for (int i = 0; i < T.length; i++) {
                        if (A[i] == mark) {
                            sum += Math.pow(T[i], k);
                            count++;
                        }
                    }
                    M[c][j] = (count > 0) ? sum / count : 0.0;

                    if (norm == 1) {
                        int totalCount = A.length;
                        int classCount = 0;
                        for (int a : A) {
                            if (a == mark) classCount++;
                        }
                        if (classCount > 0) {
                            M[c][j] = M[c][j] * (double) totalCount / classCount;
                        }
                    }
                }
            }
        }

        return new Matrix(M);
    }

    public static Matrix mtrace_moment(double[] T, int[] A, int[] orders) {
        return mtrace_moment(T, A, orders, 0, 0);
    }

    public static Matrix mtrace_moment(double[] T, int[] A, int[] orders, int after) {
        return mtrace_moment(T, A, orders, after, 0);
    }
}
