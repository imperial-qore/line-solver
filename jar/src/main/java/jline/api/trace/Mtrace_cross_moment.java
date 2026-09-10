package jline.api.trace;

import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

public final class Mtrace_cross_moment {
    private Mtrace_cross_moment() {}

    /**
     * Computes the k-th order moment of the inter-arrival time between an event
     * of class i and an event of class j, for all possible pairs of classes.
     *
     * @param T inter-arrival times
     * @param L class labels
     * @param k order of the moment
     * @return matrix where element (i,j) is the k-th order moment of the inter-arrival
     *         time between an event of class i and an event of class j
     */
    public static Matrix mtrace_cross_moment(double[] T, int[] L, int k) {
        TreeSet<Integer> distinctSet = new TreeSet<Integer>();
        for (int v : L) {
            distinctSet.add(v);
        }
        List<Integer> marks = new ArrayList<Integer>(distinctSet);
        Collections.sort(marks);
        int C = marks.size();

        double[][] MC = new double[C][C];
        int[][] count = new int[C][C];

        for (int t = 1; t < T.length; t++) {
            for (int i = 0; i < C; i++) {
                for (int j = 0; j < C; j++) {
                    if (L[t - 1] == marks.get(i).intValue() && L[t] == marks.get(j).intValue()) {
                        MC[i][j] += Math.pow(T[t], k);
                        count[i][j]++;
                    }
                }
            }
        }

        // Normalize by counts
        for (int i = 0; i < C; i++) {
            for (int j = 0; j < C; j++) {
                if (count[i][j] > 0) {
                    MC[i][j] /= count[i][j];
                } else {
                    MC[i][j] = Double.NaN;
                }
            }
        }

        return new Matrix(MC);
    }
}
