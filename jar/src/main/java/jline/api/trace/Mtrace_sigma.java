package jline.api.trace;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import jline.util.matrix.Matrix;

public final class Mtrace_sigma {
    private Mtrace_sigma() {}

    /**
     * Computes the empirical probability of observing a specific 2-element
     * sequence of events, i.e. the one-step class transition probabilities.
     *
     * @param T the inter-arrival times
     * @param L the event labels
     * @return matrix where (i,j) element is the probability of observing
     *         an event of class i followed by an event of class j
     */
    public static Matrix mtrace_sigma(double[] T, int[] L) {
        Set<Integer> set = new HashSet<Integer>();
        for (int v : L) set.add(v);
        List<Integer> marks = new ArrayList<Integer>(set);
        Collections.sort(marks);
        int C = marks.size();

        double[][] sigma = new double[C][C];

        for (int i = 0; i < C; i++) {
            for (int j = 0; j < C; j++) {
                int count = 0;
                for (int t = 0; t < L.length - 1; t++) {
                    if (L[t] == marks.get(i) && L[t + 1] == marks.get(j)) {
                        count++;
                    }
                }
                sigma[i][j] = (double) count / (L.length - 1);
            }
        }

        return new Matrix(sigma);
    }
}
