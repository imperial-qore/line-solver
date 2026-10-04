package jline.api.trace;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public final class Mtrace_sigma2 {
    private Mtrace_sigma2() {}

    /**
     * Computes the empirical probability of observing a specific 3-element
     * sequence of events, i.e. the two-step class transition probabilities.
     *
     * @param T the inter-arrival times
     * @param L the event labels
     * @return 3D array where (i,j,h) element is the probability of observing
     *         an event of class i followed by class j followed by class h
     */
    public static double[][][] mtrace_sigma2(double[] T, int[] L) {
        Set<Integer> set = new HashSet<Integer>();
        for (int v : L) set.add(v);
        List<Integer> marks = new ArrayList<Integer>(set);
        Collections.sort(marks);
        int C = marks.size();

        double[][][] sigma = new double[C][C][C];

        for (int i = 0; i < C; i++) {
            for (int j = 0; j < C; j++) {
                for (int h = 0; h < C; h++) {
                    int count = 0;
                    for (int t = 0; t < L.length - 2; t++) {
                        if (L[t] == marks.get(i) && L[t + 1] == marks.get(j) && L[t + 2] == marks.get(h)) {
                            count++;
                        }
                    }
                    sigma[i][j][h] = (double) count / (L.length - 2);
                }
            }
        }

        return sigma;
    }
}
