package jline.api.trace;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;

public final class Mtrace_split {
    private Mtrace_split() {}

    /**
     * Given a multi-class trace with inter-arrivals T and labels L,
     * creates the separate per-class traces. In each per-class trace
     * the inter-arrivals are inter-arrivals between events of the same class.
     *
     * @param T vector of inter-arrival times
     * @param L vector of class labels
     * @return array where the c-th element contains the vector of inter-arrival
     *         times for class c
     */
    public static double[][] mtrace_split(double[] T, int[] L) {
        TreeSet<Integer> uniq = new TreeSet<Integer>();
        for (int v : L) uniq.add(v);
        Integer[] labels = uniq.toArray(new Integer[0]);
        int C = labels.length;

        double[][] TL = new double[C][];
        for (int i = 0; i < C; i++) TL[i] = new double[0];

        // Compute cumulative sum of T
        double[] TCUM = new double[T.length + 1];
        TCUM[0] = 0.0;
        for (int i = 0; i < T.length; i++) {
            TCUM[i + 1] = TCUM[i] + T[i];
        }

        for (int c = 0; c < C; c++) {
            int label = labels[c].intValue();

            List<Integer> indices = new ArrayList<Integer>();
            indices.add(0); // Add initial 0
            for (int i = 0; i < L.length; i++) {
                if (L[i] == label) {
                    indices.add(i + 1); // +1 because TCUM is offset by 1
                }
            }

            if (indices.size() > 1) {
                double[] classInterArrivals = new double[indices.size() - 1];
                for (int i = 1; i < indices.size(); i++) {
                    classInterArrivals[i - 1] = TCUM[indices.get(i)] - TCUM[indices.get(i - 1)];
                }
                TL[c] = classInterArrivals;
            }
        }

        return TL;
    }
}
