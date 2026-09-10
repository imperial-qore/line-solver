package jline.api.trace;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import jline.util.Pair;

public final class Mtrace_merge {
    private Mtrace_merge() {}

    /**
     * Merges two traces in a single marked (multiclass) trace.
     *
     * @param t1 inter-arrival times of the first trace
     * @param t2 inter-arrival times of the second trace
     * @return pair of (inter-arrival times of marked process, labels of marked process)
     */
    public static Pair<double[], int[]> mtrace_merge(double[] t1, double[] t2) {
        // Create cumulative sums starting from 0
        double[] cum1 = new double[t1.length + 1];
        double[] cum2 = new double[t2.length + 1];

        cum1[0] = 0.0;
        for (int i = 0; i < t1.length; i++) {
            cum1[i + 1] = cum1[i] + t1[i];
        }

        cum2[0] = 0.0;
        for (int i = 0; i < t2.length; i++) {
            cum2[i + 1] = cum2[i] + t2[i];
        }

        // see _kb/03-api-layer.md for rationale
        List<Pair<Double, Integer>> allTimes = new ArrayList<Pair<Double, Integer>>();
        allTimes.add(new Pair<Double, Integer>(0.0, 0));

        // Events of the first trace (indices 1 to t1.length)
        for (int i = 1; i < cum1.length; i++) {
            allTimes.add(new Pair<Double, Integer>(cum1[i], i));
        }

        // Events of the second trace (indices t1.length+2 to t1.length+t2.length+1)
        for (int i = 1; i < cum2.length; i++) {
            allTimes.add(new Pair<Double, Integer>(cum2[i], t1.length + 1 + i));
        }

        // Sort by time
        Collections.sort(allTimes, new Comparator<Pair<Double, Integer>>() {
            @Override
            public int compare(Pair<Double, Integer> a, Pair<Double, Integer> b) {
                return Double.compare(a.getLeft(), b.getLeft());
            }
        });

        // Compute inter-arrival times
        double[] T = new double[allTimes.size() - 1];
        for (int i = 1; i < allTimes.size(); i++) {
            T[i - 1] = allTimes.get(i).getLeft() - allTimes.get(i - 1).getLeft();
        }

        // Determine labels
        int[] L = new int[T.length];
        for (int i = 1; i < allTimes.size(); i++) {
            int idx = allTimes.get(i).getRight();
            if (idx >= 1 && idx <= t1.length) {
                L[i - 1] = 1;
            } else if (idx >= t1.length + 2 && idx <= t1.length + 1 + t2.length) {
                L[i - 1] = 2;
            } else {
                L[i - 1] = 0;
            }
        }

        return new Pair<double[], int[]>(T, L);
    }
}
