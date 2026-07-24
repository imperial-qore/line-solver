package jline.api.trace;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.TreeSet;

import jline.util.matrix.Matrix;

public final class Mtrace_pc {
    private Mtrace_pc() {}

    /**
     * Computes the probabilities of arrival for each class.
     *
     * @param T the inter-arrival times (ignored, for orthogonality with other APIs)
     * @param C the class labels
     * @return column vector of probabilities for each class
     */
    public static Matrix mtrace_pc(double[] T, int[] C) {
        TreeSet<Integer> distinct = new TreeSet<Integer>();
        for (int c : C) {
            distinct.add(c);
        }
        List<Integer> labels = new ArrayList<Integer>(distinct);
        Collections.sort(labels);
        int m = labels.size();

        double[] pc = new double[m];

        for (int i = 0; i < m; i++) {
            int label = labels.get(i);
            int count = 0;
            for (int c : C) {
                if (c == label) {
                    count++;
                }
            }
            pc[i] = (double) count / C.length;
        }

        return new Matrix(pc);
    }
}
