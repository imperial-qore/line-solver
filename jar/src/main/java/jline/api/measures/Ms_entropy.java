package jline.api.measures;

import java.util.HashMap;
import java.util.Map;

import jline.util.matrix.Matrix;

public final class Ms_entropy {
    private Ms_entropy() {}

    /**
     * Compute entropy z=H(x) of a discrete variable x.
     *
     * @param x a matrix representing discrete values
     * @return entropy z=H(x)
     */
    public static double ms_entropy(Matrix x) {
        if (x.isEmpty()) return 0.0;

        // Count occurrences of each unique value
        Map<Integer, Integer> counts = new HashMap<Integer, Integer>();
        for (int i = 0; i < x.length(); i++) {
            int value = (int) x.get(i);
            Integer prev = counts.get(value);
            counts.put(value, (prev == null ? 0 : prev) + 1);
        }

        double n = (double) x.length();
        double entropy = 0.0;

        // Calculate entropy
        for (int count : counts.values()) {
            double p = count / n;
            if (p > 0) {
                entropy -= p * (Math.log(p) / Math.log(2.0));
            }
        }

        return Math.max(0.0, entropy);
    }
}
