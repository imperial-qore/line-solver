package jline.api.measures;

import jline.util.matrix.Matrix;
import jline.util.Pair;

import java.util.HashMap;
import java.util.Map;

public final class Ms_jointentropy {
    private Ms_jointentropy() {}

    /**
     * Compute joint entropy z=H(x,y) of two discrete variables x and y.
     *
     * @param x first matrix
     * @param y second matrix of the same length
     * @return joint entropy z=H(x,y)
     */
    public static double ms_jointentropy(Matrix x, Matrix y) {
        if (x.length() != y.length()) {
            throw new IllegalArgumentException("Input matrices must have the same length");
        }

        if (x.isEmpty()) {
            return 0.0;
        }

        // Convert to integer values and find min
        int len = x.length();
        int[] xVals = new int[len];
        int[] yVals = new int[len];
        int minVal = Integer.MAX_VALUE;

        for (int i = 0; i < len; i++) {
            int xi = (int) x.get(i);
            int yi = (int) y.get(i);
            xVals[i] = xi;
            yVals[i] = yi;
            if (xi < minVal) minVal = xi;
            if (yi < minVal) minVal = yi;
        }

        // Shift values to start from 1 (like in MATLAB)
        // Count joint occurrences
        Map<Pair<Integer, Integer>, Integer> jointCounts = new HashMap<Pair<Integer, Integer>, Integer>();
        for (int i = 0; i < len; i++) {
            Pair<Integer, Integer> pair = new Pair<Integer, Integer>(xVals[i] - minVal + 1, yVals[i] - minVal + 1);
            Integer prev = jointCounts.get(pair);
            jointCounts.put(pair, prev == null ? 1 : prev + 1);
        }

        double n = (double) len;
        double jointEntropy = 0.0;

        // Calculate joint entropy
        double log2 = Math.log(2.0);
        for (Integer count : jointCounts.values()) {
            double p = count / n;
            if (p > 0) {
                jointEntropy -= p * (Math.log(p) / log2);
            }
        }

        return Math.max(0.0, jointEntropy);
    }
}
