/**
 * @file Relative entropy (Kullback-Leibler divergence) for discrete variables
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import java.util.ArrayList;
import java.util.List;

import jline.GlobalConstants;
import jline.util.matrix.Matrix;

public final class Ms_relatentropy {
    private Ms_relatentropy() {}

    /**
     * Compute relative entropy (KL divergence) z=KL(p(x)||p(y)) of two discrete variables x and y.
     */
    public static double ms_relatentropy(Matrix x, Matrix y) {
        if (x.length() != y.length()) {
            throw new IllegalArgumentException("Input matrices must have the same length");
        }

        if (x.isEmpty()) return 0.0;

        List<Integer> xVals = new ArrayList<Integer>();
        List<Integer> yVals = new ArrayList<Integer>();
        int minVal = Integer.MAX_VALUE;
        int maxVal = Integer.MIN_VALUE;

        for (int i = 0; i < x.length(); i++) {
            int xi = (int) x.get(i);
            int yi = (int) y.get(i);
            xVals.add(xi);
            yVals.add(yi);
            minVal = Math.min(minVal, Math.min(xi, yi));
            maxVal = Math.max(maxVal, Math.max(xi, yi));
        }

        int range = maxVal - minVal + 1;

        int[] countsX = new int[range];
        int[] countsY = new int[range];

        for (Integer value : xVals) {
            countsX[value - minVal]++;
        }
        for (Integer value : yVals) {
            countsY[value - minVal]++;
        }

        double n = (double) x.length();
        double klDivergence = 0.0;

        double inv_log2 = 1.0 / Math.log(2.0);
        for (int i = 0; i < range; i++) {
            double px = countsX[i] / n;
            double py = countsY[i] / n;

            if (px > 0 && py > 0) {
                klDivergence += px * (Math.log(px) * inv_log2 - Math.log(py) * inv_log2);
            } else if (px > 0 && py == 0.0) {
                return GlobalConstants.Inf;
            }
        }

        return Math.max(0.0, klDivergence);
    }
}
