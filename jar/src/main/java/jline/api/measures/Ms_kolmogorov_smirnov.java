/**
 * @file Kolmogorov-Smirnov goodness-of-fit test statistic
 *
 * Implements the Kolmogorov-Smirnov test for determining if a sample follows
 * a hypothesized distribution. Measures the maximum difference between empirical
 * and theoretical cumulative distribution functions for distribution comparison.
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import jline.util.matrix.Matrix;

public final class Ms_kolmogorov_smirnov {
    private Ms_kolmogorov_smirnov() {}

    /**
     * Kolmogorov-Smirnov distance between two empirical distributions.
     */
    public static double ms_kolmogorov_smirnov(Matrix XX, Matrix YY) {
        if (XX.getNumCols() > 1 || YY.getNumCols() > 1) {
            if (XX.getNumCols() != YY.getNumCols()) {
                throw new IllegalArgumentException("Matrices must have the same number of columns");
            }

            double maxKS = 0.0;
            for (int col = 0; col < XX.getNumCols(); col++) {
                Matrix xCol = XX.getColumn(col);
                Matrix yCol = YY.getColumn(col);
                double ks = ms_kolmogorov_smirnov_1D(xCol, yCol);
                maxKS = Math.max(maxKS, ks);
            }
            return maxKS;
        }

        return ms_kolmogorov_smirnov_1D(XX, YY);
    }

    private static double ms_kolmogorov_smirnov_1D(Matrix XX, Matrix YY) {
        List<Double> X = new ArrayList<Double>();
        List<Double> Y = new ArrayList<Double>();

        for (int i = 0; i < XX.length(); i++) {
            double xi = XX.get(i);
            if (!Double.isNaN(xi)) X.add(xi);
        }

        for (int i = 0; i < YY.length(); i++) {
            double yi = YY.get(i);
            if (!Double.isNaN(yi)) Y.add(yi);
        }

        if (X.isEmpty() || Y.isEmpty()) return 0.0;

        int nx = X.size();
        int ny = Y.size();

        List<Double> xSorted = new ArrayList<Double>(X);
        Collections.sort(xSorted);
        List<Double> ySorted = new ArrayList<Double>(Y);
        Collections.sort(ySorted);

        List<Double> combined = new ArrayList<Double>(X);
        combined.addAll(Y);
        Collections.sort(combined);

        double eCDF = 0.0;
        double fCDF = 0.0;
        double maxDist = 0.0;

        int xi = 0;
        int yi = 0;

        for (double value : combined) {
            while (xi < nx && xSorted.get(xi) <= value) {
                eCDF += 1.0 / nx;
                xi++;
            }

            while (yi < ny && ySorted.get(yi) <= value) {
                fCDF += 1.0 / ny;
                yi++;
            }

            double dist = Math.abs(eCDF - fCDF);
            maxDist = Math.max(maxDist, dist);
        }

        return maxDist;
    }
}
