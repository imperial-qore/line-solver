/**
 * @file Kuiper statistical distance test
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import jline.util.matrix.Matrix;

public final class Ms_kuiper {
    private Ms_kuiper() {}

    /**
     * Kuiper distance between two empirical distributions.
     */
    public static double ms_kuiper(Matrix XX, Matrix YY) {
        if (XX.getNumCols() > 1 || YY.getNumCols() > 1) {
            if (XX.getNumCols() != YY.getNumCols()) {
                throw new IllegalArgumentException("Matrices must have the same number of columns");
            }
            double maxKuiper = 0.0;
            for (int col = 0; col < XX.getNumCols(); col++) {
                Matrix xCol = XX.getColumn(col);
                Matrix yCol = YY.getColumn(col);
                double k = ms_kuiper_1D(xCol, yCol);
                if (k > maxKuiper) maxKuiper = k;
            }
            return maxKuiper;
        }
        return ms_kuiper_1D(XX, YY);
    }

    private static double ms_kuiper_1D(Matrix XX, Matrix YY) {
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
        int n = nx + ny;

        // (value, X_indicator, Y_indicator)
        List<double[]> combined = new ArrayList<double[]>();
        for (Double x : X) {
            combined.add(new double[]{x, 1.0 / nx, 0.0});
        }
        for (Double y : Y) {
            combined.add(new double[]{y, 0.0, 1.0 / ny});
        }

        Collections.sort(combined, new Comparator<double[]>() {
            @Override
            public int compare(double[] a, double[] b) {
                return Double.compare(a[0], b[0]);
            }
        });

        double eCDF = 0.0;
        double fCDF = 0.0;
        double maxUp = 0.0;
        double maxDown = 0.0;

        for (int i = 0; i < n - 1; i++) {
            eCDF += combined.get(i)[1];
            fCDF += combined.get(i)[2];

            if (combined.get(i + 1)[0] != combined.get(i)[0]) {
                double height = fCDF - eCDF;
                if (height > maxUp) maxUp = height;
                if (height < maxDown) maxDown = height;
            }
        }

        return Math.abs(maxDown) + Math.abs(maxUp);
    }
}
