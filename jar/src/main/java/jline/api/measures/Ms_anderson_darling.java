package jline.api.measures;

import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public final class Ms_anderson_darling {
    private Ms_anderson_darling() {}

    /**
     * Anderson-Darling distance between two empirical distributions.
     */
    public static double ms_anderson_darling(Matrix XX, Matrix YY) {
        if (XX.getNumCols() > 1 || YY.getNumCols() > 1) {
            if (XX.getNumCols() != YY.getNumCols()) {
                throw new IllegalArgumentException("Matrices must have the same number of columns");
            }
            double maxAD = 0.0;
            for (int col = 0; col < XX.getNumCols(); col++) {
                Matrix xCol = XX.getColumn(col);
                Matrix yCol = YY.getColumn(col);
                double ad = ms_anderson_darling_1D(xCol, yCol);
                if (ad > maxAD) maxAD = ad;
            }
            return maxAD;
        }
        return ms_anderson_darling_1D(XX, YY);
    }

    private static final class Triple {
        final double value;
        final double xInd;
        final double yInd;
        Triple(double v, double x, double y) { this.value = v; this.xInd = x; this.yInd = y; }
    }

    private static double ms_anderson_darling_1D(Matrix XX, Matrix YY) {
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

        List<Triple> combined = new ArrayList<Triple>(n);
        for (int i = 0; i < nx; i++) {
            combined.add(new Triple(X.get(i), 1.0 / nx, 0.0));
        }
        for (int i = 0; i < ny; i++) {
            combined.add(new Triple(Y.get(i), 0.0, 1.0 / ny));
        }

        Collections.sort(combined, new Comparator<Triple>() {
            @Override
            public int compare(Triple a, Triple b) {
                return Double.compare(a.value, b.value);
            }
        });

        double result = 0.0;
        double eCDF = 0.0;
        double fCDF = 0.0;
        double gCDF = 0.0;
        double power = 2.0;

        for (int i = 0; i < n - 1; i++) {
            eCDF += combined.get(i).xInd;
            fCDF += combined.get(i).yInd;
            gCDF += 1.0 / n;

            double sd = Math.sqrt(n * gCDF * Math.abs(1.0 - gCDF));
            double height = Math.abs(fCDF - eCDF);

            if (combined.get(i).value != combined.get(i + 1).value) {
                if (sd > 0) {
                    result += Math.pow(height / sd, power);
                }
            }
        }

        return result;
    }
}
