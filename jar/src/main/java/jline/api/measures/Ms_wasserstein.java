/**
 * @file Wasserstein (Earth Mover's) distance
 *
 * @since LINE 3.0
 */
package jline.api.measures;

import jline.util.matrix.Matrix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public final class Ms_wasserstein {
    private Ms_wasserstein() {}

    /**
     * Wasserstein distance (Earth Mover's Distance) between two empirical distributions.
     */
    public static double ms_wasserstein(Matrix XX, Matrix YY) {
        if (XX.getNumCols() > 1 || YY.getNumCols() > 1) {
            if (XX.getNumCols() != YY.getNumCols()) {
                throw new IllegalArgumentException("Matrices must have the same number of columns");
            }
            double maxWS = 0.0;
            for (int col = 0; col < XX.getNumCols(); col++) {
                Matrix xCol = XX.getColumn(col);
                Matrix yCol = YY.getColumn(col);
                double ws = ms_wasserstein_1D(xCol, yCol);
                if (ws > maxWS) maxWS = ws;
            }
            return maxWS;
        }
        return ms_wasserstein_1D(XX, YY);
    }

    private static final class Pair {
        final double value;
        final double weight;
        Pair(double v, double w) { this.value = v; this.weight = w; }
    }

    private static double ms_wasserstein_1D(Matrix XX, Matrix YY) {
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

        List<Pair> combined = new ArrayList<Pair>();
        for (int i = 0; i < nx; i++) {
            combined.add(new Pair(X.get(i), 1.0 / nx));
        }
        for (int i = 0; i < ny; i++) {
            combined.add(new Pair(Y.get(i), -1.0 / ny));
        }

        Collections.sort(combined, new Comparator<Pair>() {
            @Override
            public int compare(Pair a, Pair b) {
                return Double.compare(a.value, b.value);
            }
        });

        double distance = 0.0;
        double cumulativeWeight = 0.0;
        for (int i = 0; i < combined.size() - 1; i++) {
            cumulativeWeight += combined.get(i).weight;
            double width = combined.get(i + 1).value - combined.get(i).value;
            distance += Math.abs(cumulativeWeight) * width;
        }
        return distance;
    }
}
