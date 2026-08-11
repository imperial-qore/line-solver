package jline.lib.perm;

import jline.util.Pair;
import jline.util.matrix.Matrix;

/**
 * Top-level utilities for the queueing-network permanent computations.
 */
public final class QueueingNetwork {
    private QueueingNetwork() {}

    /**
     * Make a matrix doubly stochastic via the Sinkhorn algorithm.
     *
     * @return Pair of (doubly stochastic matrix, rescaling factor)
     */
    public static Pair<Matrix, Double> preprocessingDS(Matrix M) {
        int n = M.getNumRows();
        double minValue = 2.220446049250314e-16;

        double[][] data = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                data[i][j] = Math.max(M.get(i, j), minValue);
            }
        }
        Matrix result = new Matrix(data);

        Matrix X = Matrix.eye(n);
        Matrix Y = Matrix.eye(n);

        double maxRowError = Double.POSITIVE_INFINITY;
        double maxColError = Double.POSITIVE_INFINITY;
        double tolerance = 0.001;

        while (maxRowError > tolerance || maxColError > tolerance) {
            double[] colSums = new double[n];
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int i = 0; i < n; i++) s += result.get(i, j);
                colSums[j] = s;
            }
            for (int j = 0; j < n; j++) {
                if (colSums[j] > 0) {
                    for (int i = 0; i < n; i++) {
                        result.set(i, j, result.get(i, j) / colSums[j]);
                        Y.set(i, j, Y.get(i, j) / colSums[j]);
                    }
                }
            }

            double[] rowSums = new double[n];
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                for (int j = 0; j < n; j++) s += result.get(i, j);
                rowSums[i] = s;
            }
            for (int i = 0; i < n; i++) {
                if (rowSums[i] > 0) {
                    for (int j = 0; j < n; j++) {
                        result.set(i, j, result.get(i, j) / rowSums[i]);
                        X.set(i, j, X.get(i, j) / rowSums[i]);
                    }
                }
            }

            double newMaxColError = 0.0;
            double newMaxRowError = 0.0;
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int i = 0; i < n; i++) s += result.get(i, j);
                newMaxColError = Math.max(newMaxColError, Math.abs(s - 1.0));
            }
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                for (int j = 0; j < n; j++) s += result.get(i, j);
                newMaxRowError = Math.max(newMaxRowError, Math.abs(s - 1.0));
            }
            maxColError = newMaxColError;
            maxRowError = newMaxRowError;
        }

        double rescalingFactor = 1.0;
        for (int i = 0; i < n; i++) {
            rescalingFactor *= X.get(i, i) * Y.get(i, i);
        }

        return new Pair<Matrix, Double>(result, rescalingFactor);
    }
}
