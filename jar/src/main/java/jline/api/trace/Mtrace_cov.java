package jline.api.trace;

import jline.util.matrix.Matrix;

public final class Mtrace_cov {
    private Mtrace_cov() {}

    /**
     * Computes the covariance matrix for multi-type traces.
     *
     * @param T the trace data (inter-arrival times or values)
     * @param A the type labels for each element in T (1-indexed)
     * @return a matrix array containing covariance values between different types
     */
    public static Matrix[][] mtrace_cov(double[] T, int[] A) {
        int C = 0;
        for (int v : A) {
            if (v > C) C = v;
        }
        int N = A.length;

        Matrix[][] COV = new Matrix[C][C];
        for (int i = 0; i < C; i++) {
            for (int j = 0; j < C; j++) {
                COV[i][j] = Matrix.zeros(2, 2);
            }
        }

        for (int c1 = 1; c1 <= C; c1++) {
            for (int c2 = 1; c2 <= C; c2++) {
                double[] X0c1v = new double[N - 1];
                double[] X1c2v = new double[N - 1];

                for (int i = 0; i < N - 1; i++) {
                    double X0c1 = (A[i] == c1) ? T[i] : 0.0;
                    double X1c2 = (A[i + 1] == c2) ? T[i + 1] : 0.0;
                    X0c1v[i] = X0c1;
                    X1c2v[i] = X1c2;
                }

                Matrix cov = computeCovariance(X0c1v, X1c2v);
                COV[c1 - 1][c2 - 1] = cov;
            }
        }

        return COV;
    }

    /**
     * Computes the 2x2 covariance matrix between two vectors.
     */
    private static Matrix computeCovariance(double[] x, double[] y) {
        int n = x.length;
        double sumX = 0.0;
        double sumY = 0.0;
        for (int i = 0; i < n; i++) {
            sumX += x[i];
            sumY += y[i];
        }
        double meanX = sumX / n;
        double meanY = sumY / n;

        double cov_xx = 0.0;
        double cov_xy = 0.0;
        double cov_yy = 0.0;

        for (int i = 0; i < n; i++) {
            double dx = x[i] - meanX;
            double dy = y[i] - meanY;
            cov_xx += dx * dx;
            cov_xy += dx * dy;
            cov_yy += dy * dy;
        }

        cov_xx /= (n - 1);
        cov_xy /= (n - 1);
        cov_yy /= (n - 1);

        return new Matrix(new double[][]{
                {cov_xx, cov_xy},
                {cov_xy, cov_yy}
        });
    }
}
