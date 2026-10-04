package jline.api.mc;

import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

public final class Dtmc_isfeasible {
    private Dtmc_isfeasible() {}

    /**
     * Check if a matrix represents a feasible DTMC transition matrix
     *
     * @param P Transition matrix to check
     * @return Feasibility tolerance level (0 if not feasible, higher values indicate better precision)
     */
    public static int dtmc_isfeasible(Matrix P) {
        int n = P.getNumRows();
        double[] sP = new double[n];

        // Compute row sums
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                sum += P.get(i, j);
            }
            sP[i] = sum;
        }

        double minSum = sP[0];
        double maxSum = sP[0];
        for (int i = 1; i < n; i++) {
            minSum = FastMath.min(minSum, sP[i]);
            maxSum = FastMath.max(maxSum, sP[i]);
        }

        double minElement = P.get(0, 0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                minElement = FastMath.min(minElement, P.get(i, j));
            }
        }

        int res = 0;
        for (int tol = 1; tol <= 15; tol++) {
            double tolerance = FastMath.pow(10.0, -(double) tol);
            if (minSum > 1.0 - tolerance &&
                    maxSum < 1.0 + tolerance &&
                    minElement > -tolerance) {
                res = tol;
            }
        }

        return res;
    }
}
