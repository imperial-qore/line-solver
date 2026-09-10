/**
 * @file Feasibility check for a continuous-time Markov chain generator
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import jline.util.matrix.Matrix;

public final class Ctmc_isfeasible {
    private Ctmc_isfeasible() {}

    /**
     * Checks that Q is a valid infinitesimal generator, with the default
     * tolerance of 1e-10.
     *
     * @param Q candidate generator matrix
     * @return true when Q is a valid generator
     */
    public static boolean ctmc_isfeasible(Matrix Q) {
        return ctmc_isfeasible(Q, 1e-10);
    }

    /**
     * Checks that Q is a valid infinitesimal generator: square, non-negative
     * off-diagonal entries, non-positive diagonal, and zero row sums, each up
     * to the given tolerance. Twin of the MATLAB ctmc_isfeasible and of the
     * Python api.mc.ctmc_isfeasible; note that dtmc_isfeasible instead returns
     * a precision level rather than a flag.
     *
     * @param Q         candidate generator matrix
     * @param tolerance numerical tolerance
     * @return true when Q is a valid generator
     */
    public static boolean ctmc_isfeasible(Matrix Q, double tolerance) {
        int n = Q.getNumRows();
        if (n == 0 || n != Q.getNumCols()) {
            return false;
        }
        for (int i = 0; i < n; i++) {
            double rowSum = 0;
            for (int j = 0; j < n; j++) {
                double q = Q.get(i, j);
                if (i != j && q < -tolerance) {
                    return false;
                }
                if (i == j && q > tolerance) {
                    return false;
                }
                rowSum += q;
            }
            if (Math.abs(rowSum) > tolerance) {
                return false;
            }
        }
        return true;
    }
}
