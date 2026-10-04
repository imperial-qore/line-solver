/**
 * @file Markovian Arrival Process hyperexponential distribution fitting
 *
 * Constructs MAP representations of two-phase hyperexponential renewal processes.
 * Used for modeling high-variability arrival processes with coefficient of variation greater than 1.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import org.apache.commons.math3.util.FastMath;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_hyperexp {
    private Map_hyperexp() {}

    /**
     * Fit a two-phase Hyper-exponential renewal process as a MAP.
     *
     * @param mean mean inter-arrival time of the process
     * @param scv  squared coefficient of variation of inter-arrival times
     * @param p    probability of being served in phase 1 (DEFAULT: p=0.99)
     * @return Fitted hyper-exponential process, or null if not feasible
     */
    public static MatrixCell map_hyperexp(double mean, double scv, double p) {
        MatrixCell D = new MatrixCell();
        if (p == 0.0) {
            p = 0.99;
        }

        double e2 = (1.0 + scv) * mean * mean;
        double delta = -4.0 * p * mean * mean + 4.0 * p * p * mean * mean + 2.0 * e2 * p - 2.0 * e2 * p * p;
        double mu2 = (-2.0 * mean + 2.0 * p * mean + FastMath.sqrt(delta)) / (e2 * p - 2.0 * mean * mean);
        double mu1 = mu2 * p / (p - 1.0 + mean * mu2);
        Matrix D0 = new Matrix(2, 2);
        Matrix D1 = new Matrix(2, 2);
        D0.set(0, 0, -mu1);
        D0.set(1, 1, -mu2);
        D1.set(0, 0, mu1 * p);
        D1.set(0, 1, mu1 * (1.0 - p));
        D1.set(1, 0, mu2 * p);
        D1.set(1, 1, mu2 * (1.0 - p));
        D.set(0, D0);
        D.set(1, D1);

        if (Map_isfeasible.map_isfeasible(D)) {
            return D;
        } else {
            mu2 = (-2 * mean + 2 * p * mean - FastMath.sqrt(delta)) / (e2 * p - 2 * mean * mean);
            mu1 = mu2 * p / (p - 1 + mean * mu2);
            D0.zero();
            D1.zero();
            D0.set(0, 0, -mu1);
            D0.set(1, 1, -mu2);
            D1.set(0, 0, mu1 * p);
            D1.set(0, 1, mu1 * (1.0 - p));
            D1.set(1, 0, mu2 * p);
            D1.set(1, 1, mu2 * (1.0 - p));
            D.set(0, D0);
            D.set(1, D1);
            if (Map_isfeasible.map_isfeasible(D)) {
                // the second root is the solution; returning null here would
                // discard a feasible hyperexponential
                return D;
            } else if (p > 1e-6) {
                return map_hyperexp(mean, scv, p / 10.0);
            } else {
                return null;
            }
        }
    }
}
