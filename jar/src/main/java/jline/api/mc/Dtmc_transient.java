/**
 * @file Transient distribution of a discrete-time Markov chain
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import jline.util.matrix.Matrix;

public final class Dtmc_transient {
    private Dtmc_transient() {}

    /**
     * Transient distribution of a DTMC, i.e. the rows pi(k) = pi0*P^k for
     * k = 0,...,steps. Twin of the MATLAB dtmc_transient and of the Python
     * api.mc.dtmc_transient.
     *
     * @param P     transition matrix
     * @param pi0   initial distribution, uniform when null
     * @param steps number of steps
     * @return matrix with one row per step, from step 0 to step steps
     */
    public static Matrix dtmc_transient(Matrix P, Matrix pi0, int steps) {
        int n = P.getNumRows();
        if (steps < 0) {
            throw new IllegalArgumentException("The number of steps must be non-negative.");
        }
        Matrix pik = new Matrix(1, n);
        for (int i = 0; i < n; i++) {
            pik.set(0, i, (pi0 == null || pi0.length() != n) ? 1.0 / n : pi0.get(i));
        }
        Matrix pit = new Matrix(steps + 1, n);
        for (int k = 0; k <= steps; k++) {
            for (int i = 0; i < n; i++) {
                pit.set(k, i, pik.get(0, i));
            }
            pik = pik.mult(P);
        }
        return pit;
    }
}
