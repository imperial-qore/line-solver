/**
 * @file CTMC time-reversal transformation
 *
 * Computes the infinitesimal generator of the time-reversed continuous-time
 * Markov chain using detailed balance equations. Time-reversal is fundamental
 * in queueing theory and statistical mechanics for analyzing equilibrium properties.
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import jline.util.matrix.Matrix;

public final class Ctmc_timereverse {
    private Ctmc_timereverse() {}

    /**
     * Compute the infinitesimal generator of the time-reserved CTMC
     *
     * @param Q Infinitesimal generator of the CTMC
     * @return Infinitesimal generator of the time-reversed CTMC
     */
    public static Matrix ctmc_timereverse(Matrix Q) {
        Matrix piq = Ctmc_solve.ctmc_solve(Q);
        Matrix Qrev = new Matrix(Q.getNumCols(), Q.getNumRows());
        for (int i = 0; i < Q.getNumRows(); i++) {
            for (int j = 0; j < Q.getNumCols(); j++) {
                Qrev.set(i, j, Q.get(i, j) * piq.get(i) / piq.get(j));
            }
        }
        return Qrev.transpose();
    }
}
