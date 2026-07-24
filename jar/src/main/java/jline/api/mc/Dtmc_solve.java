/**
 * @file Discrete-time Markov chain steady-state solver
 *
 * Computes the steady-state probability distribution for DTMCs by converting the
 * transition matrix problem (P-I)x = 0 into a CTMC-equivalent system and leveraging
 * the robust CTMC solver with automatic reducibility handling.
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import jline.util.matrix.Matrix;

public final class Dtmc_solve {
    private Dtmc_solve() {}

    /**
     * Returns the steady-state solution of a DTMC.
     *
     * @param P Transition matrix of the DTMC
     * @return Steady-state solution vector of the DTMC
     */
    public static Matrix dtmc_solve(Matrix P) {
        Matrix Plocal = P.copy();
        // P - eye(size(P))
        for (int i = 0; i < Plocal.getNumRows(); i++) {
            Plocal.set(i, i, Plocal.get(i, i) - 1.0);
        }
        return Ctmc_solve.ctmc_solve(Plocal);
    }
}
