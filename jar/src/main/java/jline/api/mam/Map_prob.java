/**
 * @file Markovian Arrival Process equilibrium distribution computation
 *
 * Computes equilibrium probability distribution of underlying CTMC for MAP analysis.
 * Essential for steady-state analysis and fundamental MAP performance calculations.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.api.mc.Ctmc_solve;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_prob {
    private Map_prob() {}

    /**
     * Computes the equilibrium distribution of the underlying continuous-time Markov chain for a MAP.
     *
     * This function calculates the steady-state probabilities of the continuous-time Markov chain
     * underlying a Markovian Arrival Process (MAP). The equilibrium distribution is computed
     * by solving the system Q*pi = 0, where Q = D0 + D1 is the infinitesimal generator matrix.
     *
     * @param MAP The Markovian Arrival Process stored in a MatrixCell, containing the (D0, D1) matrices
     * @return The equilibrium distribution as a Matrix (row vector)
     */
    public static Matrix map_prob(MatrixCell MAP) {
        return Ctmc_solve.ctmc_solve(Map_infgen.map_infgen(MAP));
    }

    /**
     * Computes the equilibrium distribution of the underlying continuous-time Markov chain for a MAP.
     *
     * This function calculates the steady-state probabilities of the continuous-time Markov chain
     * underlying a Markovian Arrival Process (MAP) given the D0 and D1 matrices directly.
     *
     * @param D0 The hidden transition matrix of the MAP
     * @param D1 The visible transition matrix of the MAP
     * @return The equilibrium distribution as a Matrix (row vector)
     */
    public static Matrix map_prob(Matrix D0, Matrix D1) {
        return Ctmc_solve.ctmc_solve(Map_infgen.map_infgen(D0, D1));
    }
}
