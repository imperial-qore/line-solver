/**
 * @file Markovian Arrival Process CTMC steady-state analysis
 *
 * Computes steady-state probability vector of underlying continuous-time Markov chain.
 * Essential for MAP stationary analysis and fundamental performance metric calculations.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.api.mc.Ctmc_solve;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_piq {
    private Map_piq() {}

    /**
     * Computes the steady-state vector (pi) of the Continuous-Time Markov Chain (CTMC) underlying a Markovian Arrival Process (MAP).
     *
     * This function calculates the steady-state distribution of the CTMC by solving for the stationary distribution of the infinitesimal
     * generator matrix Q, which is obtained by summing the hidden (D0) and visible (D1) transition matrices.
     *
     * @param D0 The hidden transition matrix of the MAP, representing transitions without visible events.
     * @param D1 The visible transition matrix of the MAP, representing transitions with visible events.
     * @return The steady-state vector (pi) of the CTMC.
     */
    public static Matrix map_piq(Matrix D0, Matrix D1) {
        return Ctmc_solve.ctmc_solve(Map_infgen.map_infgen(D0, D1));
    }

    /**
     * Computes the steady-state vector (pi) of the Continuous-Time Markov Chain (CTMC) underlying a Markovian Arrival Process (MAP).
     *
     * This is a convenience method that extracts the D0 and D1 matrices from a given MAP stored in a MatrixCell and computes the
     * steady-state distribution of the CTMC by solving for the stationary distribution of the infinitesimal generator matrix Q.
     *
     * @param MAP The Markovian Arrival Process stored in a MatrixCell, containing the (D0, D1) matrices.
     * @return The steady-state vector (pi) of the CTMC.
     */
    public static Matrix map_piq(MatrixCell MAP) {
        return map_piq(MAP.get(0), MAP.get(1));
    }
}
