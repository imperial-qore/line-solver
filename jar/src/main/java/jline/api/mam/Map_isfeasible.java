/**
 * @file Markovian Arrival Process feasibility checking interface
 *
 * Provides convenient interface for MAP feasibility validation with configurable tolerance.
 * Wrapper around detailed feasibility checking algorithms for ease of use.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.MatrixCell;
import org.apache.commons.math3.util.FastMath;

public final class Map_isfeasible {
    private Map_isfeasible() {}

    /**
     * Checks if the provided MAP is feasible based on the given tolerance.
     *
     * <p>This method evaluates whether the MAP transition matrices meet the necessary conditions for a
     * valid MAP, such as non-negativity of elements, proper row/column sums, etc.
     *
     * @param MAP the MatrixCell representing the MAP transition matrices
     * @param TOL the tolerance level for numerical stability checks
     * @return true if the MAP is feasible, false otherwise
     */
    public static boolean map_isfeasible(MatrixCell MAP, double TOL) {
        return Map_checkfeasible.map_checkfeasible(MAP, TOL);
    }

    /**
     * Checks if the provided MAP is feasible using a default tolerance.
     *
     * <p>This method uses a standard tolerance level to determine the feasibility of the MAP.
     *
     * @param MAP the MatrixCell representing the MAP transition matrices
     * @return true if the MAP is feasible, false otherwise
     */
    public static boolean map_isfeasible(MatrixCell MAP) {
        int TOLMAGNITUDE = 15;
        for (int k = TOLMAGNITUDE; k >= 1; k--) {
            boolean check = Map_checkfeasible.map_checkfeasible(MAP, FastMath.pow(10.0, -k));
            if (check) {
                return k > Map_feastol.map_feastol();
            }
        }
        return false;
    }
}
