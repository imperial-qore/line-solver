/**
 * @file DTMC stochastic matrix normalization
 *
 * Converts non-negative matrices into valid discrete-time Markov chain transition
 * matrices by row normalization. Ensures row stochasticity property (row sums = 1)
 * required for valid DTMC transition matrices.
 *
 * @since LINE 3.0
 */
package jline.api.mc;

import jline.util.matrix.Matrix;
import org.apache.commons.math3.util.FastMath;

public final class Dtmc_makestochastic {
    private Dtmc_makestochastic() {}

    /**
     * Normalize a given non-negative matrix into a DTMC.
     *
     * @param P nonegative matrix
     * @return Transition matrix of the DTMC
     */
    public static Matrix dtmc_makestochastic(Matrix P) {
        P = P.copy();
        int n = P.length();
        for (int i = 0; i < n; i++) {
            double rowSum = P.getRow(i).elementSum();
            if (rowSum > 0) {
                for (int j = 0; j < n; j++) {
                    P.set(i, j, P.get(i, j) / rowSum);
                }
                rowSum = P.getRow(i).elementSum();
                P.set(i, i, FastMath.min(Math.max(0.0, 1.0 - (rowSum - P.get(i, i))), 1.0));
            } else {
                for (int j = 0; j < n; j++) {
                    P.set(i, j, 0.0);
                }
                P.set(i, i, 1.0);
            }
        }
        return P;
    }
}
