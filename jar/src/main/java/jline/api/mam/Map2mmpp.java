/**
 * @file Markovian Arrival Process to Markov Modulated Poisson Process conversion
 *
 * Converts MAP representations to MMPP format for compatibility with MMPP-specific algorithms.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.Pair;
import jline.util.matrix.Matrix;

import static jline.io.InputOutput.line_warning;

public final class Map2mmpp {
    private Map2mmpp() {}

    /**
     * Convert a MAP to MMPP format by extracting generator matrix Q and rate matrix LAMBDA.
     *
     * @param MAP Input MAP as Matrix[] where MAP[0] = D0 and MAP[1] = D1
     * @return Pair containing (Q, LAMBDA) where Q = D0 + D1 and LAMBDA = D1
     */
    public static Pair<Matrix, Matrix> map2mmpp(Matrix[] MAP) {
        if (MAP.length != 2) {
            throw new IllegalArgumentException("MAP must contain exactly 2 matrices [D0, D1]");
        }

        Matrix D0 = MAP[0];
        Matrix D1 = MAP[1];

        // Check if D1 is diagonal (true MMPP requirement)
        double[] D1RowVec = new double[D1.getNumRows()];
        for (int i = 0; i < D1.getNumRows(); i++) {
            D1RowVec[i] = D1.get(i, i);
        }
        Matrix D1_diag_matrix = Matrix.diag(D1RowVec);
        Matrix diff = D1.sub(D1_diag_matrix);

        double tolerance = 1e-10;
        double maxDiff = 0.0;
        for (int i = 0; i < diff.getNumRows(); i++) {
            for (int j = 0; j < diff.getNumCols(); j++) {
                double absVal = Math.abs(diff.get(i, j));
                if (absVal > maxDiff) {
                    maxDiff = absVal;
                }
            }
        }

        if (maxDiff > tolerance) {
            line_warning("map2mmpp", "The MAP is not a MMPP, LAMBDA is not diagonal");
        }

        Matrix Q = D0.add(D1);
        Matrix LAMBDA = D1.copy();

        return new Pair<Matrix, Matrix>(Q, LAMBDA);
    }
}
