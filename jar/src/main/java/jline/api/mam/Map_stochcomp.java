/**
 * @file Markovian Arrival Process stochastic complementation
 *
 * Performs state elimination through stochastic complementation while preserving MAP properties.
 * Used for model reduction and aggregation in large-scale MAP analysis.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Map_stochcomp {
    private Map_stochcomp() {}

    /**
     * Performs stochastic complementation on a MAP by eliminating specified states.
     */
    public static MatrixCell map_stochcomp(MatrixCell MAP, int[] retainIdx) {
        return map_stochcomp(MAP.get(0), MAP.get(1), retainIdx);
    }

    /**
     * Performs stochastic complementation on a MAP by eliminating specified states.
     */
    public static MatrixCell map_stochcomp(Matrix D0, Matrix D1, int[] retainIdx) {
        int n = D0.getNumRows();

        // Validate retain indices
        for (int idx : retainIdx) {
            if (idx < 0 || idx >= n) {
                throw new IllegalArgumentException("Retain index " + idx + " is out of bounds [0, " + (n - 1) + "]");
            }
        }

        // Create the full Q matrix
        Matrix Q = D0.add(1.0, D1);

        // Find eliminated indices
        java.util.List<Integer> elim = new java.util.ArrayList<Integer>();
        for (int i = 0; i < n; i++) {
            boolean keep = false;
            for (int r : retainIdx) {
                if (r == i) { keep = true; break; }
            }
            if (!keep) elim.add(i);
        }
        int[] eliminatedIdx = new int[elim.size()];
        for (int i = 0; i < elim.size(); i++) eliminatedIdx[i] = elim.get(i);

        if (eliminatedIdx.length == 0) {
            // Nothing to eliminate, return original MAP
            return Map_normalize.map_normalize(D0, D1);
        }

        // Extract submatrices
        Matrix Q_RE = extractSubmatrix(Q, retainIdx, eliminatedIdx);
        Matrix Q_EE = extractSubmatrix(Q, eliminatedIdx, eliminatedIdx);
        Matrix Q_RR = extractSubmatrix(Q, retainIdx, retainIdx);
        Matrix Q_ER = extractSubmatrix(Q, eliminatedIdx, retainIdx);

        // Compute the new Q matrix using stochastic complementation
        // QNew = Q_RR + Q_RE * (-Q_EE)^(-1) * Q_ER
        Matrix minusQ_EE = Q_EE.scale(-1.0);
        Matrix invMinusQ_EE = minusQ_EE.inv();
        Matrix QNew = Q_RR.add(1.0, Q_RE.mult(invMinusQ_EE).mult(Q_ER));

        // Extract D1 submatrices
        Matrix D1_RR = extractSubmatrix(D1, retainIdx, retainIdx);
        Matrix D1_ER = extractSubmatrix(D1, eliminatedIdx, retainIdx);

        // Compute new D0 and D1 matrices
        Matrix D0new = QNew.sub(1.0, D1_RR);
        Matrix D1new = D1_RR.add(1.0, Q_RE.mult(invMinusQ_EE).mult(D1_ER));

        return Map_normalize.map_normalize(D0new, D1new);
    }

    /**
     * Extracts a submatrix from the given matrix based on row and column indices.
     */
    private static Matrix extractSubmatrix(Matrix matrix, int[] rowIndices, int[] colIndices) {
        Matrix result = new Matrix(rowIndices.length, colIndices.length);
        for (int i = 0; i < rowIndices.length; i++) {
            for (int j = 0; j < colIndices.length; j++) {
                result.set(i, j, matrix.get(rowIndices[i], colIndices[j]));
            }
        }
        return result;
    }
}
