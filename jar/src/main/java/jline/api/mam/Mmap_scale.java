/**
 * @file Marked Markovian Arrival Process temporal scaling operations
 *
 * Rescales MMAP inter-arrival distributions to achieve specified mean values.
 * Essential for model calibration and parameter adjustment in multiclass systems.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_scale {
    private Mmap_scale() {}

    /**
     * Changes the mean inter-arrival time of a Markovian Arrival Process with marked arrivals (MMAP).
     */
    public static MatrixCell mmap_scale(MatrixCell MMAP, Matrix M, int maxIter) {
        int C = MMAP.size() - 2;
        MatrixCell SCALED = new MatrixCell(2 + C);

        if (M.length() == 1) {
            // Single scaling factor case - uniform scaling
            double MOLD = Map_mean.map_mean(MMAP.get(0), MMAP.get(1));
            double ratio = MOLD / M.get(0);

            Matrix D0 = new Matrix(MMAP.get(0));
            D0.scaleEq(ratio);
            SCALED.set(0, D0);
            Matrix D1 = new Matrix(MMAP.get(1));
            D1.scaleEq(ratio);
            SCALED.set(1, D1);

            for (int c = 0; c < C; c++) {
                Matrix a = new Matrix(MMAP.get(2 + c));
                a.scaleEq(ratio);
                SCALED.set(2 + c, a);
            }
        } else {
            // Multiple scaling factors case - requires iterative approximation
            SCALED.set(0, MMAP.get(0).copy());
            SCALED.set(1, new Matrix(MMAP.get(0).getNumRows(), MMAP.get(0).getNumCols(),
                    MMAP.get(0).getNumRows() * MMAP.get(0).getNumCols()));
            Matrix l = Mmap_count_lambda.mmap_count_lambda(MMAP);

            // Initial heuristic approximation
            for (int c = 0; c < C; c++) {
                if (l.get(c) > 0) {
                    Matrix a = MMAP.get(2 + c).copy();
                    a.scaleEq((1 / M.get(c)) / l.get(c));
                    SCALED.set(2 + c, a);
                    SCALED.set(1, SCALED.get(1).add(1.0, a));
                } else {
                    SCALED.set(2 + c, new Matrix(MMAP.get(2 + c).getNumRows(), MMAP.get(2 + c).getNumCols(),
                            MMAP.get(2 + c).getNumRows() * MMAP.get(2 + c).getNumCols()));
                }
            }

            // MATLAB mmap_scale returns right after this normalization in the
            // per-class (vector) case: the iterative refinement below that
            // point is unreachable there, so no refinement is applied here
            // either. A previous reachable refinement produced invalid MMAPs
            // that failed the QBD (sub)stochasticity check downstream.
            SCALED = Mmap_normalize.mmap_normalize(SCALED);
        }
        return SCALED;
    }

    /**
     * Overloaded function for backward compatibility.
     */
    public static MatrixCell mmap_scale(MatrixCell MMAP, Matrix M) {
        return mmap_scale(MMAP, M, 30);
    }
}
