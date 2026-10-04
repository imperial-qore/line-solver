/**
 * @file Markovian Arrival Process marking for multiclass processes
 *
 * Creates Marked MAP (MMAP) representations by adding class labels to MAP arrivals.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

import static jline.io.InputOutput.line_warning;
import static jline.io.InputOutput.mfilename;

public final class Map_mark {
    private Map_mark() {}

    /**
     * Creates a Marked Markovian Arrival Process (MMAP) by marking a given MAP with additional
     * phases based on specified marking probabilities.
     *
     * @param MAP  The original Markovian Arrival Process stored in a MatrixCell.
     * @param prob A matrix containing the marking probabilities.
     * @return A MatrixCell representing the Marked Markovian Arrival Process (MMAP).
     */
    public static MatrixCell map_mark(MatrixCell MAP, Matrix prob) {
        Matrix prob_local = prob.copy();
        if ((prob_local.elementSum() < 1 - 1e-6) || (prob_local.elementSum() > 1 + 1e-6)) {
            line_warning(mfilename(new Object() {}),
                    "Input marking probabilities do not sum to 1. Normalizing.");
        }
        prob_local.scaleEq(prob_local.elementSum());
        int R = prob_local.getNumCols();
        MatrixCell mmap = new MatrixCell(2 + R);
        mmap.set(0, MAP.get(0).copy());
        mmap.set(1, MAP.get(1).copy());
        for (int r = 0; r < R; r++) {
            mmap.set(2 + r, Matrix.createLike(MAP.get(1)));
            Matrix a = MAP.get(1).copy();
            a.scaleEq(prob_local.get(r));
            mmap.set(2 + r, a);
        }
        return mmap;
    }
}
