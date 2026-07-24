/**
 * @file Marked Markovian Arrival Process class hiding operations
 *
 * Hides specified arrival classes in MMAP processes by removing observable events.
 * Used for model reduction and analyzing subsystems in multiclass arrival processes.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_hide {
    private Mmap_hide() {}

    /**
     * Hides specified types of arrivals in a Markovian Arrival Process with marked arrivals (MMAP).
     *
     * @param MMAP  the original MMAP
     * @param types a matrix containing the indices of the types to be hidden
     * @return a new MMAP with the specified types hidden
     */
    public static MatrixCell mmap_hide(MatrixCell MMAP, Matrix types) {
        // Copy ALL cells first (kept classes retain their matrices), then zero
        // the hidden types; copying only {D0,D1} would leave the kept classes'
        // cells null and crash mmap_normalize (MATLAB: MMAP{2+k} = 0*MMAP{2+k}).
        MatrixCell mmap = new MatrixCell();
        for (int i = 0; i < MMAP.size(); i++) {
            mmap.set(i, MMAP.get(i).copy());
        }
        // The number of types to hide is the element count of the row/col
        // vector, NOT Matrix.length() = max(rows,cols): a degenerate 1x0 vector
        // (K == 1, so setdiff(1:K,r) is empty) has length() == 1 and would
        // dereference types.get(0) out of bounds. MATLAB's `for k=types(:)'`
        // simply does not iterate over an empty set.
        int nTypes = types.getNumRows() * types.getNumCols();
        for (int i = 0; i < nTypes; i++) {
            mmap.set((int) (2 + types.get(i)), new Matrix(MMAP.get(0).getNumRows(), MMAP.get(0).getNumRows()));
        }

        return Mmap_normalize.mmap_normalize(mmap);
    }
}
