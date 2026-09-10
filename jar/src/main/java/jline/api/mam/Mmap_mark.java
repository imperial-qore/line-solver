package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_mark {
    private Mmap_mark() {}

    /**
     * Converts a Markovian Arrival Process with marked arrivals (MMAP) into a new MMAP
     * with redefined classes based on a given probability matrix.
     *
     * @param MMAP the original MMAP with K types
     * @param prob a KxR matrix describing the probability of a type-k arrival in the
     *             original MMAP being marked as a type-r arrival in the new MMAP
     * @return a new MMAP with R types (classes)
     */
    public static MatrixCell mmap_mark(MatrixCell MMAP, Matrix prob) {
        int K = prob.getNumRows();
        int R = prob.getNumCols();
        MatrixCell mmap = new MatrixCell(2 + R);
        mmap.set(0, MMAP.get(0).copy());
        mmap.set(1, MMAP.get(1).copy());
        for (int r = 0; r < R; r++) {
            mmap.set(2 + r, new Matrix(MMAP.get(0).length(), MMAP.get(0).length(),
                    MMAP.get(0).length() * MMAP.get(0).length()));
            for (int k = 0; k < K; k++) {
                Matrix a = MMAP.get(2 + k).copy();
                a.scaleEq(prob.get(k, r));
                mmap.set(2 + r, mmap.get(2 + r).add(1.0, a));
            }
        }

        return mmap;
    }
}
