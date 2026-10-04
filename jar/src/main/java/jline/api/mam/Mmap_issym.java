package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_issym {
    private Mmap_issym() {}

    /**
     * Checks if an MMAP is symmetric.
     * An MMAP is symmetric if all its matrices are symmetric.
     *
     * @param mmap The MMAP to check
     * @param tolerance Tolerance for numerical comparison
     * @return True if the MMAP is symmetric, false otherwise
     */
    public static boolean mmap_issym(MatrixCell mmap, double tolerance) {
        // Check each matrix in the MMAP
        for (int i = 0; i < mmap.size(); i++) {
            Matrix matrix = mmap.get(i);

            // Check if matrix is square
            if (matrix.getNumRows() != matrix.getNumCols()) {
                return false;
            }

            // Check if matrix is symmetric
            for (int row = 0; row < matrix.getNumRows(); row++) {
                for (int col = 0; col < matrix.getNumCols(); col++) {
                    double diff = Math.abs(matrix.get(row, col) - matrix.get(col, row));
                    if (diff > tolerance) {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    public static boolean mmap_issym(MatrixCell mmap) {
        return mmap_issym(mmap, 1e-12);
    }
}
