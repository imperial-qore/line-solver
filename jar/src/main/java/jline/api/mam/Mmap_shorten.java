package jline.api.mam;

import jline.util.matrix.MatrixCell;

public final class Mmap_shorten {
    private Mmap_shorten() {}

    /**
     * Converts an MMAP representation from M3A format to BUTools format.
     *
     * <p>In the M3A format, an MMAP is represented as a MAP followed by multiple D1 matrices for different markings (D0, D1, D1a, D1b, ...).
     * In the BUTools format, the MMAP representation skips the initial D1 matrix, directly listing the marking matrices (D0, D1a, D1b, ...).
     * This method reorders the matrices accordingly.
     *
     * @param mmap the MatrixCell containing the MMAP representation in M3A format
     * @return a MatrixCell representing the MMAP in BUTools format
     */
    public static MatrixCell mmap_shorten(MatrixCell mmap) {
        MatrixCell result = new MatrixCell();
        for (int i = 0; i < mmap.size(); i++) {
            if (i == 0) {
                result.set(0, mmap.get(0));
            } else if (i != 1) {
                result.set(i - 1, mmap.get(i));
            }
        }
        return result;
    }
}
