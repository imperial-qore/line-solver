/**
 * @file Multi-class Absorbing Phase-type distribution MMAP-based fitting
 *
 * Fits MAPH(2,m) by approximating characteristics of input MMAP processes.
 * Used for converting multiclass arrival processes to phase-type service representations.
 *
 * @since LINE 3.0
 */
package jline.api.mam;

import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Maph2m_fit_mmap {
    private Maph2m_fit_mmap() {}

    /**
     * Fits an MMAP[m] with a second-order MAPH[m] that matches the class
     * probabilities (always fitted exactly) and the backward moments.
     *
     * @param mmap MMAP to fit (of arbitrary order)
     * @return Fitted second-order MAPH[m]
     */
    public static MatrixCell maph2m_fit_mmap(MatrixCell mmap) {
        double M1 = Map_moment.map_moment(mmap, 1);
        double M2 = Map_moment.map_moment(mmap, 2);
        double M3 = Map_moment.map_moment(mmap, 3);

        Matrix P = Mmap_pc.mmap_pc(mmap);
        Matrix moments = new Matrix(1, 1);
        moments.set(0, 0, 1.0);
        Matrix B = Mmap_backward_moment.mmap_backward_moment(mmap, moments);

        return Maph2m_fit.maph2m_fit(M1, M2, M3, P, B);
    }
}
