package jline.api.mam;

import jline.io.Ret;
import jline.util.matrix.Matrix;
import jline.util.matrix.MatrixCell;

public final class Mmap_mixture_fit_mmap {
    private Mmap_mixture_fit_mmap() {}

    /**
     * Fits a mixture of Markovian Arrival Processes (MMAPs) to match the given moments.
     *
     * This method first computes the second-order cross-moments (`P2`) and the first three cross-moments (`M1`, `M2`, `M3`)
     * of the given MMAP. It then fits a mixture of MMAPs to match these moments, returning the parameters of the fitted mixture.
     *
     * @param mmap the original MMAP to be fitted
     * @return a `mamMMAPMixtureFit` containing the fitted mixture parameters
     */
    public static Ret.mamMMAPMixtureFit mmap_mixture_fit_mmap(MatrixCell mmap) {
        MatrixCell P2 = Mmap_sigma2.mmap_sigma2_cell(mmap);
        Matrix M1 = Mmap_cross_moment.mmap_cross_moment(mmap, 1);
        Matrix M2 = Mmap_cross_moment.mmap_cross_moment(mmap, 2);
        Matrix M3 = Mmap_cross_moment.mmap_cross_moment(mmap, 3);

        return Mmap_mixture_fit.mmap_mixture_fit(P2, M1, M2, M3);
    }
}
