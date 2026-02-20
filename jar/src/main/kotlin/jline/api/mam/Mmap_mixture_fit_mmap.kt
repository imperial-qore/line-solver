package jline.api.mam

import jline.io.Ret
import jline.util.matrix.MatrixCell

/**
 * Fits a mixture of Markovian Arrival Processes (MMAPs) to match the given moments.
 *
 *
 * This method first computes the second-order cross-moments (`P2`) and the first three cross-moments (`M1`, `M2`, `M3`)
 * of the given MMAP. It then fits a mixture of MMAPs to match these moments, returning the parameters of the fitted mixture.
 *
 * @param mmap the original MMAP to be fitted
 * @return a `mmap_mixture_fit_return_type` containing the fitted mixture parameters
 */
fun mmap_mixture_fit_mmap(mmap: MatrixCell): Ret.mamMMAPMixtureFit {
    val P2 = mmap_sigma2_cell(mmap)
    val M1 = mmap_cross_moment(mmap, 1)
    val M2 = mmap_cross_moment(mmap, 2)
    val M3 = mmap_cross_moment(mmap, 3)

    return mmap_mixture_fit(P2, M1, M2, M3)
}
/**
 * MMAP mixture fit mmap algorithms
 */
@Suppress("unused")
class MmapMixtureFitMmapAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}