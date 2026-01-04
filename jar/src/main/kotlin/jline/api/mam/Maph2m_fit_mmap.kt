/**
 * @file Multi-class Absorbing Phase-type distribution MMAP-based fitting
 * 
 * Fits MAPH(2,m) by approximating characteristics of input MMAP processes.
 * Used for converting multiclass arrival processes to phase-type service representations.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell

/**
 * Fits an MMAP[m] with a second-order MAPH[m] that matches the class
 * probabilities (always fitted exactly) and the backward moments.
 *
 * @param mmap MMAP to fit (of arbitrary order)
 * @return Fitted second-order MAPH[m]
 */
fun maph2m_fit_mmap(mmap: MatrixCell): MatrixCell {
    val M1 = map_moment(mmap, 1)
    val M2 = map_moment(mmap, 2)
    val M3 = map_moment(mmap, 3)
    
    val P = mmap_pc(mmap) // Class probabilities
    val moments = jline.util.matrix.Matrix(1, 1)
    moments[0, 0] = 1.0
    val B = mmap_backward_moment(mmap, moments) // Backward moments
    
    return maph2m_fit(M1, M2, M3, P, B)
}
/**
 * MAPH 2m fit mmap algorithms
 */
@Suppress("unused")
class Maph2mFitMmapAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}