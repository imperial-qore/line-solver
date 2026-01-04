/**
 * @file Markovian Arrival MAP with Marked arrivals gamma forward-backward MMAP fitting
 * 
 * Fits MAMAP with autocorrelation control using forward-backward moments from MMAP input.
 * Combines correlation structure preservation with statistical moment matching.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Fits a second-order acyclic MMAP[m] to match the characteristics of the input MMAP.
 * This function extracts moments, class probabilities, and forward/backward moments
 * from the input MMAP and returns a fitted approximation.
 *
 * @param mmap Input MMAP (Markovian Arrival Process with marked arrivals) stored as MatrixCell
 * @return Fitted second-order MAMAP that matches the input characteristics
 */
fun mamap2m_fit_gamma_fb_mmap(mmap: MatrixCell): MatrixCell {
    // Extract moments and characteristics from input MMAP
    val M1 = map_moment(mmap[0], mmap[1], 1)
    val M2 = map_moment(mmap[0], mmap[1], 2)
    val M3 = map_moment(mmap[0], mmap[1], 3)
    val GAMMA = map_gamma(mmap)

    val P = mmap_pc(mmap)
    val moments = Matrix(1, 1)
    moments[0, 0] = 1.0
    val F = mmap_forward_moment(mmap, moments)
    val B = mmap_backward_moment(mmap, moments)
    
    return mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P.toArray1D(), F.toArray1D(), B.toArray1D())
}

/**
 * Computes the second-order MAMAP[m] fitting the given ordinary moments,
 * autocorrelation decay rate, class probabilities, forward moments, and backward moments.
 *
 * @param M1 First moment of inter-arrival times
 * @param M2 Second moment of inter-arrival times
 * @param M3 Third moment of inter-arrival times
 * @param GAMMA Auto-correlation decay rate
 * @param P Class probabilities
 * @param F First-order forward moments
 * @param B First-order backward moments
 * @return Fitted MAMAP[m]
 */
fun mamap2m_fit_gamma_fb(M1: Double, M2: Double, M3: Double, GAMMA: Double, 
                         P: DoubleArray, F: DoubleArray, B: DoubleArray): MatrixCell {
    
    // Ensure F and B are column vectors
    val Fcol = F.copyOf()
    val Bcol = B.copyOf()
    
    // Fit underlying AMAP(2) - up to four equivalent representations may be found
    val (_, MAPS) = amap2_fit_gamma(M1, M2, M3, GAMMA)
    
    // Handle marked Poisson process case
    if (MAPS.size == 1 && MAPS[0][0].numRows == 1) {
        val MAP = MAPS[0]
        val m = P.size
        
        // Fit class probabilities
        val result = MatrixCell(2 + m)
        result[0] = MAP[0].copy()
        result[1] = MAP[1].copy()
        
        for (c in 0 until m) {
            result[2 + c] = result[1].scale(P[c])
        }
        
        return result
    }
    
    // If no valid AMAPs found, return a simple approximation
    if (MAPS.isEmpty()) {
        return createFallbackMMAP(M1, P)
    }
    
    // Use the underlying AMAP(2) form which produces the least error
    val MMAPS = mutableListOf<MatrixCell>()
    val ERRORS = mutableListOf<Double>()
    
    for (j in MAPS.indices) {
        try {
            val (fittedMMAP, fF, fB) = mamap2m_fit_fb_multiclass(MAPS[j], P, Fcol, Bcol)
            MMAPS.add(fittedMMAP)
            
            // Compute fitting error
            val forwardError = computeRelativeError(fF, Fcol)
            val backwardError = computeRelativeError(fB, Bcol)
            ERRORS.add(forwardError + backwardError)
            
        } catch (e: Exception) {
            // If fitting fails for this AMAP, skip it
            MMAPS.add(createFallbackMMAP(M1, P))
            ERRORS.add(Double.MAX_VALUE)
        }
    }
    
    // Return the MMAP with minimum error
    val bestIndex = ERRORS.indices.minByOrNull { ERRORS[it] } ?: 0
    return MMAPS[bestIndex]
}

/**
 * Creates a fallback MMAP when fitting fails
 */
private fun createFallbackMMAP(M1: Double, P: DoubleArray): MatrixCell {
    val m = P.size
    val result = MatrixCell(2 + m)
    
    // Create simple Poisson process
    result[0] = Matrix(1, 1)
    result[0][0, 0] = -1.0 / M1
    result[1] = Matrix(1, 1)
    result[1][0, 0] = 1.0 / M1
    
    // Add class probabilities
    for (c in 0 until m) {
        result[2 + c] = result[1].scale(P[c])
    }
    
    return result
}

/**
 * Computes relative error between fitted and target moments
 */
private fun computeRelativeError(fitted: DoubleArray, target: DoubleArray): Double {
    var error = 0.0
    for (i in fitted.indices) {
        if (target[i] != 0.0) {
            val relError = (fitted[i] / target[i] - 1.0)
            error += relError * relError
        }
    }
    return error
}
/**
 * MAMAP 2m fit gamma fb mmap algorithms
 */
@Suppress("unused")
class Mamap2mFitGammaFbMmapAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}