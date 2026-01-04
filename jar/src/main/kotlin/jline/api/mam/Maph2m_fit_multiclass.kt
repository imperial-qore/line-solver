/**
 * @file Multi-class Absorbing Phase-type distribution multiclass fitting
 * 
 * Fits MAPH(2,m) models to multiclass characteristics with class-specific parameters.
 * Specialized fitting for complex multiclass phase-type service time distributions.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell

/**
 * Fits a multi-class MAPH(2,m) model to given multi-class characteristics.
 *
 * @param M1 First moment of inter-arrival times
 * @param M2 Second moment of inter-arrival times
 * @param M3 Third moment of inter-arrival times
 * @param classProbs Array of class probabilities
 * @param classRates Array of class-specific rates (optional)
 * @param backwardMoments Array of backward moments for each class (optional)
 * @return Fitted MAPH(2,m) model
 */
fun maph2m_fit_multiclass(
    M1: Double,
    M2: Double, 
    M3: Double,
    classProbs: DoubleArray,
    classRates: DoubleArray? = null,
    backwardMoments: DoubleArray? = null
): MatrixCell {
    
    val m = classProbs.size
    
    // Ensure class probabilities sum to 1
    val normalizedProbs = classProbs.clone()
    val probSum = normalizedProbs.sum()
    if (Math.abs(probSum - 1.0) > 1e-6) {
        for (i in normalizedProbs.indices) {
            normalizedProbs[i] /= probSum
        }
    }
    
    // If class rates not provided, distribute evenly based on probabilities
    val rates = classRates ?: normalizedProbs.map { it / M1 }.toDoubleArray()
    
    // If backward moments not provided, use uniform distribution
    val backMoments = backwardMoments ?: DoubleArray(m) { M1 }
    
    // Convert arrays to matrices
    val probMatrix = jline.util.matrix.Matrix(1, normalizedProbs.size)
    for (i in normalizedProbs.indices) {
        probMatrix[0, i] = normalizedProbs[i]
    }
    
    val backMatrix = jline.util.matrix.Matrix(backMoments.size, 1)
    for (i in backMoments.indices) {
        backMatrix[i, 0] = backMoments[i]
    }
    
    return maph2m_fit(M1, M2, M3, probMatrix, backMatrix)
}
/**
 * MAPH 2m fit multiclass algorithms
 */
@Suppress("unused")
class Maph2mFitMulticlass {
    companion object {
        // Class documentation marker for Dokka
    }
}