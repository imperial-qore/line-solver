/**
 * @file Markovian Arrival Process kurtosis computation
 * 
 * Computes kurtosis of MAP inter-arrival times measuring tail heaviness and distribution shape.
 * Important for characterizing extreme behavior and heavy-tailed properties in arrival processes.
 * 
 * @since LINE 3.0
 */
package jline.api.mam

import jline.util.matrix.MatrixCell

/**
 * Computes the kurtosis of the inter-arrival times in a Markovian Arrival Process (MAP).
 *
 * The kurtosis is computed using the formula:
 * KURT = (m4 - 4*m3*m1 + 6*m2*m1^2 - 3*m1^4) / Var(X)^2
 *
 * where mi is the i-th moment and Var(X) is the variance of the inter-arrival times.
 *
 * @param MAP The Markovian Arrival Process stored in a MatrixCell, containing the (D0, D1) matrices
 * @return The kurtosis of the inter-arrival times
 */
fun map_kurt(MAP: MatrixCell): Double {
    // Compute the first four moments
    val m1 = map_moment(MAP, 1)
    val m2 = map_moment(MAP, 2)
    val m3 = map_moment(MAP, 3)
    val m4 = map_moment(MAP, 4)
    
    // Compute variance
    val variance = map_var(MAP)
    
    // Compute kurtosis using the formula
    val numerator = m4 - 4.0 * m3 * m1 + 6.0 * m2 * m1 * m1 - 3.0 * m1 * m1 * m1 * m1
    val denominator = variance * variance
    
    return numerator / denominator
}
/**
 * MAP kurt algorithms
 */
@Suppress("unused")
class MapKurtAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}