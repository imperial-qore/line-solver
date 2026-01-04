package jline.api.mam

import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import java.util.*

/**
 * APH random generation algorithms.
 * 
 * Provides methods for generating random Acyclic Phase-type (APH) distributions.
 * APH distributions are special cases of phase-type distributions with acyclic structure,
 * commonly used to model service times and inter-arrival times in queueing systems.
 *
 * @since LINE 3.0
 */

/**
 * Generates a random Acyclic Phase-type (APH) distribution with K phases.
 *
 * An APH distribution is a special case of a phase-type distribution where the
 * underlying Markov chain has an acyclic structure (no cycles). This function
 * generates a random APH by creating matrices with the appropriate structure.
 *
 * @param K the number of phases
 * @return a MatrixCell representing the random APH distribution
 */
fun aph_rand(K: Int): MatrixCell {
    val random = Random()
    
    // Generate random D1 matrix (full matrix)
    val D1 = Matrix(K, K)
    for (i in 0 until K) {
        for (j in 0 until K) {
            D1.set(i, j, random.nextDouble())
        }
    }
    
    // Generate random D0 matrix with acyclic structure
    val D0 = Matrix(K, K)
    for (i in 0 until K) {
        for (j in 0 until K) {
            if (j < i) {
                // Set elements below diagonal to 0 for acyclic structure
                D0.set(i, j, 0.0)
            } else {
                D0.set(i, j, random.nextDouble())
            }
        }
    }
    
    // Create temporary MAP
    val tempMAP = MatrixCell()
    tempMAP[0] = D0
    tempMAP[1] = D1
    
    // Apply renewal process and normalize
    val renewalMAP = map_renewal(tempMAP)
    return map_normalize(renewalMAP)
}

/**
 * Generates a random Acyclic Phase-type (APH) distribution with 2 phases.
 *
 * @return a MatrixCell representing the random APH distribution
 */
fun aph_rand(): MatrixCell {
    return aph_rand(2)
}
/**
 * APH  rand algorithms
 */
@Suppress("unused")
class AphRandAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}