/**
 * @file CTMC random generator matrix construction
 * 
 * Generates random infinitesimal generator matrices for continuous-time Markov chains 
 * with specified properties. Used for testing CTMC algorithms and creating benchmark 
 * problems with controllable statistical characteristics.
 * 
 * @since LINE 3.0
 */
package jline.api.mc

import jline.util.matrix.Matrix

/**
 * Form a random infinitesimal generator of a CTMC
 *
 * @param length size of random matrix
 * @return Infinitesimal generator of CTMC
 */

fun ctmc_rand(length: Int): Matrix {
    val rand_matrix = Matrix(length, length)
    rand_matrix.randMatrix(length)
    return ctmc_makeinfgen(rand_matrix)
}
/**
 * CTMC rand algorithms
 */
@Suppress("unused")
class CtmcRandAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}