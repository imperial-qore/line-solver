/**
 * @file DTMC random transition matrix generation
 * 
 * Generates random stochastic transition matrices for discrete-time Markov chains 
 * with controllable properties. Essential for algorithm testing, benchmarking, 
 * and creating synthetic DTMC models with specified characteristics.
 * 
 * @since LINE 3.0
 */
package jline.api.mc

import jline.util.matrix.Matrix

/**
 * Form a random infinitesimal generator of a DTMC
 *
 * @param length size of random matrix
 * @return Infinitesimal generator of CTMC
 */

fun dtmc_rand(length: Int): Matrix {
    val rand_matrix = Matrix(length, length)
    rand_matrix.randMatrix(length)
    return dtmc_makestochastic(rand_matrix)
}
/**
 * DTMC rand algorithms
 */
@Suppress("unused")
class DtmcRandAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}