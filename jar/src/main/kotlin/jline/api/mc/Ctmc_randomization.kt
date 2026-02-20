/**
 * @file CTMC randomization (uniformization) method
 * 
 * Converts a continuous-time Markov chain into an equivalent discrete-time chain
 * using the randomization technique. This transformation enables the use of DTMC
 * algorithms for CTMC analysis and is fundamental to numerical CTMC solution methods.
 *
 * @since LINE 3.0
 */
package jline.api.mc

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath
import kotlin.random.Random

/**
 * Convert a CTMC to a DTMC using randomization technique
 *
 * @param Q Infinitesimal generator matrix of the CTMC
 * @param q Optional uniformization rate. If not provided, it will be computed automatically
 * @return Pair containing the transition matrix P and the uniformization rate q
 */
fun ctmc_randomization(Q: Matrix, q: Double? = null): Pair<Matrix, Double> {
    val uniformizationRate = q ?: (Q.elementMaxAbs() + Random.nextDouble())
    
    val n = Q.numRows
    val I = Matrix.eye(n)
    val P = Q.scale(1.0 / uniformizationRate).add(I)
    
    val stochasticP = dtmc_makestochastic(P)
    
    return Pair(stochasticP, uniformizationRate)
}
/**
 * CTMC randomization algorithms
 */
@Suppress("unused")
class CtmcRandomizationAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}