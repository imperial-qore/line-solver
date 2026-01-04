/**
 * @file Upper Geometric Bound computation for queue lengths in single-class networks
 * 
 * Computes the upper Geometric Bound (GB) for queue lengths in closed single-class queueing
 * networks. Provides tight upper bounds for performance analysis and bounding techniques
 * in approximate MVA algorithms.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.util.Maths
import jline.util.matrix.Matrix
import kotlin.math.pow

/**
 * Computes the upper Geometric Bound (GB) for the queue length of the given closed single-class
 * queueing networks
 *
 * @param L - service demand matrix
 * @param N - population
 * @param Z - think time
 * @param i - station index
 * @return - the upper GB for the queue length
 */

fun pfqn_qzgbup(L: Matrix, N: Double, Z: Double, i: Int): Double {
    val sumL = L.elementSum()
    val sumLSq = L.elementMult(L, null).elementSum()
    val sigma = sumLSq / sumL
    val Yi = L[i] * Maths.min(1.0 / L.elementMax(), N / (Z + sumL + sigma * (N - 1 - Z * pfqn_xzabaup(L, N - 1, Z))))
    val Qgb = if (Yi < 1) {
        Yi / (1 - Yi) - (Yi.pow((N + 1)) / (1 - Yi))
    } else {
        N
    }
    return Qgb
}
/**
 * PFQN qzgbup algorithms
 */
@Suppress("unused")
class PfqnQzgbupAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}