/**
 * @file Lower Geometric Bound computation for queue lengths in single-class networks
 * 
 * Computes the lower Geometric Bound (GB) for queue lengths in closed single-class queueing
 * networks. Provides tight lower bounds for performance analysis and serves as input for
 * other approximation algorithms requiring queue length bounds.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.util.matrix.Matrix
import kotlin.math.pow

/**
 * Computes the lower Geometric Bound (GB) for the queue length of the given closed single-class
 * queueing networks
 *
 * @param L - service demand matrix
 * @param N - population
 * @param Z - think time
 * @param i - station index
 * @return - the lower GB for the queue length
 */

fun pfqn_qzgblow(L: Matrix, N: Double, Z: Double, i: Int): Double {
    val Yi = N * L[i] / (Z + L.elementSum() + L.elementMax() * N)
    val Qgb = Yi / (1 - Yi) - (Yi.pow(N + 1) / (1 - Yi))
    return Qgb
}
/**
 * PFQN qzgblow algorithms
 */
@Suppress("unused")
class PfqnQzgblowAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}