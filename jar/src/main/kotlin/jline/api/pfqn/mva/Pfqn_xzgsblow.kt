/**
 * @file Lower Geometric Square-Root Bound (GSB) for throughput in single-class networks
 * 
 * Computes the lower GSB for throughput in closed single-class queueing networks using
 * geometric mean analysis and square-root approximations. Provides tighter bounds than
 * basic asymptotic approximations through advanced geometric techniques.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Computes the lower Geometric Square-Root Bound (GSB) for the throughput of the given closed single-class
 * queueing networks
 *
 * @param L - service demand matrix
 * @param N - population
 * @param Z - think time
 * @return - the lower GSB for the throughput
 */

fun pfqn_xzgsblow(L: Matrix, N: Double, Z: Double): Double {
    val M = L.length()
    val maxL = L.elementMax()
    var R = Z + L.elementSum() + maxL * (N - 1)
    for (i in 0..<M) {
        if (L[i] < maxL) {
            R += (L[i] - maxL) * pfqn_qzgblow(L, N - 1, Z, i)
        }
    }
    return 2 * N / (R + FastMath.sqrt(R * R - 4 * Z * maxL * (N - 1)))
}
/**
 * PFQN xzgsblow algorithms
 */
@Suppress("unused")
class PfqnXzgsblowAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}