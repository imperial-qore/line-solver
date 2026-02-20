/**
 * @file Upper Geometric Square-Root Bound (GSB) for throughput in single-class networks
 * 
 * Computes the upper GSB for throughput in closed single-class queueing networks using
 * geometric mean analysis and square-root approximations. Provides enhanced upper bounds
 * through advanced geometric bounding techniques.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.util.matrix.Matrix
import org.apache.commons.math3.util.FastMath

/**
 * Computes the upper Geometric Square-Root Bound (GSB) for the throughput of the given closed single-class
 * queueing networks
 *
 * @param L - service demand matrix
 * @param N - population
 * @param Z - think time
 * @return - the upper GSB for the throughput
 */

fun pfqn_xzgsbup(L: Matrix, N: Double, Z: Double): Double {
    val M = L.length()
    val maxL = L.elementMax()
    var R = Z + L.elementSum() + maxL * (N - 1)
    for (i in 0..<M) {
        if (L[i] < maxL) {
            R += (L[i] - maxL) * pfqn_qzgbup(L, N - 1, Z, i)
        }
    }
    return 2 * N / (R + FastMath.sqrt(R * R - 4 * Z * maxL * N))
}
/**
 * PFQN xzgsbup algorithms
 */
@Suppress("unused")
class PfqnXzgsbupAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}