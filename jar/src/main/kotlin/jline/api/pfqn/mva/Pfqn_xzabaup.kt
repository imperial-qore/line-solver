/**
 * @file Upper Asymptotic Bound Approximation (ABA) for throughput in single-class networks
 * 
 * Computes the upper ABA bound for throughput in closed single-class queueing networks.
 * Provides fundamental upper bounds based on bottleneck analysis and system capacity
 * constraints for performance evaluation.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.mva

import jline.util.Maths
import jline.util.matrix.Matrix

/**
 * Computes the upper ABA for the throughput of the given closed single-class queueing networks
 *
 * @param L - service demand matrix
 * @param N - population
 * @param Z - think time
 * @return - the upper ABA for the throughput
 */

fun pfqn_xzabaup(L: Matrix, N: Double, Z: Double): Double {
    val e1 = 1 / L.elementMax()
    val e2 = N / (L.elementSum() + Z)
    return Maths.min(e1, e2)
}
/**
 * PFQN xzabaup algorithms
 */
@Suppress("unused")
class PfqnXzabaupAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}