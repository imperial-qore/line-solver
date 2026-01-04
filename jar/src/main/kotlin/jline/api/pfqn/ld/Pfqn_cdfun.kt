/**
 * @file Class-dependent scaling function evaluator for load-dependent queueing systems
 * 
 * Provides functionality to evaluate class-dependent scaling functions in load-dependent queueing networks.
 * Calculates scaling factors based on queue-length dependent service rates, supporting state-dependent
 * service mechanisms where service capacity varies with the number of jobs present.
 * 
 * @since LINE 3.0
 */
package jline.api.pfqn.ld

import jline.util.SerializableFunction
import jline.util.matrix.Matrix

/**
 * Evaluate class-dependent (CD) scaling function
 *
 * @param nvec      Per-class queue-length values. The values can be continuous.
 * @param cdscaling CD functions indexed by station index, or null if none
 * @param M         Number of stations
 * @return matrix of scaling factors
 */

fun pfqn_cdfun(nvec: Matrix,
               cdscaling: List<SerializableFunction<Matrix?, Double>?>?,
               M: Int): Matrix {
    val r = Matrix(M, 1)
    r.fill(1.0)
    if (!(cdscaling == null || cdscaling.isEmpty())) {
        for (i in 0..<M) {
            val func = cdscaling.getOrNull(i)
            if (func != null) {
                r[i, 0] = 1.0 / func.apply(Matrix.extractRows(nvec, i, i + 1, null))
            }
        }
    }
    return r
}
/**
 * PFQN cdfun algorithms
 */
@Suppress("unused")
class PfqnCdfunAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}