/**
 * @file NPFQN Traffic Splitting with Class Switching
 * 
 * Implements traffic splitting algorithms for non-product-form queueing networks
 * with class switching capabilities. Handles the decomposition of traffic streams
 * while accounting for class transitions in NPFQN analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.npfqn

import jline.api.mam.mmap_normalize
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell
import org.apache.commons.math3.util.FastMath

/**
 * NPFQN traffic splitting with class switching algorithms.
 * 
 * Provides methods for splitting MMAP traffic flows in non-product-form queueing networks
 * with class switching capabilities. This enables modeling of complex routing behaviors
 * where customers can change service classes as they move through the network.
 *
 * @since LINE 3.0
 */

/**
 * Splits MMAP traffic flows with class switching.
 *
 * @param MMAP MMAP to be split
 * @param P    Class switching matrix
 * @return Map of split MMAPs
 */
fun npfqn_traffic_split_cs(MMAP: MatrixCell, P: Matrix): Map<Int, MatrixCell?> {
    MMAP.size()
    val R = P.numRows
    val J = P.numCols
    val M = FastMath.round(J / R.toFloat())
    val SMMAP: MutableMap<Int, MatrixCell?> = HashMap()
    for (jst in 0..<M) {
        SMMAP[jst] = MatrixCell()
        SMMAP[jst]!![0] = MMAP[0].add(1.0, MMAP[1])
        SMMAP[jst]!![1] = Matrix(MMAP[1].numRows, MMAP[1].numCols)
        for (s in 0..<R) {
            SMMAP[jst]!![2 + s] = Matrix(SMMAP[jst]!![0].numRows, SMMAP[jst]!![0].numCols)
            for (r in 0..<R) {
                val a = MMAP[2 + r].copy()
                a.scaleEq(P[r, (jst - 1) * R + s])
                SMMAP[jst]!![2 + s] = SMMAP[jst]!![2 + s].add(1.0, a)
                SMMAP[jst]!![1] = SMMAP[jst]!![1].add(1.0, a)
                SMMAP[jst]!![0] = SMMAP[jst]!![0].add(1.0, a)
            }
        }
        SMMAP[jst] = mmap_normalize(SMMAP[jst]!!)
    }
    return SMMAP
}