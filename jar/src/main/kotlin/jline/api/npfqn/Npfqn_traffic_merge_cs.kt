/**
 * @file NPFQN Traffic Merging with Class Switching
 * 
 * Implements traffic merging algorithms for non-product-form queueing networks
 * with class switching capabilities. Handles the aggregation of multiple traffic
 * streams while accounting for class transitions in NPFQN analysis.
 * 
 * @since LINE 3.0
 */
package jline.api.npfqn

import jline.api.mam.mmap_mark
import jline.api.mam.mmap_super
import jline.util.matrix.Matrix
import jline.util.matrix.MatrixCell

/**
 * Merges MMAP traffic flows with class switching using default parameters.
 *
 * @param MMAPs Map of MMAP traffic flows
 * @param prob  Probability matrix for class switching
 * @return Merged and normalized MMAP traffic flow
 */
fun npfqn_traffic_merge_cs(MMAPs: Map<Int, MatrixCell>, prob: Matrix, config: String = "default"): MatrixCell? {
    val n = MMAPs.size
    val R = prob.numCols
    val MMAP_copy: MutableMap<Int, MatrixCell> = HashMap()
    for ((key, innerMap) in MMAPs) {
        MMAP_copy[key] = MatrixCell(innerMap)
    }
    for (i in 0..<n) {
        val P = Matrix(R, R, R * R)
        for (r in 0..<R) {
            for (s in 0..<R) {
                P[r, s] = prob[(i - 1) * R + r, s]
            }
        }
        MMAP_copy[i] = mmap_mark(MMAP_copy[i]!!, P)
    }
    var SMMAP: MatrixCell? = MatrixCell()
    if (n == 1) {
        SMMAP = MMAPs[0]
    } else {
        if (config == "default" || config == "super") {
            SMMAP = MMAPs[0]
            for (j in 1..<n) {
                SMMAP = mmap_super(SMMAP!!, MMAP_copy[j]!!, "match")
            }
        }
    }
    return SMMAP
}