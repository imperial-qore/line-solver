package jline.api.sn

import jline.lang.NetworkStruct

/**
 * Checks if the network uses fork and/or join nodes
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return boolean
 */
fun snHasForkJoin(sn: NetworkStruct): Boolean {
    for (i in 0..<sn.fj.numRows) {
        for (j in 0..<sn.fj.numCols) {
            if (sn.fj[i, j] > 0) return true
        }
    }
    return false
}
/**
 * Stochastic network HasForkJoin algorithms
 */
@Suppress("unused")
class SnhasforkjoinAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}