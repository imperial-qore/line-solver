package jline.api.sn

import jline.lang.NetworkStruct

/**
 * Checks if the network uses class priorities
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return boolean
 */
fun snHasPriorities(sn: NetworkStruct): Boolean {
    for (i in 0..<sn.classprio.numRows) {
        for (j in 0..<sn.classprio.numCols) {
            if (sn.classprio[i, j] > 0) return true
        }
    }
    return false
}
/**
 * Stochastic network HasPriorities algorithms
 */
@Suppress("unused")
class SnhasprioritiesAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}