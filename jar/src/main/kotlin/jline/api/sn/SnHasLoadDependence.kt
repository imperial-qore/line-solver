package jline.api.sn

import jline.lang.NetworkStruct

/**
 * Checks if the network has a station with load-dependent service process
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return boolean
 */

fun snHasLoadDependence(sn: NetworkStruct): Boolean {
    return sn.lldscaling.numCols > 0
}
/**
 * Stochastic network HasLoadDependence algorithms
 */
@Suppress("unused")
class SnhasloaddependenceAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}