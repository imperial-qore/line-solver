package jline.api.sn

import jline.lang.NetworkStruct

/**
 * Checks if the network has closed     classes with non-integer populations
 *
 * @param sn       - NetworkStruct object for the queueing network model
 * @return boolean
 */

fun snHasFractionalPopulations(sn: NetworkStruct): Boolean {
    return !sn.njobs.isInteger
}
/**
 * Stochastic network HasFractionalPopulations algorithms
 */
@Suppress("unused")
class SnhasfractionalpopulationsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}