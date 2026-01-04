package jline.api.sn

import jline.lang.NetworkStruct

/**
 * Checks if the network has one or more open classes
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return boolean
 */

fun snHasOpenClasses(sn: NetworkStruct): Boolean {
    return sn.njobs.hasInfinite()
}
/**
 * Stochastic network HasOpenClasses algorithms
 */
@Suppress("unused")
class SnhasopenclassesAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}