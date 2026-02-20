package jline.api.sn

import jline.lang.NetworkStruct

/**
 * Checks if the network has one or more closed classes
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return boolean
 */
fun snHasClosedClasses(sn: NetworkStruct): Boolean {
    return sn.njobs.hasFinite()
}
/**
 * Stochastic network HasClosedClasses algorithms
 */
@Suppress("unused")
class SnhasclosedclassesAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}