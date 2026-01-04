package jline.api.sn

import jline.lang.NetworkStruct

/**
 * Checks if the network has one or more closed classes
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return boolean
 */
fun snHasMultipleClosedClasses(sn: NetworkStruct): Boolean {
    return sn.njobs.hasMultipleFinite()
}
/**
 * Stochastic network HasMultipleClosedClasses algorithms
 */
@Suppress("unused")
class SnhasmultipleclosedclassesAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}