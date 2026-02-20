package jline.api.sn

import jline.lang.NetworkStruct

/**
 * Checks if the network uses class-switching
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return boolean
 */

fun snHasClassSwitching(sn: NetworkStruct): Boolean {
    return sn.nclasses != sn.nchains
}
/**
 * Stochastic network HasClassSwitching algorithms
 */
@Suppress("unused")
class SnhasclassswitchingAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}