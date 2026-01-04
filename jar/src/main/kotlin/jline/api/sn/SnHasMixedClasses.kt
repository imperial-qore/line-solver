package jline.api.sn

import jline.lang.NetworkStruct

/**
 * Checks if the network has both open and closed classes
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return boolean
 */
fun snHasMixedClasses(sn: NetworkStruct): Boolean {
    return sn.njobs.hasFinite() && sn.njobs.hasInfinite()
}
/**
 * Stochastic network HasMixedClasses algorithms
 */
@Suppress("unused")
class SnhasmixedclassesAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}