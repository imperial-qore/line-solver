package jline.api.sn

import jline.lang.NetworkStruct

/**
 * Checks if the network is a mixed model.
 *
 * @param sn the NetworkStruct object for the queueing network model
 * @return true if the network has mixed classes, false otherwise
 */
fun snIsMixedModel(sn: NetworkStruct): Boolean {
    return snHasMixedClasses(sn)
}
/**
 * Stochastic network IsMixedModel algorithms
 */
@Suppress("unused")
class SnismixedmodelAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}