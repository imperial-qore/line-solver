package jline.api.sn

import jline.lang.NetworkStruct

fun snHasMultiServer(sn: NetworkStruct): Boolean {
    return sn.nservers.elementMax() > 1
}
/**
 * Stochastic network HasMultiServer algorithms
 */
@Suppress("unused")
class SnhasmultiserverAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}