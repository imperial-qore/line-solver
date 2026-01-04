package jline.api.sn

import jline.lang.NetworkStruct

fun snHasMultiChain(sn: NetworkStruct): Boolean {
    return sn.nchains > 1
}
/**
 * Stochastic network HasMultiChain algorithms
 */
@Suppress("unused")
class SnhasmultichainAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}