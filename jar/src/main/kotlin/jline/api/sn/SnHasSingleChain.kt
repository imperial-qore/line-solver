package jline.api.sn

import jline.lang.NetworkStruct

fun snHasSingleChain(sn: NetworkStruct): Boolean {
    return sn.nchains == 1
}
/**
 * Stochastic network HasSingleChain algorithms
 */
@Suppress("unused")
class SnhassinglechainAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}