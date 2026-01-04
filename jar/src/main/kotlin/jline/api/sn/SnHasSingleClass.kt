package jline.api.sn

import jline.lang.NetworkStruct

fun snHasSingleClass(sn: NetworkStruct): Boolean {
    return sn.nclasses == 1
}
/**
 * Stochastic network HasSingleClass algorithms
 */
@Suppress("unused")
class SnhassingleclassAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}