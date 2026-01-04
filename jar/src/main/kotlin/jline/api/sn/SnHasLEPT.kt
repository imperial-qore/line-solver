package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasLEPT(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.LEPT)
}
/**
 * Stochastic network HasLEPT algorithms
 */
@Suppress("unused")
class SnhasleptAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}