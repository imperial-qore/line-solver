package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasGPS(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.GPS)
}
/**
 * Stochastic network HasGPS algorithms
 */
@Suppress("unused")
class SnhasgpsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}