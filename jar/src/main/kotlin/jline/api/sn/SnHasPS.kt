package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasPS(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.PS)
}
/**
 * Stochastic network HasPS algorithms
 */
@Suppress("unused")
class SnhaspsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}