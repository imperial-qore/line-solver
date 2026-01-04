package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasPolling(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.POLLING)
}
/**
 * Stochastic network HasPolling algorithms
 */
@Suppress("unused")
class SnhaspollingAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}