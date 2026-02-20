package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasINF(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.INF)
}
/**
 * Stochastic network HasINF algorithms
 */
@Suppress("unused")
class SnhasinfAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}