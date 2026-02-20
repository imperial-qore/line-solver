package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasSEPT(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.SEPT)
}
/**
 * Stochastic network HasSEPT algorithms
 */
@Suppress("unused")
class SnhasseptAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}