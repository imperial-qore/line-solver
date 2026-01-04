package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasSJF(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.SJF)
}
/**
 * Stochastic network HasSJF algorithms
 */
@Suppress("unused")
class SnhassjfAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}