package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasLJF(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.LJF)
}
/**
 * Stochastic network HasLJF algorithms
 */
@Suppress("unused")
class SnhasljfAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}