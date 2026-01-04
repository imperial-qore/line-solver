package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasHOL(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.HOL)
}
/**
 * Stochastic network HasHOL algorithms
 */
@Suppress("unused")
class SnhasholAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}