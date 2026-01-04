package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasDPSPRIO(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.DPSPRIO)
}
/**
 * Stochastic network HasDPSPRIO algorithms
 */
@Suppress("unused")
class SnhasdpsprioAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}