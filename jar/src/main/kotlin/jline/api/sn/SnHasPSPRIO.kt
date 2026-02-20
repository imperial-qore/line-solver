package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasPSPRIO(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.PSPRIO)
}
/**
 * Stochastic network HasPSPRIO algorithms
 */
@Suppress("unused")
class SnhaspsprioAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}