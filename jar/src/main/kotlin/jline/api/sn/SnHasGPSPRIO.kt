package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasGPSPRIO(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.GPSPRIO)
}
/**
 * Stochastic network HasGPSPRIO algorithms
 */
@Suppress("unused")
class SnhasgpsprioAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}