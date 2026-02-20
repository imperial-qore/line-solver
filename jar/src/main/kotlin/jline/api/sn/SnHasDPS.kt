package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasDPS(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.DPS)
}
/**
 * Stochastic network HasDPS algorithms
 */
@Suppress("unused")
class SnhasdpsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}