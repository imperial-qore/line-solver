package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasLCFS(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.LCFS)
}
/**
 * Stochastic network HasLCFS algorithms
 */
@Suppress("unused")
class SnhaslcfsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}