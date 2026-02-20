package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasLCFSPR(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.LCFSPR)
}
/**
 * Stochastic network HasLCFSPR algorithms
 */
@Suppress("unused")
class SnhaslcfsprAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}