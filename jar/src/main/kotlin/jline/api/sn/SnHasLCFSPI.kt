package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasLCFSPI(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.LCFSPI)
}
/**
 * Stochastic network HasLCFSPI algorithms
 */
@Suppress("unused")
class SnhaslcfspiAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}