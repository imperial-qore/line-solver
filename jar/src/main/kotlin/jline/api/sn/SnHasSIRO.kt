package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasSIRO(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.SIRO)
}
/**
 * Stochastic network HasSIRO algorithms
 */
@Suppress("unused")
class SnhassiroAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}