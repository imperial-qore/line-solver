package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun sn_has_srpt(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.SRPT)
}

/**
 * Stochastic network has SRPT scheduling check
 */
@Suppress("unused")
class SnHasSrptAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}
