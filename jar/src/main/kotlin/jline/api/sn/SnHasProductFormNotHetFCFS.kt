package jline.api.sn

import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.lang.constant.SchedStrategy

/**
 * Checks if the network satisfies product-form assumptions (does not have heterogeneous FCFS)
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return boolean
 */

fun snHasProductFormNotHetFCFS(sn: NetworkStruct): Boolean {
    var ret = true
    for (i in 0..<sn.sched.size) {
        ret =
            ret && (sn.sched[sn.stations[i]] == SchedStrategy.INF || sn.sched[sn.stations[i]] == SchedStrategy.PS || sn.sched[sn.stations[i]] == SchedStrategy.FCFS || sn.sched[sn.stations[i]] == SchedStrategy.LCFSPR || sn.sched[sn.stations[i]] == SchedStrategy.EXT)
    }
    ret = ret && !snHasPriorities(sn)
    ret = ret && !snHasForkJoin(sn)
    ret = ret && !snHasSDRouting(sn)

    for (i in 0..<sn.sched.size) {
        if (sn.sched[sn.stations[i]] == SchedStrategy.FCFS) {
            for (j in 0..<sn.scv.numCols) {
                if (java.lang.Double.isFinite(sn.scv[i, j]) && sn.scv[i, j] > 0) {
                    ret =
                        ret && (sn.scv[i, j] > 1 - GlobalConstants.FineTol) && (sn.scv[i, j] < 1 + GlobalConstants.FineTol)
                }
            }
        }
    }

    return ret
}
/**
 * Stochastic network HasProductFormNotHetFCFS algorithms
 */
@Suppress("unused")
class SnHasProductFormNotHetFCFS {
    companion object {
        // Class documentation marker for Dokka
    }
}