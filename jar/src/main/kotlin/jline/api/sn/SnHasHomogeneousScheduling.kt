package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

/**
 * Checks if the network uses an identical scheduling strategy at every station
 *
 * @param sn       - NetworkStruct object for the queueing network model
 * @param strategy - Scheduling strategy
 * @return boolean
 */

fun snHasHomogeneousScheduling(sn: NetworkStruct, strategy: SchedStrategy): Boolean {
    var stratCount = 0
    for (i in 0..<sn.sched.size) {
        if (sn.sched[sn.stations[i]] == strategy) stratCount++
    }
    return stratCount == sn.sched.size
}
/**
 * Stochastic network HasHomogeneousScheduling algorithms
 */
@Suppress("unused")
class SnhashomogeneousschedulingAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}