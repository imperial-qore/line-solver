package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy
import jline.util.matrix.Matrix

/**
 * Checks if the network has one or more stations with multiclass heterogeneous FCFS
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return boolean
 */

fun snHasMultiClassHeterFCFS(sn: NetworkStruct): Boolean {
    for (i in 0..<sn.sched.size) {
        if (sn.sched[sn.stations[i]] != SchedStrategy.FCFS) {
            continue
        }
        val row = Matrix.extractRows(sn.rates, i, i + 1, null)
        if (row.elementMax() - row.elementMin() > 0) {
            return true
        }
    }
    return false
}
/**
 * Stochastic network HasMultiClassHeterFCFS algorithms
 */
@Suppress("unused")
class SnhasmulticlassheterfcfsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}