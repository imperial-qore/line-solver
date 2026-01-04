package jline.api.sn

import jline.lang.NetworkStruct
import jline.GlobalConstants
import jline.lang.constant.SchedStrategy
import jline.util.matrix.Matrix

/**
 * Checks if the network has one or more stations with multiclass heterogeneous FCFS
 * and exponential service
 *
 * @param sn - NetworkStruct object for the queueing network model
 * @return boolean
 */
fun snHasMultiClassHeterExpFCFS(sn: NetworkStruct): Boolean {
    for (i in 0..<sn.sched.size) {
        if (sn.sched[sn.stations[i]] != SchedStrategy.FCFS) {
            continue
        }
        val row = Matrix.extractRows(sn.rates, i, i + 1, null)
        if (row.elementMax() - row.elementMin() > 0) {
            val scvs = Matrix.extractRows(sn.scv, i, i + 1, null)
            if ((scvs.elementMax() < 1 + GlobalConstants.FineTol) && (scvs.elementMin() < 1 - GlobalConstants.FineTol)) {
                return true
            }
        }
    }
    return false
}
/**
 * Stochastic network HasMultiClassHeterExpFCFS algorithms
 */
@Suppress("unused")
class SnhasmulticlassheterexpfcfsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}