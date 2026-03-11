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
        // Compute max and min ignoring NaN values to match MATLAB's max/min behavior
        var rateMax = Double.NEGATIVE_INFINITY
        var rateMin = Double.POSITIVE_INFINITY
        var hasRate = false
        for (j in 0..<row.numElements) {
            val v = row.get(j)
            if (!java.lang.Double.isNaN(v)) {
                if (v > rateMax) rateMax = v
                if (v < rateMin) rateMin = v
                hasRate = true
            }
        }
        if (hasRate && rateMax - rateMin > 0) {
            val scvs = Matrix.extractRows(sn.scv, i, i + 1, null)
            // Compute max and min of SCV ignoring NaN, matching MATLAB:
            // max(scvs) < 1+FineTol && min(scvs) > 1-FineTol
            var scvMax = Double.NEGATIVE_INFINITY
            var scvMin = Double.POSITIVE_INFINITY
            for (j in 0..<scvs.numElements) {
                val v = scvs.get(j)
                if (!java.lang.Double.isNaN(v)) {
                    if (v > scvMax) scvMax = v
                    if (v < scvMin) scvMin = v
                }
            }
            if (scvMax < 1 + GlobalConstants.FineTol && scvMin > 1 - GlobalConstants.FineTol) {
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