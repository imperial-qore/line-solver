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
        // Compute max and min ignoring NaN values to match MATLAB's max/min behavior.
        // EJML's elementMax/elementMin propagate NaN when the first element is NaN,
        // which incorrectly masks heterogeneous service rates in LN layer models
        // where non-visiting classes have NaN (Disabled) rates.
        var max = Double.NEGATIVE_INFINITY
        var min = Double.POSITIVE_INFINITY
        var hasValue = false
        for (j in 0..<row.numElements) {
            val v = row.get(j)
            if (!java.lang.Double.isNaN(v)) {
                if (v > max) max = v
                if (v < min) min = v
                hasValue = true
            }
        }
        if (hasValue && max - min > 0) {
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