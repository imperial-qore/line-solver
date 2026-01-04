package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasMultiClassFCFS(sn: NetworkStruct): Boolean {
    for ((key, value) in sn.sched) {
        if (value == SchedStrategy.FCFS) {
            var nnz_rates = 0
            for (j in 0..<sn.nclasses) {
                if (sn.rates[key.stationIdx, j] > 0) {
                    nnz_rates++
                }
            }
            if (nnz_rates > 1) {
                return true
            }
        }
    }
    return false
}
/**
 * Stochastic network HasMultiClassFCFS algorithms
 */
@Suppress("unused")
class SnhasmulticlassfcfsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}