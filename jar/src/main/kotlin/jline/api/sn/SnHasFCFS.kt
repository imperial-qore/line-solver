/**
 * @file Stochastic network FCFS scheduling detector
 * 
 * Identifies queueing networks using First-Come-First-Served scheduling disciplines. 
 * FCFS detection is crucial for algorithm selection and performance analysis since 
 * FCFS networks often admit efficient solution techniques.
 * 
 * @since LINE 3.0
 */
package jline.api.sn

import jline.lang.NetworkStruct
import jline.lang.constant.SchedStrategy

fun snHasFCFS(sn: NetworkStruct): Boolean {
    return sn.sched.containsValue(SchedStrategy.FCFS)
}
/**
 * Stochastic network HasFCFS algorithms
 */
@Suppress("unused")
class SnhasfcfsAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}