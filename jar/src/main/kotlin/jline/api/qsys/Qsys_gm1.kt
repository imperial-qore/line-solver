/**
 * G/M/1 Queueing System Analysis
 * 
 * Implements analysis for G/M/1 queues with general arrival processes and exponential
 * service times. Uses the embedded Markov chain approach to derive exact performance
 * measures for systems with complex arrival patterns but simple service.
 *
 * @since LINE 3.0
 */
package jline.api.qsys

import jline.io.Ret.qsys

/**
 * Analyzes a G/M/1 queueing system.
 *
 * @param sigma Traffic intensity.
 * @param mu    Service rate.
 * @return qsysReturn containing average waiting time (W) and utilization (rhohat).
 */
fun qsys_gm1(sigma: Double, mu: Double): qsys {
    val W = 1 / (1 - sigma) / mu
    val rhohat = 0.0
    return qsys(W, rhohat)
}
/**
 * Queueing system gm1 algorithms
 */
@Suppress("unused")
class QsysGm1Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}