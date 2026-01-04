/**
 * @file M/M/1 queueing system analysis
 * 
 * Implements exact analytical solutions for the M/M/1 queue (Poisson arrivals, exponential
 * service times, single server). This fundamental queueing model provides closed-form
 * expressions for key performance metrics used in capacity planning and system design.
 *
 * @since LINE 3.0
 */
package jline.api.qsys

import jline.io.Ret.qsys

/**
 * Analyzes an M/M/1 queueing system.
 *
 * @param lambda Arrival rate.
 * @param mu     Service rate.
 * @return qsysReturn containing average waiting time (W) and utilization (rho).
 */
fun qsys_mm1(lambda: Double, mu: Double): qsys {
    val rho = lambda / mu
    val W = rho / (1 - rho) / lambda
    return qsys(W, rho)
}
/**
 * Queueing system mm1 algorithms
 */
@Suppress("unused")
class QsysMm1Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}