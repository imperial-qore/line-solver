/**
 * @file M/G/1 queueing system analysis
 * 
 * Implements the Pollaczek-Khinchine formula for M/G/1 queues with Poisson arrivals
 * and general service time distributions. Uses the first two moments of service time
 * to compute exact performance measures via transform methods.
 *
 * @since LINE 3.0
 */
package jline.api.qsys

import jline.io.Ret.qsys

/**
 * Analyzes an M/G/1 queueing system.
 *
 * @param lambda Arrival rate.
 * @param mu     Service rate.
 * @param cs     Coefficient of variation of the service time.
 * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
 */
fun qsys_mg1(lambda: Double, mu: Double, cs: Double): qsys {
    val rho = lambda / mu
    val Q = rho + rho * rho / (2 * (1 - rho)) + lambda * lambda * cs * cs / (mu * mu) / (2 * (1 - rho))
    val W = Q / lambda
    val rhohat = Q / (1 + Q)
    return qsys(W, rhohat)
}
/**
 * Queueing system mg1 algorithms
 */
@Suppress("unused")
class QsysMg1Algo {
    companion object {
        // Class documentation marker for Dokka
    }
}