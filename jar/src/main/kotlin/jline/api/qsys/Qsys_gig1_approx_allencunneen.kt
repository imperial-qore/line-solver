/**
 * @file Allen-Cunneen approximation for G/G/1 queues
 * 
 * Implements the widely-used Allen-Cunneen approximation for general G/G/1 queueing
 * systems. This two-moment approximation provides excellent accuracy for most practical
 * applications and is considered one of the best general-purpose G/G/1 approximations.
 *
 * @since LINE 3.0
 */
package jline.api.qsys

import jline.io.Ret.qsys

/**
 * Analyzes a G/G/1 queueing system using the Allen-Cunneen approximation.
 *
 * @param lambda Arrival rate.
 * @param mu     Service rate.
 * @param ca     Coefficient of variation of the arrival process.
 * @param cs     Coefficient of variation of the service time.
 * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
 */
fun qsys_gig1_approx_allencunneen(lambda: Double, mu: Double, ca: Double, cs: Double): qsys {
    val rho = lambda / mu
    val W = (rho / (1 - rho)) / mu * ((cs * cs + ca * ca) / 2) + 1 / mu
    val rhohat = W * lambda / (1 + W * lambda)
    return qsys(W, rhohat)
}
/**
 * Queueing system gig1 approx allencunneen algorithms
 */
@Suppress("unused")
class QsysGig1ApproxAllencunneenAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}