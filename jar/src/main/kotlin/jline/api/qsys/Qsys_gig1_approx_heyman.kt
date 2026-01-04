package jline.api.qsys

import jline.io.Ret.qsys

/**
 * Analyzes a G/G/1 queueing system using Heyman's approximation.
 *
 * @param lambda Arrival rate.
 * @param mu     Service rate.
 * @param ca     Coefficient of variation of the arrival process.
 * @param cs     Coefficient of variation of the service time.
 * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
 */
fun qsys_gig1_approx_heyman(lambda: Double, mu: Double, ca: Double, cs: Double): qsys {
    val rho = lambda / mu
    val W = rho / (1 - rho) / mu * (ca * ca + cs * cs) / 2 + 1.0 / mu
    val rhohat = W * lambda / (1 + W * lambda)
    return qsys(W, rhohat)
}
/**
 * Queueing system gig1 approx heyman algorithms
 */
@Suppress("unused")
class QsysGig1ApproxHeymanAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}