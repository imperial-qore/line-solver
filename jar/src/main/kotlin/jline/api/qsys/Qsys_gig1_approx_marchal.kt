package jline.api.qsys

import jline.io.Ret.qsys

/**
 * Analyzes a G/G/1 queueing system using Marchal's approximation.
 *
 * @param lambda Arrival rate.
 * @param mu     Service rate.
 * @param ca     Coefficient of variation of the arrival process.
 * @param cs     Coefficient of variation of the service time.
 * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
 */
fun qsys_gig1_approx_marchal(lambda: Double, mu: Double, ca: Double, cs: Double): qsys {
    val rho = lambda / mu
    val Wmm1 = rho / (1 - rho)
    val W = Wmm1 * (1 + cs * cs) / 2 / mu * (ca + rho * rho * cs * cs) / (1 + rho * rho * cs * cs) + 1 / mu
    val rhohat = W * lambda / (1 + W * lambda)
    return qsys(W, rhohat)
}
/**
 * Queueing system gig1 approx marchal algorithms
 */
@Suppress("unused")
class QsysGig1ApproxMarchalAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}