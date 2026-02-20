package jline.api.qsys

import jline.io.Ret.qsys

/**
 * Calculates an upper bound on the waiting time for a G/G/1 system using Kingman's formula.
 *
 * @param lambda Arrival rate.
 * @param mu     Service rate.
 * @param ca     Coefficient of variation of the arrival process.
 * @param cs     Coefficient of variation of the service time.
 * @return qsysReturn containing upper bound on average waiting time (W) and modified utilization (rhohat).
 */
fun qsys_gig1_ubnd_kingman(lambda: Double, mu: Double, ca: Double, cs: Double): qsys {
    val rho = lambda / mu
    val W = rho / (1 - rho) * (ca * ca + cs * cs) / 2 * (1 / mu) + (1 / mu)
    val rhohat = W * lambda / (1 + W * lambda)
    return qsys(W, rhohat)
}
/**
 * Queueing system gig1 ubnd kingman algorithms
 */
@Suppress("unused")
class QsysGig1UbndKingmanAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}