package jline.api.qsys

import jline.io.Ret.qsys

/**
 * Analyzes a G/G/k queueing system using Kingman's approximation.
 *
 * @param lambda Arrival rate.
 * @param mu     Service rate.
 * @param ca     Coefficient of variation of the arrival process.
 * @param cs     Coefficient of variation of the service time.
 * @param k      Number of servers.
 * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
 */
fun qsys_gigk_approx_kingman(lambda: Double, mu: Double, ca: Double, cs: Double, k: Int): qsys {
    qsys_mmk(lambda, mu, k)
    val W = (ca * ca + cs * cs) / 2 * (qsys.W - 1 / mu) + 1 / mu
    val rhohat = W * lambda / (1 + W * lambda)
    return qsys(W, rhohat)
}
/**
 * Queueing system gigk approx kingman algorithms
 */
@Suppress("unused")
class QsysGigkApproxKingmanAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}