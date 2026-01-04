package jline.api.qsys

import jline.io.Ret.qsys
import org.apache.commons.math3.util.FastMath
import kotlin.math.pow

/**
 * Analyzes a G/G/k queueing system using an approximation method.
 *
 * @param lambda Arrival rate.
 * @param mu     Service rate.
 * @param ca     Coefficient of variation of the arrival process.
 * @param cs     Coefficient of variation of the service time.
 * @param k      Number of servers.
 * @return qsysReturn containing average waiting time (W) and modified utilization (rhohat).
 */
fun qsys_gigk_approx(lambda: Double, mu: Double, ca: Double, cs: Double, k: Int): qsys {
    val rho = lambda / (mu * k)
    var alpha = 0.0
    alpha = if (rho > 0.7) {
        (rho.pow(k.toDouble()) + rho) / 2
    } else {
        FastMath.pow(rho, (k + 1) / 2.0)
    }
    val W = (alpha / mu) * (1 / (1 - rho)) * (ca * ca + cs * cs) / (2 * k) + 1 / mu
    val rhohat = W * lambda / (1 + W * lambda)
    return qsys(W, rhohat)
}
/**
 * Queueing system gigk approx algorithms
 */
@Suppress("unused")
class QsysGigkApproxAlgo {
    companion object {
        // Class documentation marker for Dokka
    }
}